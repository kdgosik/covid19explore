library(tidyverse)
library(tidytext)
library(widyr)
library(igraph)
library(visNetwork)
library(irlba)
library(pubmed.mineR)
library(broom)
library(stringr)
## Helper Functions #####

slide_windows <- function(tbl, doc_var, window_size) {
  # each word gets a skipgram (window_size words) starting on the first
  # e.g. skipgram 1 starts on word 1, skipgram 2 starts on word 2
  
  each_total <- tbl %>% 
    group_by(!!doc_var) %>% 
    mutate(doc_total = n(),
           each_total = pmin(doc_total, window_size, na.rm = TRUE)) %>%
    pull(each_total)
  
  rle_each <- rle(each_total)
  counts <- rle_each[["lengths"]]
  counts[rle_each$values != window_size] <- 1
  
  # each word get a skipgram window, starting on the first
  # account for documents shorter than window
  id_counts <- rep(rle_each$values, counts)
  window_id <- rep(seq_along(id_counts), id_counts)
  
  
  # within each skipgram, there are window_size many offsets
  indexer <- (seq_along(rle_each[["values"]]) - 1) %>%
    map2(rle_each[["values"]] - 1,
         ~ seq.int(.x, .x + .y)) %>% 
    map2(counts, ~ rep(.x, .y)) %>%
    flatten_int() +
    window_id
  
  tbl[indexer, ] %>%
    bind_cols(data_frame(window_id)) %>%
    group_by(window_id) %>%
    filter(n_distinct(!!doc_var) == 1) %>%
    ungroup
}


nearest_synonyms <- function(df, token) {
  df %>%
    widely(~ . %*% (.[token, ]), sort = TRUE)(item1, dimension, value) %>%
    select(-item2)
}


analogy <- function(df, token1, token2, token3) {
  df %>%
    widely(~ . %*% (.[token1, ] - .[token2, ] + .[token3, ]), sort = TRUE)(item1, dimension, value) %>%
    select(-item2)
  
}


## more stop words
data("stop_words")
stop_words <- rbind(stop_words, c("article", "mine"))
stop_words <- rbind(stop_words, c("publisher", "mine"))
stop_words <- rbind(stop_words, c("report", "mine"))
stop_words <- rbind(stop_words, c("doi", "mine"))

## read experimental data ###########

crohns <- readabs("~/Downloads/pubmed_result.txt")

## creating crohns tibble
crohns_tibble <- tibble(PMID = crohns@PMID,
                        text = crohns@Abstract) %>%
  mutate(text = coalesce(text),
         text = str_replace_all(text, "&#x27;|&quot;|&#x2F;", "'"),   ## weird encoding
         text = str_replace_all(text, "<a(.*?)>", " "),               ## links 
         text = str_replace_all(text, "&gt;|&lt;|&amp;", " "),        ## html yuck
         text = str_replace_all(text, "&#[:digit:]+;", " "),          ## html yuck
         text = str_replace_all(text, "<[^>]*>", " "),                ## mmmmm, more html yuck
         text = str_replace_all(text, "'s", ""),
         text = tolower(text),
         text = str_replace_all(text, "inflammatory bowel disease", "ibd"),
         text = str_replace_all(text, "ulcerative colitis" = "uc"),
         text = str_replace_all(text, "crohn disease", "cd"))      ## accronyms


## creating pmi of words in a 8 word window
tidy_pmi <- crohns_tibble %>%
  unnest_tokens(word, text) %>%
  anti_join(stop_words) %>%
  add_count(word) %>%
  filter(n >= 20) %>%
  select(-n) %>%
  slide_windows(quo(PMID), 8) %>%
  pairwise_pmi(word, window_id)


tidy_word_vectors <- tidy_pmi %>%
  widely_svd(item1, item2, pmi, nv = 256, maxit = 1000)


tidy_word_vectors %>%
  nearest_synonyms("crohn")


tidy_word_vectors %>%
  analogy("crohn", "risk", "colitis")


tidy_word_vectors %>%
  filter(dimension <= 12) %>%
  group_by(dimension) %>%
  top_n(12, abs(value)) %>%
  ungroup %>%
  mutate(item1 = reorder(item1, value)) %>%
  group_by(dimension, item1) %>%
  arrange(desc(value)) %>%
  ungroup %>%
  mutate(item1 = factor(paste(item1, dimension, sep = "__"), 
                        levels = rev(paste(item1, dimension, sep = "__"))),
         dimension = factor(paste0("Dimension ", dimension),
                            levels = paste0("Dimension ", as.factor(1:24)))) %>%
  ggplot(aes(item1, value, fill = dimension)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~dimension, scales = "free_y", ncol = 4) +
  scale_x_discrete(labels = function(x) gsub("__.+$", "", x)) +
  coord_flip() +
  labs(x = NULL, y = "Value",
       title = "First 12 principal components of the Crohn's Disease corpus",
       subtitle = "Top words contributing to the components that explain the most variation")
