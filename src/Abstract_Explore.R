library(pubmed.mineR)
library(rentrez)
library(magrittr)
library(igraph)
library(dplyr)
library(purrr)
library(networkD3)
library(foreach)
library(doParallel)
library(stringdist)
library(ggplot2)
rm(list = ls()); gc(reset = TRUE)

## helper functions ####

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

suggest_word <- function( word, list ) {
  
  list[which.min(stringdist(
    stemDocument(word), 
    stemDocument(list), method = "lv"))]
  
}

remove_dups <- function( word_list ){
  i <- 1
  while( i <= length(word_list) ) {
    
    remove_idx <- -which(stemDocument(tolower(word_list[i])) == stemDocument(tolower(word_list)))[-1]
    
    if( length(remove_idx) > 0 ) word_list <- word_list[remove_idx]
    i <- i + 1
  }
    
  word_list
}


doid <- read.csv("data/DOID.csv", stringsAsFactors = FALSE)

# first 400 abstracts
# abstracts <- readabs("pubmed_result_crohns_search1.txt")
abstracts <- readabs("data/pubmed_result_colorectalcancer.txt")
pmids <- abstracts@PMID[1:200]


cl <- makeCluster(4)
registerDoParallel(cl)

pubmed_list <- foreach( i = pmids, .export = "pubtator_function" ) %dopar% {
  pubtator_function(i)
}

stopCluster(cl)

disease_data <- pubmed_list %>%
  transpose %$%
  #Map(function(a, b) data.frame(cbind(id = a, disease = remove_dups(b))), PMID, Diseases) %>%
  Map(function(a, b) data.frame(cbind(id = a, disease = b)), PMID, Diseases) %>%
  Filter(function(x) NROW(x) != 0, .) %>%
  Filter(function(x) NCOL(x) == 2, .) %>%
  do.call(rbind, .) %>%
  mutate(disease = tolower(disease), 
         disease = gsub("diseases", "disease", disease)) %>%
  unique

# keeping diseases over 3 letters long, (non-abbreviated)
disease_data <- disease_data[nchar(disease_data$disease) > 3, ]


genetic_data <- pubmed_list %>%
  transpose %$%
  #Map(function(a,b) data.frame(cbind(id = a, gene = remove_dups(b))), PMID, Genes) %>%
  Map(function(a, b) data.frame(cbind(id = a, gene = b)), PMID, Diseases) %>%
  Filter(function(x) NROW(x) != 0, .) %>%
  Filter(function(x) NCOL(x) == 2, .) %>%
  do.call(rbind, .) %>%
  mutate(gene = as.character(gene)) %>%
  unique

data("HGNCdata")
HGNCdata %<>% 
  mutate(Approved.Symbol = as.character(Approved.Symbol),
         Approved.Name = as.character(Approved.Name),
         Previous.Symbols = as.character(Previous.Symbols),
         Previous.Names = as.character(Previous.Names),
         Synonyms = as.character(Synonyms))


for( i in 1 : nrow(HGNCdata) ) {

  name_list <- as.character(HGNCdata[i, c("Approved.Name", "Previous.Symbols", "Previous.Names", "Synonyms")])
  genetic_data$gene[tolower(genetic_data$gene) %in% tolower(name_list)] <- paste(HGNCdata[i, "Approved.Symbol"])

}

for( i in 1 : nrow(HGNCdata) ) {

  genetic_data$gene[tolower(genetic_data$gene) %in% tolower(genetic_data$gene)[i]] <- genetic_data$gene[i]

}

tnf_list <- c("TNF", "TNFa", "TNF-a", "Tumor necrosis factor (TNF)-alpha",
              "TNF-alpha", "TNF alpha", "TNFA", "tumor necrosis factor",
              "tumor necrosis factor-a", "tumor necrosis factor alpha",  
              "tumor necrosis factor-alpha", "tumor necrosis factor a",
              "Tumor Necrosis Factor Alpha", "tumour necrosis factor",
              "tumour necrosis factor-a")
genetic_data$gene[genetic_data$gene %in% tnf_list] <- "TNF"

anti_tumor <- c("anti-tumour necrosis factor a", "anti-tumour necrosis factor",
                "anti-tumor necrosis factor", "Anti-tumor necrosis factor-a",
                "anti-tumour necrosis factor")
genetic_data$gene[genetic_data$gene %in% anti_tumor] <- "anti-tumour necrosis factor"

genetic_data$gene[genetic_data$gene %in% "interleukin-1"] <- "IL-1"
genetic_data$gene[genetic_data$gene %in% "interleukin (IL)-6"] <- "IL-6"
genetic_data$gene[genetic_data$gene %in% c("interleukin-10", "Interleukin-10")] <- "IL-10"
genetic_data$gene[genetic_data$gene %in% "interleukin 17"] <- "IL-17"
genetic_data$gene[genetic_data$gene %in% "Interleukin-18"] <- "IL-18"
genetic_data$gene[genetic_data$gene %in% c("interleukin-23", "interleukin 23", "IL23", "Interleukin 23")] <- "IL-23"
genetic_data$gene[genetic_data$gene %in% c("interleukin 28A", "IL28A", "IL-28a")] <- "IL-28a"
genetic_data$gene[genetic_data$gene %in% "Interleukin (IL)-32"] <- "IL-32"
genetic_data$gene[genetic_data$gene %in% "TLR)2"] <- "TLR2"
genetic_data$gene[genetic_data$gene %in% "VDR"] <- "vitamin D receptor"
genetic_data$gene[genetic_data$gene %in% c("stat1", "Stat1IEC-KO")] <- "STAT1"
genetic_data$gene[genetic_data$gene %in% "signal transducer and activator of transcription 3"] <- "STAT3"
genetic_data$gene[genetic_data$gene %in% c("stat1", "Stat1IEC-KO")] <- "STAT1"
genetic_data$gene[genetic_data$gene %in% "P-gp"] <- "P-glycoprotein" 
genetic_data$gene[genetic_data$gene %in% "CRP  <  3"] <- "CRP"

genetic_data <- unique(genetic_data)


chemical_data <- pubmed_list %>%
  transpose %$%
  #Map(function(a,b) data.frame(cbind(id = a, chemicals = remove_dups(b))), PMID, Chemicals) %>%
  Map(function(a,b) data.frame(cbind(id = a, chemicals = b)), PMID, Chemicals) %>%
  Filter(function(x) NROW(x) != 0, .) %>%
  Filter(function(x) NCOL(x) == 2, .) %>%
  do.call(rbind, .) %>%
  mutate(chemicals = tolower(chemicals)) %>%
  unique



# net <- graph_from_data_frame(genetic_data)
# plot(net)
simpleNetwork(disease_data, zoom = TRUE)
simpleNetwork(genetic_data, zoom = TRUE)
simpleNetwork(chemical_data, zoom = TRUE)


disease_gene <- full_join(disease_data, genetic_data) %>%
  select(-id) %>%
  na.omit

disease_gene %>% 
  group_by(disease, gene) %>% 
  summarise(N = n()) %>% 
  arrange(-N)

disease_chemical <- full_join(disease_data, chemical_data) %>%
  select(-id) %>%
  na.omit

chemical_gene <- full_join(chemical_data, genetic_data) %>%
  select(-id) %>%
  na.omit

### non - weighted ######
# tmp <- disease_gene %>%
#   select(from = disease, to = gene) %>%
#   rbind(. , {disease_chemical %>% select(from = disease, to = chemicals)}) %>%
#   rbind(. , {chemical_gene %>% select(from = chemicals, to = gene)})


### weighted ############
tmp <- disease_gene %>% 
  group_by(disease, gene) %>% 
  summarise(N = n()) %>% 
  arrange(-N) %>%
  select(from = disease, to = gene, weight = N) %>%
  rbind(. , {disease_chemical %>% 
      group_by(disease, chemicals) %>% 
      summarise(N = n()) %>% 
      arrange(-N) %>%
      select(from = disease, to = chemicals, weight = N)
    
  }) %>%
  rbind(. , {chemical_gene %>% 
      group_by(chemicals, gene) %>% 
      summarise(N = n()) %>% 
      arrange(-N) %>%
      select(from = chemicals, to = gene, weight = N)
    
  }) %>%
  mutate( weight = 1/weight)

nodes <- do.call(rbind,
                 list(
                   data.frame(
                     label = unique(disease_data$disease)[unique(disease_data$disease) %in% unique(c(tmp$from, tmp$to))], 
                     group = "Disease",
                     stringsAsFactors = FALSE),
                   data.frame(
                     label = unique(genetic_data$gene)[unique(genetic_data$gene) %in% unique(c(tmp$from, tmp$to))], 
                     group = "Gene",
                     stringsAsFactors = FALSE),
                   data.frame(
                     label = unique(chemical_data$chemicals)[unique(chemical_data$chemicals) %in% unique(c(tmp$from, tmp$to))], 
                     group = "Chemical",
                     stringsAsFactors = TRUE) 
                 )
) %>%
  mutate(id = 1 : n() - 1,
         label = factor(label),
         group = factor(group)) %>%
  select(id, everything())

# nodes <- do.call(rbind,
#                  list(
#                    data.frame(
#                      label = unique(disease_data$disease), 
#                      group = "Disease",
#                      stringsAsFactors = FALSE),
#                    data.frame(
#                      label = unique(genetic_data$gene), 
#                      group = "Gene",
#                      stringsAsFactors = FALSE),
#                    data.frame(
#                      label = unique(chemical_data$chemicals), 
#                      group = "Chemical",
#                      stringsAsFactors = FALSE) 
#                  )
#                  ) %>%
#   mutate(id = 1 : n() - 1) %>%
#   select(id, everything())

edges <- data.frame(
               source = match(tmp$from, nodes$label) - 1,
               target = match(tmp$to, nodes$label) - 1,
               value = tmp$weight,
               stringsAsFactors = FALSE
               ) %>%
  unique
  
forceNetwork(Links = edges, Nodes = nodes,
             Source = "source", Target = "target", Value = "value", 
             NodeID = "label", Group = "group", 
             linkDistance = JS("function(d){return d.value * 10}"),
             opacity = 0.4, zoom = TRUE)