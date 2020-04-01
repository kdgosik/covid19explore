## TM Analysis Function #############################

TM.Analysis<-function(corpus){

require(tm)
require(SnowballC)
require(ggplot2)
require(wordcloud)
require(fpc)
require(qgraph)

  corpus <- tm_map(corpus, content_transformer(tolower))
  corpus <- tm_map(corpus, content_transformer(removePunctuation))
  corpus <- tm_map(corpus, content_transformer(removeNumbers))
  corpus <- tm_map(corpus, content_transformer(stripWhitespace))
  
  cat("Removing Stop Words ... \n")
  myStopwords <- c(stopwords("english"),"please","thanks","thank","pm", "mg","cc")
  corpus <- tm_map(corpus, removeWords, myStopwords)
  
  freq <- ceiling((0.15 * length(corpus)))
  cat("Calculating Term Document Matrix ... \n")
  tdm <- TermDocumentMatrix(corpus, control = list(wordLengths = c(3, Inf)))
  
  cat("Finding Frequent Terms ... \n")
  freq.terms <- findFreqTerms(tdm, lowfreq = freq)
  term.freq <- rowSums(as.matrix(tdm))
  term.freq <- subset(term.freq, term.freq >= freq)
  term.freq <- sort(term.freq, decreasing = TRUE)
  
  cat("Creating Frequency Data Frame ... \n")
  freq.df <- data.frame(term = names(term.freq), freq = term.freq)
  
  cat("Finding Term Correlations ... \n")
  term.cor <- sapply(freq.terms, function(x){findAssocs(tdm, x, 0.3)})
  
  m <- as.matrix(tdm) 
  word.freq <- sort(rowSums(m), decreasing = TRUE)
  
  cat("Clustering Documents ... \n")
  m2 <- t(m)
  colnames(m2) <- rownames(tdm)
  clust <- floor(sqrt(length(corpus)/log(length(corpus))))
  kmeansResult <- kmeans(m2, clust)
  round(kmeansResult$centers, digits = 3)
  
  grp.list <- list()
  for( i in 1 : clust ) {
    clust.name <- paste("Cluster", i, sep = "")
    s <- sort(kmeansResult$centers[i,], decreasing = TRUE)
    grp.list[[i]] <- match(names(s)[1 : (10*freq)], colnames(m2))
  }
  
  cat("Saving Clustering Results ... \n")
  file.txt <- paste(deparse(substitute(corpus)), "tm.analysis.txt", sep = ".")
  sink(file.txt)
  print(term.cor)
  for( i in 1 : clust ){
    cat(paste("cluster ",i,": ",sep=""))
    s<-sort(kmeansResult$centers[i,],decreasing=T)
    cat(names(s)[1:6],"\n")
  }
  sink()
  
  cat("Saving PDF of Analysis ... \n")
  file.pdf <- paste(deparse(substitute(corpus)), "tm.analysis.pdf", sep = ".")
  
  pdf(file.pdf,onefile = TRUE, paper = "letter")
  
  ggplot(freq.df, aes(x = term, y = freq)) + 
    geom_bar(stat = "identity") + 
    xlab("Terms") + 
    ylab("Count") + 
    coord_flip()
  
  wordcloud(words = names(word.freq), 
            freq = word.freq, 
            min.freq = freq, 
            random.order = FALSE)
  
  qgraph(cor(m2), 
         threshold = 0.3, 
         cut = 0.65, 
         layout = "spring", 
         groups = grp.list, 
         vsize = 3, 
         labels = colnames(m2))
  
  dev.off()
}
