library(Rtsne)
library(impute)
ty<-read.csv("path_to/example/tsne_test.csv.csv",header = T)
tsne_out <- Rtsne(
  t(ty[,-(1:2)]),
  dims = 2,
  check_duplicates = FALSE,
  perplexity = 2301,
)