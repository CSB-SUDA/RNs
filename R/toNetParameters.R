#' @importFrom igraph graph
getNetParameters <- function(tsvfile="string_interactions.tsv"){
  networkMat <- read.table(file = tsvfile,sep = "",fill =T,header=F)
  g =  make_graph(as.matrix(networkMat[,1:2], directed=FALSE))
  black_degree = degree(g)
  black_degree = as.data.frame(black_degree)
  black_betweenness = betweenness(g,normalized = T)
  black_betweenness = as.data.frame(black_betweenness)
  black_closeness = closeness(g,normalized = T)
  black_closeness = as.data.frame(black_closeness)
  averageShortestPath = 1/black_closeness
}
