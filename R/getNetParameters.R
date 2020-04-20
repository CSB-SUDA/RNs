#' @import igraph
#' @export getNetParameters

getNetParameters <- function(x,directed=FALSE,weights=NULL,shortest.paths=FALSE){
  options(warn =-1)
  g <- graph_from_data_frame(x,directed = directed)
  result.net.vertices <- as_data_frame(g,what = "both")$vertices$name
  result.net.edges <- as_data_frame(g,what = "both")$edges
  if(length(weights)==nrow(result.net.edges) && is.numeric(weight))
    {
        weight.tmp <-weights
    }else if(is.null(weights)){
        weight.tmp <- rep(1,nrow(result.net.edges))
    }else stop("weights input is not numeric or not equal the rows of edges")
  rm(weights)
  result.degree <- degree(g)
  result.clossness.norm <- closeness(g,normalized = T)
  result.betweenness.norm <- betweenness(g,normalized = T)
  result.clossness <- closeness(g,normalized = F)
  result.betweenness <- betweenness(g,normalized = F)
  result.degree.centralization <- (sum(max(degree(g))-degree(g)))/((vcount(g)-1)*(vcount(g)-2))
  result.clossness.centralization <- (2*vcount(g)-3)*(sum(max(closeness(g))-closeness(g)))/((vcount(g)-1)*(vcount(g)-2))
  result.betweenss.centralization <- 2*betweenness(g)/((vcount(g)-1)*(vcount(g)-2))
  if(shortest.paths == TRUE){
    result.shortest.paths.mat <- matrix(rep(NA,length(result.net.vertices)^2),ncol=length(result.net.vertices),nrow=length(result.net.vertices))
    for(rowline in seq(1,length(result.net.vertices))){
      for(colline in seq(1,length(result.net.vertices))){
        tempBox <-length(shortest_paths(g,from=result.net.vertices[rowline],to=result.net.vertices[colline],weights = weight.tmp,output = "both")$vpath[[1]])
        if(tempBox==1)tempBox=NA else tempBox = tempBox-1
        result.shortest.paths.mat[rowline,colline]<-tempBox
        cat(rowline)
      }
    }
  }else{
    result.shortest.paths.mat<-matrix(NA)
  }
  result.Leverage.centrality <- lvcent(g)
  result <- new("netParameters",
                net.vertices = result.net.vertices,
                net.edges = result.net.edges,
                degree = result.degree,
                clossness.norm = result.clossness.norm,
                betweenness.norm = result.betweenness.norm,
                clossness = result.clossness,
                betweenness = result.betweenness,
                degree.centralization = result.degree.centralization,
                clossness.centralization = result.clossness.centralization,
                betweenss.centralization = result.betweenss.centralization,
                shortest.paths.mat = result.shortest.paths.mat,
                Leverage.centrality = result.Leverage.centrality)
  options(warn =0)
  return(result)
}
