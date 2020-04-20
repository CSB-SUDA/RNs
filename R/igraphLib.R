
lvcent <- function(graph){
    k <- degree(graph)
    n <- vcount(graph)
    return(sapply(1:n, function(v) { mean((k[v]-k[neighbors(graph,v)]) / (k[v]+k[neighbors(graph,v)])) }))
}
