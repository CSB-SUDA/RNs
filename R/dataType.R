#' @importFrom methods  setClass setGeneric setMethod setValidity

setClass("profileR",
         slots = list(
              data = "data.frame",
              p.value.cutoff = "numeric",
              fold.change.cutoff = "numeric",
              p.value = "numeric",
              p.adjust = "numeric",
              fold.change = "numeric",
              selected.gene = "logical",
              up.or.down = "character",
              sample.label = "numeric"))

setClass("profileData",
         slots = list(
             data = "data.frame",
             p.value.cutoff = "numeric",
             fold.change.cutoff = "numeric",
             sample.label = "numeric",
             p.adjust.method = "character"),
         prototype = list(
             p.value.cutoff = 0.05,
             fold.change.cutoff = 1.5,
             p.adjust.method = "BH"))

setClass("netParameters",
         slots =list(
              net.vertices = "character",
              net.edges = "data.frame",
              degree = "numeric",
              clossness.norm = "numeric",
              betweenness.norm = "numeric",
              clossness = "numeric",
              betweenness = "numeric",
              degree.centralization = "numeric",
              clossness.centralization = "numeric",
              betweenss.centralization = "numeric",
              shortest.paths.mat = "matrix",
              Leverage.centrality = "numeric"),
          prototype = list(
              degree = NULL,
              clossness.norm = NULL,
              betweenness.norm = NULL,
              clossness = NULL,
              betweenness = NULL,
              degree.centralization = NULL,
              clossness.centralization = NULL,
              betweenss.centralization = NULL,
              shortest.paths.mat = NULL,
              Leverage.centrality = NULL)
)
