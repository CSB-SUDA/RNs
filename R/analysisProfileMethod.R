#' @export analysisIt
analysisIt <- function(data = "data.frame",
                       p.value.cutoff = 0.05,
                       fold.change.cutoff = 1.5,
                       sample.label = "",
                       p.adjust.method = "BH"){
    if(!p.adjust.method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"))
        stop("rror in the p - value adjustment method")
    itProfile <- new("profileData",
                   data = data,
                   p.value.cutoff = p.value.cutoff,
                   fold.change.cutoff = fold.change.cutoff,
                   sample.label = sample.label,
                   p.adjust.method = p.adjust.method)

    return(analysisIt_intenal(itProfile))
}
