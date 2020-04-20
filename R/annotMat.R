annotMat <- function(geneMat, gplPlatform) {
  data("matDb")
  geneMat <- na.omit(geneMat)
  annotation_name <- sub("\\s+", "", as.character(matDb[matDb[, 1] == gplPlatform, 3]))
  if (length(annotation_name) != 0) {
          des.annotation.db <- paste(annotation_name, ".db", sep = "")
          if (!requireNamespace(des.annotation.db)) {
          BiocManager::install(des.annotation.db)
      }
      suppressMessages(library(des.annotation.db,character.only = T))
      annotation_db_SYMBOL <- paste(annotation_name, "SYMBOL", sep = "")
      annotation_db_cmd <- paste("ids<-toTable(", annotation_db_SYMBOL, ")")
      eval(parse(text = annotation_db_cmd))
      affySetExprs <- merge(ids, geneMat, by.x = "probe_id", by.y = "ID_REF")

      exprSet <- affySetExprs[affySetExprs$probe_id %in% ids$probe_id, ]
      # check no annnotation about probes
      ids <- ids[ids$probe_id %in% exprSet$probe_id, ]
      # acoording to order of befor dataset,sorting probe_id
      ids <- ids[match(exprSet$probe_id, ids$probe_id), ]
      rownames(exprSet) <- exprSet[, 1]
      tmp <- by(exprSet[, 3:ncol(exprSet)], exprSet$symbol, function(x) rownames(x)[which(rowMeans(x) %in% mean_one(rowMeans(x)))])
      probes <- as.character(tmp)
      expeSet <- exprSet[rownames(exprSet) %in% probes, ]
      ids <- ids[ids$probe_id %in% probes, ]
      rownames(expeSet) <- ids[match(expeSet$probe_id, ids$probe_id), 2]
      expeSet <- expeSet[, -1]
      expeSet <- expeSet[, -1]
      cat("dealed sucessfully!!")
      # clear name space
      return(expeSet)
  }
}

mean_one <- function(number_array) {
    number_array_distance <- abs(number_array - mean(number_array))
    return(number_array[number_array_distance %in% min(number_array_distance)][1])
}
