#' A function that USES entropy to clean data
#'
#' @description The rows of gene expression profiles with low entropy,
#' all 0 and null values were deleted.Give the gene expression
#' profile (data.frame) and return the cleaned data (data.frame).
#'
#' @param radduoData data.frame,The row name is the gene and the column name is the sample (data.frame).
#' @param filtercutoff numeric,The default value is 0.1,The entropy value less than this parameter is removed.
#' @param na.dealed logical,The default value is TRUE,Whether to remove all 0 rows.
#' @param all.zero.remove logical,The default value is TRUE,Whether to remove all NA rows.
#' @export
#' @return data.frame,The gene expression profile after removing the low entropy value
#' @author Xingyi Liu
#'
geneentropyfilter <- function(radduoData,filtercutoff=0.1,na.dealed=TRUE,all.zero.remove=TRUE){
    p <- NULL
    #I'm going to get rid of all 0's
    if (na.dealed == TRUE) radduoData <- na.omit(radduoData)
    if (all.zero.remove) radduoData <- radduoData[which(rowSums(abs(radduoData))>0),]
    #Remove genes with low entropy value and copy matlab function geneentropyfilter()
    for (rowline in seq(1,nrow(radduoData))){
        rowdata<-as.numeric(as.character(radduoData[rowline,]))
        linep <- hist(rowdata,breaks=seq(min(rowdata),max(rowdata),(max(rowdata)-min(rowdata))/ceiling(ncol(radduoData)/2)))$counts
        linep <- linep/ncol(radduoData)
        linep[linep==0]=1
        entropy = -sum(linep*log2(linep))
        p<-c(p,entropy)
        cat(paste0(rowline," line are dealed!!!\n"))
    }
    radduoData <- radduoData[p>quantile(p,filtercutoff),]
    return(radduoData)
}
