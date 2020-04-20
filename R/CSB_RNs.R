#' Compute the function of CSB_RNs
#'
#' @description The CSB_RNs parameters were calculated according to
#' the expression spectrum and network topology parameters.
#'
#' @param data profileR,The data type provided in the package
#' @param netParameterCsvfile The exported CSV network file from
#' the STRING database is entered into cytoscape, and the exported
#' file is calculated by the plug-in networkAnalysis
#' @return data.frame Returns the calculated result data
#' @importFrom e1071 svm
#' @export CSB_RNs
#' @author Xingyi Liu
#' @seealso \code{\link{profileR}}


CSB_RNs <- function(data,netParameterMatrix){
    if(!class(netParameterMatrix)[1] == "netParameters")stop("the type of the parameter--netParameterMatrix is Error!!")
    if(!class(data)[1] == "profileR")stop("the type of the data--profileR is Error!!")
    pTable <-data@data
    pTable <- pTable[rownames(pTable)%in% netParameterMatrix@net.vertices,]
    pTable<-t(pTable)
    featureRankedList = svmrfeFeatureRanking(pTable,data@sample.label)
    geneOrder<-data.frame(cbind(rank=featureRankedList,ID=colnames(pTable)))
    geneOrder<-geneOrder[order(featureRankedList),]
    geneOrder<-cbind(geneOrder,score=round(((1+length(geneOrder[,1]))-as.numeric(as.vector(geneOrder$rank)))/length(geneOrder[,1]),5))
    other_score<-data.frame(AverageShortestPathLength=1/netParameterMatrix@clossness.norm,Degree=netParameterMatrix@degree,name=netParameterMatrix@net.vertices)
    lastScore<-merge(other_score,geneOrder[,c(2,3)],by.x="name",by.y="ID",all.x=F,all.y=F)
    resultmat<-data.frame(lastScore,CSB_RN=lastScore[,3]*lastScore[,4]/lastScore[,2])
    resultmat<-resultmat[order(resultmat$CSB_RN,decreasing = T),]
    return(resultmat)
}
