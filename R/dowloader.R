##' dowload gene profiler from ncbi ftp website and return data.frame about that
##'
##' @name get.GEO.data
##' @docType function
##' @rdname get.GEO.data-function
##' @title get.GEO.data function
##' @param GEO.acession GEO acession id to dowload file
##' @param raw.data if TRUE,dowload raw data
##' @param result.data gene profiler data
##' @export get.GEO.data
##' @aliases get.GEO.data dowloadGEO
##' @author Xingyi Liu

get.GEO.data <- function(GEO.acession="28735",platfrom="GPL6244"){
    GSEGetted <- paste0("GSE",GEO.acession)
    #from a url based on input
    if(nchar(GSEGetted)==9){
        dirNnn<-paste(substring(GSEGetted,1,6),"nnn/",sep="")
    }else{
        dirNnn<-paste(substring(GSEGetted,1,5),"nnn/",sep="")
    }
    prefixUrl<-"ftp://ftp.ncbi.nlm.nih.gov/geo/series/"
    matName<-paste(GSEGetted,"_series_matrix.txt.gz",sep="")
    url<-paste(prefixUrl,dirNnn,GSEGetted,"/matrix/",matName,sep="")
    cat("Start downloading\n")
    download.file(url,destfile=paste0("./tmp/",matName),quiet = TRUE)
    cat("download successfully\n")
    data <- read.table(paste0("./tmp/",matName),
                       sep = "\t",
                       header=T,
                       fill= TRUE,
                       stringsAsFactors = FALSE,
                       comment.char = "!")
    unlink(paste0("./tmp/",matName))
    data <- annotMat(data,platfrom)
    return(data)
}

