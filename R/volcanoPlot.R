#' Draw a volcano plot
#'
#' @description Based on the results of the analysisIt function, a volcano map is drawn
#'
#' @param object profileR,The data type is getted from analysisIt function.
#' @param title character,The main title of the picture.
#' @param xlab character,Notes on the horizontal axis of the picture.
#' @param ylab character,Note the vertical axis of the drawing.
#' @param points numeric,the default is 10,The number of store labels up and down.
#'
#' @importFrom ggrepel geom_text_repel
#' @import ggplot2
#' @author Xingyi Liu


setGeneric("volcanoPlot",function(object,xlab="",ylab="",title="",points=10)standardGeneric("volcanoPlot"))
setMethod("volcanoPlot",signature(object="profileR"),function(object,xlab="",ylab="",title="",points=10){
    if (title=="") title <- paste0(" Volcano Plot (PDR<",object@p.value.cutoff,"and FC >",object@fold.change.cutoff,")")
    if (xlab == "") xlab <- "log2(fold change)"
    if (ylab == "") ylab <- "-log10(FDR)"
    matFC_P<-object@data[object@p.adjust!=0,]
    volcanoThreshold<-object@up.or.down
    windowsFonts(TNR = windowsFont("Times New Roman"))
    #testMat<-cbind(volcanoThreshold,matFC_P)
    #table(testMat[testMat$volcanoThreshold=="up",]$P_value>0.05)
    volcanoXLog2FC<-log2(object@fold.change)
    volcanoYNegalog10P<--log10(object@p.adjust)
    volcanoMat<-data.frame(cbind(geneSYMOL=as.vector(rownames(object@data)),Threshold=volcanoThreshold,
                                 XLog2FC=volcanoXLog2FC,
                                 YNegalog10P=volcanoYNegalog10P))

    volcanoMatLabel<-labelHotspot(volcanoMat,points)
    y_cutoff <- max(as.numeric(as.character(object@p.adjust[object@selected.gene])))
    volcanoplot<-ggplot(data=volcanoMat,aes(x=volcanoXLog2FC,
                                          y=volcanoYNegalog10P,colour=volcanoThreshold),family="TNR")+
    geom_point(size=2.0,na.rm = TRUE)+
    scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
    geom_vline(xintercept = c(-log2(object@fold.change.cutoff),log2(object@fold.change.cutoff)),lty=4,col="black",lwd=0.5)+
    geom_hline(yintercept=-log10(y_cutoff),lty=4,col="black",lwd=0.5)+
    labs(x= xlab ,y=ylab ,title=title)+
    xlim(min(volcanoXLog2FC)-1,max(volcanoXLog2FC)+1)+
    theme(plot.title=element_text(hjust = 0.5,vjust = 2),axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),axis.text = element_text(size=16))+
    geom_text_repel(label=volcanoMatLabel,colour="black",
                    size=4,na.rm = TRUE,nudge_x = 0.2)

    print(volcanoplot)
})


labelHotspot<-function(volcanoMat,topNUM){
    volcanoMat_up<-volcanoMat[volcanoMat$Threshold=="up",]
    volcanoMat_up_gene<-tail(volcanoMat_up[order(as.numeric(as.vector(volcanoMat_up$YNegalog10P))),],topNUM)
    volcanoMat_down<-volcanoMat[volcanoMat$Threshold=="down",]
    volcanoMat_down_gene<-tail(volcanoMat_down[order(as.numeric(as.vector(volcanoMat_down$YNegalog10P))),],topNUM)
    volcanoMat$labelGHY<-""
    topNUMgene<-c(as.vector(volcanoMat_up_gene$geneSYMOL),as.vector(volcanoMat_down_gene$geneSYMOL))
    volcanoMat$labelGHY[match(topNUMgene,volcanoMat$geneSYMOL)]<-topNUMgene
    return(volcanoMat$labelGHY)
}
