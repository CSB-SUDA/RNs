setGeneric("analysisIt_intenal",function(object)standardGeneric("analysisIt_intenal"))
setMethod("analysisIt_intenal",signature(object="profileData"),function(object){
  gseMat <- object@data
  sampleCategory <- object@sample.label
  tumor<-data.frame(gseMat[,sampleCategory=="1"])
  print(data.frame(cbind(colnames(tumor),"The control group")))
  notumor<-data.frame(gseMat[,sampleCategory=="0"])
  print(data.frame(cbind(colnames(notumor),"The experimental group")))
  matFC<-NULL
  matP<-NULL
  matP.adjust<-NULL
  FCupOrDown<-NULL
  for(i in seq(1,nrow(gseMat))){
    p<-t.test(notumor[i,],tumor[i,])
    matP<-c(matP,p$p.value)
    matFC<-c(matFC,2^(sum(tumor[i,])/ncol(tumor[i,])-sum(notumor[i,])/ncol(notumor[i,])))
    if(matFC[i]>1) FCupOrDown=c(FCupOrDown,"up") else FCupOrDown=c(FCupOrDown,"down")
    if(i%%round(nrow(gseMat)/100)==0){
      cat(paste0(round(i/round(nrow(gseMat)/100)),"% expression profile analysis has been completed, NO.",i," line!!!!\n"))
    }
    if(i==nrow(gseMat)){
      cat(paste0("100% expression profile analysis has been completed，NO.",i," line!!!!\n"))
    }
  }
  matP.adjust <- p.adjust(matP,method =  object@p.adjust.method,n=length(matP))
  FC_P <-  matP.adjust < object@p.value.cutoff & (matFC > object@fold.change.cutoff | 1/matFC >object@fold.change.cutoff)
  FCupOrDown[FC_P == FALSE]="nosig"
  return(new("profileR",
             data = object@data,
             p.value.cutoff = object@p.value.cutoff,
             fold.change.cutoff = object@fold.change.cutoff,
             p.value = matP,
             p.adjust = matP.adjust,
             fold.change = matFC,
             selected.gene = FC_P,
             up.or.down = FCupOrDown,
             sample.label = object@sample.label))
})

