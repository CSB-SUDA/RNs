#####################################
setwd("E:/Study/Master/subject/data/GEO")
library(affy)
data.raw <- ReadAffy(celfile.path = "GSE19804_RAW")
data.rma <- rma(data.raw)
exp.rma <- exprs(data.rma)
write.table(exp.rma,file="GSE19804_rma.txt",sep='\t',quote=F,col.names =NA)
#####################################
library(limma)
library('gplots')
foldChange=1
padj=0.001
exprSet<-read.table("GSE19804.txt",header=T) 
group <- read.table("19804type.txt",header=T)
design <- model.matrix(~0+factor(group[,2]))
colnames(design)=levels(factor(group[,2]))
rownames(design)=colnames(exprSet) 
cont.matrix<-makeContrasts(paste0(unique(group[,2]),collapse = "-"),levels = design) 
fit <- lmFit(exprSet,design) 
fit2 <- contrasts.fit(fit, cont.matrix) 
fit2 <- eBayes(fit2) 
tempOutput = topTable(fit2, coef=1, n=Inf,adjust="BH")
nrDEG = na.omit(tempOutput)
diff =nrDEG
diffSig = diff[(diff$adj.P.Val < padj & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]
write.table(diffSig, "diffSig19804.txt",sep="\t",quote=F)  
#####################################
setwd("E:/Study/Master/subject/data")
library(clusterProfiler)
library(org.Hs.eg.db)
target_gene_id<-read.table("geneEntrezid.txt",header=T)
display_number = c(50, 50, 50)
ego_MF <- enrichGO(OrgDb="org.Hs.eg.db",
             gene = target_gene_id$EntrezGeneID,
             pvalueCutoff = 0.05,
             ont = "MF",
             readable=TRUE)
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[1], ]
ego_CC <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_id$EntrezGeneID,
                   pvalueCutoff = 0.05,
                   ont = "CC",
                   readable=TRUE)
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
ego_BP <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_id$EntrezGeneID,
                   pvalueCutoff = 0.05,
                   ont = "BP",
                   readable=TRUE)
ego_result_BP <- na.omit(as.data.frame(ego_BP)[1:display_number[3], ])
def<-rbind(as.data.frame(ego_MF)[1:50, ],as.data.frame(ego_CC)[1:50, ],as.data.frame(ego_BP)[1:50, ])
write.table(def,"GO_top50.txt",sep="\t",quote=F,col.names =NA)  
require(DOSE)
ekk <- enrichKEGG(gene=target_gene_id$EntrezGeneID,organism="human",pvalueCutoff=0.05)
ekk<-setReadable(ekk,OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write.table(ekk@result,"KEGG.txt",row.names =F,quote=F,sep="\t") 