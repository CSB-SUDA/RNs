###############################
library(igraph)
library(dnet)
options(scipen=200)
nodes <- read.table("nodes.txt")
links <- read.table("links.txt")
net <- graph_from_data_frame(d = links,vertices = nodes,directed = F)
Seeds<-read.table("seeds.txt")
Seeds[1:12,1]=1
rownames(Seeds) <- nodes[,1]
PTmatrix_all<- dRWR(g=net, normalise="laplacian", setSeeds=Seeds, restart=0.8, parallel=FALSE) 
result=matrix(nrow = 190, ncol = 1)
result[,1]<-PTmatrix_all[,1]
colnames(result)="all_Affinity"
rownames(result)=nodes[,1]
result<-as.data.frame(result)
library(dplyr) # install “assertthat” first 
result2 <- mutate(result,rn=row.names(result))
sort<- arrange(result2,desc(all_Affinity))
rownames(sort)<- sort$rn
sort <- select(sort,-rn)
write.table(sort,"allseeds_0.8_sort.txt",sep="\t",quote=F,col.names =NA)
id2<-as.data.frame(sort[1:30,1])
colnames(id2)="V1"
a<-merge(id2,links,all.x=TRUE,sort=TRUE)
links2=links
links2[,1]<-links[,2]
links2[,2]<-links[,1]
b<-merge(id2,links2,all.x=TRUE,sort=TRUE)
b2=b
b[,1]<-b2[,2]
b[,2]<-b2[,1]
c<-rbind(a,b)
c2<-c[complete.cases(c),]
c3<-c2[!duplicated(c2, fromLast=TRUE), ]
write.table(c3,"top10%net.txt",quote=F,sep="\t",row.names=F,col.names=F)
bet<-edge_betweenness(net, e = E(net), directed = F) 
bet2<-cbind(links,bet)
colnames(c3)=("name1","name2")
colnames(bet2)=("name1","name2","bet")
subbet<-merge(bet2,c3,by=c("name1","name2"),sort=F)
write.table(bet2,"edge_betweenness.txt",sep="\t",quote=F,row.names=F,col.names=F)
library(GOSemSim)
hsGO2 <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont="BP", computeIC=FALSE)
aa<-mgeneSim(nodes[,1], semData=hsGO2, measure="Wang", combine="BMA", verbose=FALSE)
ut <- upper.tri(aa) 
test<-data.frame(row = rownames(aa)[row(aa)[ut]],column = rownames(aa)[col(aa)[ut]], cor =(aa)[ut] )  
test2=test
test2[,1]=test[,2]
test2[,2]=test[,1]
bind<-rbind(test,test2)
GOnet<-merge(bind,c3,by=c("name1","name2"),sort=F)  
data<-merge(GOnet,subbet,by=c("name1","name2"),sort=F) 
k1 = 99/(max(data$GO)-min(data$GO))
data[,5]=1+k1*(data$GO-min(data$GO))
colnames(data)[5]="GO1_100"
k2 = 99/(max(data$bet)-min(data$bet))
data[,6]=1+k2*(data$bet-min(data$bet))
colnames(data)[6]="bet1_100"
data[,7]=data$GO1_100+data$bet1_100
colnames(data)[7]="sumgobet"
g = graph_from_data_frame(data, directed=FALSE)
E(g)$weight = data$sumgobet
GN<-edge.betweenness.community(g,weights = E(g)$weight)
RW<-walktrap.community(g,weights = E(g)$weight)
LP<-label.propagation.community(g,weights=E(g)$weight)
SG<-spinglass.community(g,weights=E(g)$weight)
write.table(GN$membership,"top10%GN.txt",sep="\t",quote=F,col.names="community",row.names=get.vertex.attribute(g)[[1]])
write.table(RW$membership,"top10%RW.txt",sep="\t",quote=F,col.names="community",row.names=get.vertex.attribute(g)[[1]]) 
write.table(LP$membership,"top10%LP.txt",sep="\t",quote=F,col.names="community",row.names=get.vertex.attribute(g)[[1]])
write.table(SG$membership,"top10%SG.txt",sep="\t",quote=F,col.names="community",row.names=get.vertex.attribute(g)[[1]])
module<-read.table("top10%net.txt",header=F) 
m1<-read.table("top10%GN.txt") 
m1<-m1[-1,]
m11<-m1[m1$V2==1,]   
m1<-as.data.frame(m11[,1])
colnames(m1)="V1"
aa<-merge(m1,module,all.x=TRUE,sort=TRUE) 
module2=module
module2[,1]<-module[,2]
module2[,2]<-module[,1]
bb<-merge(m1,module2,all.x=TRUE,sort=TRUE)
bb2=bb
bb[,1]<-bb2[,2]
bb[,2]<-bb2[,1]
cc<-rbind(aa,bb)
cc2<-cc[complete.cases(cc),] 
cc3<-cc2[!duplicated(cc2, fromLast=TRUE), ]
write.table(cc3,"top10%GN_m1.txt",sep="\t",quote=F,row.names=F,col.names=c("name1","name2"))
###############################