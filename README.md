# CSB_RNs
Procedure flow for drug target discovery
## 一、example codes
### 1、Download the chip data in the GEO database and annotate it
```
profilerFrame <- get.GEO.data(GEO.acession = "28735",platfrom = "GPL6244")
#profilerFrame is the data.frame data type, behavioral gene, and gene expression profile data listed in the sample.
#There are not many chip platforms that can be annotated automatically, so you can experiment and do it yourself.tabel.1
profilerFrame
```
<div align=center>
  <img width="1000" src="https://github.com/liuxingyi/bioinformatics/blob/master/pictures/data_1.png"/>
</div>

### 2、Plot the boxplot plot to observe the data outline

```
boxplot(profilerFrame,col="gray",xlab="samples",ylab="express label",main="expression quentity boxplt")
```
<div align=center>
  <img width="1000" src="https://github.com/liuxingyi/bioinformatics/blob/master/pictures/data_2.png"/>
</div>

### 3、Filter the data according to the entropy value, remove the blank rows, all 0 rows
```
profilerFramecleaned<-geneentropyfilter(profilerFrame)
```
### 4、The expression profile data were analyzed according to the classification
```
profileData <- analysisIt(profilerFrame,sample.label = c(rep(c(1,0),times=45)))
```
### 5、Draw a volcano map according to the threshold

```
volcanoPlot(profileData)
```
<div align=center>
  <img width="700" src="https://github.com/liuxingyi/bioinformatics/blob/master/pictures/data_3.png"/>
</div>

### 5、Plot the barplot of enrichment analysis

```
geneSet <- rownames(profileData@data)[profileData@selected.gene]
plotGOtermsbarplot(geneSet = geneSet)
```

<div align=center>
  <img width="700" src="https://github.com/liuxingyi/bioinformatics/blob/master/pictures/data_4.png"/>
</div>

### 6、Calculate the parameters of the protein network exported from the string database

```
net <- read.table("./string_interactions.tsv",header = T)
netParameters <-getNetParameters(net)
CSB_RNs(profileData,netParameters)
```
<div align=center>
  <img width="700" src="https://github.com/liuxingyi/bioinformatics/blob/master/pictures/data_5.png"/>
</div>


## 二、installation
```
install.packages("devtools")
library(devtools)
install_github("liuxingyi/CSB_RNs")

#load installed packages
library("targetFinder")
```

## **appendix**

### table.1 Chip platform that can be annotated

gpl |Organism | bioc_package
----|----------|---------
GPL32| "Mus musculus"| mgu74a
GPL33| "Mus musculus"| mgu74b
GPL34 |"Mus musculus"| mgu74c
GPL74 |"Homo sapiens"| hcg110
GPL75 |"Mus musculus"| mu11ksuba
GPL76 |"Mus musculus"| mu11ksubb
GPL77 |"Mus musculus"| mu19ksuba
GPL78 |"Mus musculus" |mu19ksubb
GPL79 |"Mus musculus" |mu19ksubc
GPL80 |"Homo sapiens" |hu6800
GPL81 |"Mus musculus"| mgu74av2
GPL82| "Mus musculus" |mgu74bv2
GPL83 | "Mus musculus" |mgu74cv2
GPL85 |"Rattus norvegicus"| rgu34a
GPL86 |"Rattus norvegicus" |rgu34b
GPL87 |"Rattus norvegicus" |rgu34c
GPL88 |"Rattus norvegicus" |rnu34
GPL89 |"Rattus norvegicus" |rtu34
GPL91 |"Homo sapiens" |hgu95av2
GPL92 |"Homo sapiens" |hgu95b
GPL93 |"Homo sapiens" |hgu95c
GPL94 |"Homo sapiens" |hgu95d
GPL95 |"Homo sapiens" |hgu95e
GPL96 |"Homo sapiens" |hgu133a
GPL97 |"Homo sapiens" |hgu133b
GPL98 |"Homo sapiens"| hu35ksuba
GPL99 |"Homo sapiens"| hu35ksubb
GPL100 |"Homo sapiens" |hu35ksubc
GPL101 |"Homo sapiens" |hu35ksubd
GPL201 |"Homo sapiens"| hgfocus
GPL339 |"Mus musculus" |moe430a
GPL340 |"Mus musculus" |mouse4302
GPL341 |"Rattus norvegicus" |rae230a
GPL342 |"Rattus norvegicus" |rae230b
GPL570 |"Homo sapiens" |hgu133plus2
GPL571 |"Homo sapiens" |hgu133a2
GPL886 |"Homo sapiens" |hgug4111a
GPL887 |"Homo sapiens" |hgug4110b
GPL1261 |"Mus musculus" |mouse430a2
GPL1352 |"Homo sapiens" |u133x3p
GPL1355 |"Rattus norvegicus" |rat2302
GPL1708 |"Homo sapiens" |hgug4112a
GPL2891 |"Homo sapiens" |h20kcod
GPL2898 |"Rattus norvegicus" |adme16cod
GPL3921 |"Homo sapiens"| hthgu133a
GPL4191 |"Homo sapiens" |h10kcod
GPL5689 |"Homo sapiens" |hgug4100a
GPL6097 |"Homo sapiens" |illuminaHumanv1
GPL6102 |"Homo sapiens" |illuminaHumanv2
GPL6244 |"Homo sapiens" |hugene10sttranscriptcluster
GPL6947 |"Homo sapiens" |illuminaHumanv3
GPL8300 |"Homo sapiens" |hgu95av2
GPL8490 |"Homo sapiens" |IlluminaHumanMethylation27k
GPL10558| "Homo sapiens" |illuminaHumanv4
GPL11532 |"Homo sapiens" |hugene11sttranscriptcluster
GPL13497| "Homo sapiens"| HsAgilentDesign026652
GPL13534 |"Homo sapiens" |IlluminaHumanMethylation450k
GPL13667 |"Homo sapiens" |hgu219
GPL15380 |"Homo sapiens" |GGHumanMethCancerPanelv1
GPL15396 |"Homo sapiens" |hthgu133b
GPL17897 |"Homo sapiens" |hthgu133a


