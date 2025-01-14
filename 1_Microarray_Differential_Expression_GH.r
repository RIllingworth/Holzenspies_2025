#This R script reads in raw microarray quantitation files, normalises them and analyses for differntial gene expression based on the specificed paramaters
##the raw microarray files can be sourced from GEO accession - GSE286195 


##script start
options(scipen=20)
rootdir<-"./top/level/path/to/analysis/location"

#load the required R libraries
library(limma)
library(Biobase)

#set working directory for the analysis
setwd(paste(rootdir,"/Analysis",sep=""))

##read sample key including sample identification data to read in the raw files
target<-read.delim("../Target_Files/TargetAll_Marray.txt",header=TRUE)

##Ceeate an object containing all your expression data and annotation - recogniseable  
Expr<-read.maimages(files=as.character(unlist(target$File)), source="agilent", names=unlist(target[,2]),verbose=TRUE,columns=list(E="gMedianSignal",Eb="gBGMedianSignal"), other.columns=NULL, annotation=list(Control="ControlType",Name="ProbeName",ID="SystematicName",FeatureNum="FeatureNum"), green.only=TRUE,sep="\t")
colnames(Expr$E)<-colnames(Expr$Eb)<-target[,2]

############################################################################################
#Remove cotrol probes and collapse replicates based on gene name (symbol) annotated using the supplier's (Agilent) GAL file. annotation can also be perormed using the GEO platform associated with the GEO series accession for this study (GSE286195) 
#GAL (annotation files) is ordered based on array format, not as per the analysis output so must be reorderred accordingly

Ano <- readGAL("../Annotation_Files/028005_D_20130207.gal")
Ano<-Ano[which(Ano$Block==1),]
a<-order(Ano$RefNumber)
Expr$genes<-cbind(Expr$genes,Ano[a,])
a<-which(Expr$genes$ControlType == 0)
Expr<-Expr[a,]
   

###save out a normalised expression matrix for GEO
####Normalise the data including a nominal offset value of 10 (normalisation converts values to a log2 scale)
GEO_BS<-backgroundCorrect(Expr, method="normexp", offset=10)
GEO_Norm <-normalizeBetweenArrays(GEO_BS, method="quantile")
GEO<-cbind(GEO_Norm$genes[,c(2,4)],GEO_Norm$E)
GEO<-GEO[,c(1:41,45:47)]
write.table(GEO,file="./top/level/path/to/analysis/location/Analysis/Arrays/GEO/NormalisedTable.txt",sep="\t",col.names=TRUE,row.names=FALSE)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#Collapse data to individual array probe values(duplicated probes are collapsed to a single mean value for both E and Eb values)
Expr$genes$Collapse<-Expr$genes$ProbeName
tmpE<-aggregate(data.frame(Expr$E),by=list(Expr$genes$Collapse), FUN=mean)
tmpEb<-aggregate(data.frame(Expr$Eb),by=list(Expr$genes$Collapse), FUN=mean)
a<-which(duplicated(Expr$genes$Collapse))
b<-1:length(Expr$genes$Collapse)
c<-setdiff(b,a)

Expr<-Expr[c,]
a<-match(Expr$genes$Collapse,tmpE[,1])
tmpE<-tmpE[a,1:ncol(tmpE)]
tmpEb<-tmpEb[a,1:ncol(tmpEb)]

Expr$E<-tmpE[,2:(ncol(tmpE))]
Expr$Eb<-tmpEb[,2:(ncol(tmpEb))]

#rename columns as per target file annotation
colnames(Expr)<-target[,2]

#give columns unique numeric values
row.names(Expr$E)<-c(1:nrow(Expr$E))
row.names(Expr$Eb)<-c(1:nrow(Expr$Eb))
row.names(Expr$genes)<-c(1:nrow(Expr$E))

#Remove refseq genes which are not also analysed in the later 4SU data - used a 4SU dataset to perform the trim
Subset<-read.table("./top/level/path/to/analysis/location/Annotation_Files/RefseqGenes_TrimSet.txt",header=TRUE)
a<-match(Expr$genes$SystematicName,Subset[,1])
b<-which(a>0)
Expr<-Expr[b,]

############################################################################################
####Normalise the data including a nominal offset value of 10 - normalisation also collapses values to a log2 scale
ExprBS<-backgroundCorrect(Expr, method="normexp", offset=10)
ExprNorm <-normalizeBetweenArrays(ExprBS, method="quantile")

####QC Plots
##create a boxplot of normalisation
setwd("./top/level/path/to/analysis/Analysis/Arrays/QC")
jpeg(filename = "boxplot.jpg")
par(mfcol=c(3,1))   
boxplot(log2(Expr$E),ylim=c(0,20))
boxplot(log2(ExprBS$E),ylim=c(0,20))
boxplot(ExprNorm$E,ylim=c(0,20))
dev.off()

####clustering pre and post normalisation
jpeg(filename = "Clustering.jpg",height=600,width=1800)
par(mfcol=c(2,1))   
Raw <- dist(t(as.matrix(log2(Expr$E),byrow=TRUE)))
plot(hclust(Raw), labels=colnames(Expr$E))
Norm <- dist(t(as.matrix(ExprNorm$E,byrow=TRUE)))
plot(hclust(Norm), labels=colnames(ExprNorm$E))
dev.off()

############################################################################################QC Plot candidate genes of interest across all replicates
dir<-paste(getwd(),"Genes",sep="/")
dir.create(dir)
setwd(dir)

X<-c("Sox2","Gsc","Nanog","Gata4","Gata6","Hhex","T","Pou5f1","Sox6","Klf4","Dab2","Actb","Gapdh","Rex1","Dnmt3a","Hoxa1","Hoxa9","Hoxb1","Hoxb9","Snai1","Cxcr4","Eomes","Cdx2")

for(i in 1:length(X)){
a<-grep(paste("^",X[i],"$",sep=""), ExprNorm$genes$GeneName)

if(length(a)==1){
jpeg(filename = paste(X[i],".jpeg",sep=""))
plot(ExprNorm$E[a,],main=X[i],pch=16,ylab="log2(Normalised Expression)",ylim=c(0,18))
par(new=TRUE)
plot(ExprNorm$E[a,],main="",xlab="",ylab="",axes=FALSE,ylim=c(0,18))
#legend("topright",legend=unique(target[,4]),col=c("white","grey90","grey80","grey70","grey60","grey50","grey40","grey30","grey20","grey10","black"),bty="n",cex=0.7)
dev.off()}
else {NULL}
}

setwd(paste(rootdir,"/Analysis",sep=""))

######### Create Limma 'Fit' object
a<-unique(target[,4])
#create design matrix
design<-matrix(data=0,nrow=nrow(target),ncol=length(a))
colnames(design)<-a
for(i in 1:length(a)){
tmp<-grep(unique(a[i]),target[,4])
design[tmp,i]<-1}
colnames(design)<-sub(" ","_",colnames(design))
colnames(design)<-sub(" ","_",colnames(design))
fit<-lmFit(ExprNorm,design)

##########plot combined genes from 'Fit' object
setwd(paste(rootdir,"/Analysis/Arrays/",sep=""))
   
dir<-paste(getwd(),"Genes_CombinedReps",sep="/")
dir.create(dir)
setwd(dir)

X<-c("Sox2","Gsc","Nanog","Gata4","Gata6","Hhex","T","Pou5f1","Sox6","Klf4","Dab2","Actb","Gapdh","Rex1","Dnmt3a","Hoxa1","Hoxa2","Hoxa3","Hoxa5","Hoxa7","Hoxa9","Hoxb1","Hoxb9","Snai1","Cdh1","Cdh2","Snai1","Snai2","Zeb1", "Zeb2", "Klf8","Cxcr4","Eomes","Cdx2","Shh")

for(i in 1:length(X)){
a<-grep(paste("^",X[i],"$",sep=""), fit$genes$GeneName)

if(length(a)>=1){
jpeg(filename = paste(X[i],".jpeg",sep=""),height=400,width=1000)
tmp<-colMeans(rbind(fit$coefficients[a,],fit$coefficients[a,]))
plot(tmp,main=X[i],col=c("grey50"),pch=16,axes=FALSE,ylab="log2(Normalised Expression)",ylim=c(0,18))
axis(1,at=c(1:length(tmp)),labels=c(names(tmp)),cex=0.7)
axis(2,at=c(0,9,18))
par(new=TRUE)
plot(tmp,main="",xlab="",ylab="",axes=FALSE,ylim=c(0,18))
#legend("topright",legend=unique(target[,4]),col=c("white","grey90","grey80","grey70","grey60","grey50","grey40","grey30","grey20","grey10","black"),bty="n",cex=0.7)
dev.off()}
else {NULL}
}

      

##############################Differential gene Expression Analysis# # # # # # # # # # # #
##plot MA plots of each comparison as a single .pdf files
setwd(paste(rootdir,"/Analysis/Arrays/Differential_Expression_Analysis",sep=""))
  
#read in the target file with the required contrast compariosns and output names
Contrasts<-read.table("./top/level/path/to/analysis/location/Target_Files/Contrasts_Marray.txt",header=TRUE,colClasses="character")

## ## ## thresholds for differential expression analysis ## ## ## ## 
#adjusted p-value threshold  
co<-0.01
#log2 Fold Change threshold
diff<-1.5
dir<-paste(getwd(),"TopTables",sep="/")
dir.create(dir)
setwd(dir)

Stats<-matrix(0,ncol=6,nrow=nrow(Contrasts),dimnames=list(Contrasts[,2],c("UP Probes","DOWN Probes","UP GeneNames","DOWN GeneNames","UP Refseq","DOWN Refseq")))
for(k in 1:nrow(Contrasts)){
x<-Contrasts[k,1]
name<-Contrasts[k,2]
cont.matrix <- makeContrasts(x,levels=colnames(design))
fit2 <- contrasts.fit(fit, cont.matrix)

#############limma based analysis applying the contrasts and the thresholds specificed above
fit2 <- eBayes(fit2)
TT<-topTable(fit2,adjust.method="BH",number=nrow(fit2),sort.by="none")

#####
pdf(paste(name,".pdf",sep=""),useDingbats=FALSE)
par(mar=c(5.1, 5.1, 4.1, 2.1), bty="l")
smoothScatter(TT$AveExpr,TT$logFC,xlab="Mean Expression (Log2)",ylab="Differential Expression (Log2)",main=name,cex.main=2,cex.lab=2,ylim=c(-10,10),nbin = 500,bandwidth=c(0.05,0.1),nrpoints = 0, colramp = colorRampPalette(c("white", "black")),axes=FALSE)
axis(1)
axis(2,at=c(-10,-5,0,5,10))
a<-which((TT$logFC>=diff)&(TT$adj.P.Val<=co))
b<-which((TT$logFC<=(-diff))&(TT$adj.P.Val<=co))
c<-length(unique(TT$GeneName[a]))
d<-length(unique(TT$GeneName[b]))
e<-length(unique(TT$SystematicName[a]))
f<-length(unique(TT$SystematicName[b]))
points(TT$AveExpr[a],TT$logFC[a],pch=21, bg="#d95f02")
points(TT$AveExpr[b],TT$logFC[b],pch=21, bg="#7570b3")
xl<-sub("-.*","",x)
yl<-sub(".*-","",x)
legend("topright",legend=c(paste("Up in ",xl," = ",length(a)," (",e,")",sep=""),paste("Up in ",yl," = ",length(b)," (",f,")",sep="")),bty="n",cex=0.7)
dev.off()
###############
Stats[k,1]<-length(a)
Stats[k,2]<-length(b)
Stats[k,3]<-c
Stats[k,4]<-d
Stats[k,5]<-e
Stats[k,6]<-f

##whole TT
write.table(TT,file=paste(name,"TopTable.txt",sep="_"),sep="\t",col.names=TRUE,row.names=FALSE)
##UP
write.table(TT[a,c(3,12,16)],file=paste(name,"pval",co,"diff",diff,"UP.txt",sep="_"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(TT[b,c(3,12,16)],file=paste(name,"pval",co,"diff",diff,"DOWN.txt",sep="_"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}
write.table(Stats,paste("DifferentialExpression_ManualComparisons_minlogFC=",diff,"_minp=",co,".txt",sep=""),sep="\t",quote=FALSE)

#####Make a summary table with all of the information
b<-list.files()
a<-grep("_TopTable.txt",b)
files<-b[a]
names<-sub("_TopTable.txt","",files)

TT<-read.table(files[1],header=TRUE)  
a<-order(TT[,2])
Summary<-TT[a,c(2,3,12)]

for(i in 1:length(files)){
tmp<-read.table(files[i],header=TRUE)  
a<-order(tmp[,2])
tmp<-data.frame(tmp[a,c(16,20)])
colnames(tmp)<-paste(names[i],colnames(tmp),sep="_")
Summary<-cbind(Summary,tmp)}

a<-match(fit$gene$ProbeName,Summary$ProbeName)

ExprMerged<-cbind(fit$genes[,c(2,3,12)],fit$coefficients,Summary[a,])
ExprMerged<-ExprMerged[,c(1:22,26:ncol(ExprMerged))] 
write.table(ExprMerged,"Summary_Marray_Expression_Table.txt",sep="\t",quote=FALSE)  

coord<-list()
a<-grep("adj.P.Val",colnames(ExprMerged))
for(i in 1:length(a)){
coord[[i]]<-which(ExprMerged[,a[i]]<=0.01)}
coord<-unique(unlist(coord))
write.table(ExprMerged[coord,],"Significant_Marray_Expression_Table.txt",sep="\t",quote=FALSE)  
 
#Repeat analysis and make split plots for high resolution publishable figures        
           
##############################Differential gene Expression Analysis.
##plot single pdf files
setwd(paste(rootdir,"/Analysis/Arrays/Differential_Expression_Analysis",sep=""))
#read in file with required contrasts and names
Contrasts<-read.table("../../../Target_Files/Contrasts_Marray.txt",header=TRUE,colClasses="character")

dir<-paste(getwd(),"Split_Plots",sep="/")
dir.create(dir)
setwd(dir)

######### Create 'Fit' object
for(k in 1:nrow(Contrasts)){
x<-Contrasts[k,1]
name<-Contrasts[k,2]
cont.matrix <- makeContrasts(x,levels=colnames(design))
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
TT<-topTable(fit2,adjust.method="BH",number=nrow(fit2),sort.by="none")

###########################tiff MA plot
tiff(paste(name,".tiff",sep=""),width=7,height=7,units="in",res=150)
par(mar=c(5.1, 5.1, 4.1, 2.1), bty="n")
par(xaxt="n",yaxt="n",col.axis="white")
smoothScatter(TT$AveExpr,TT$logFC,xlab="",ylab="",main="",cex.main=2,cex.lab=2,ylim=c(-10,10),nbin = 500,bandwidth=c(0.05,0.1),nrpoints = 0, colramp = colorRampPalette(c("white", "black")),xlim=c(4,18))
a<-which((TT$logFC>=diff)&(TT$adj.P.Val<=co))
b<-which((TT$logFC<=(-diff))&(TT$adj.P.Val<=co))
c<-length(unique(TT$GeneName[a]))
d<-length(unique(TT$GeneName[b]))
e<-length(unique(TT$SystematicName[a]))
f<-length(unique(TT$SystematicName[b]))
points(TT$AveExpr[a],TT$logFC[a],pch=21, bg="#d95f02",cex=0.7)
points(TT$AveExpr[b],TT$logFC[b],pch=21, bg="#7570b3",cex=0.7)
dev.off()

###########################pdf MA plot
pdf(paste(name,".pdf",sep=""),useDingbats=FALSE,width=7,height=7)
par(mar=c(5.1, 5.1, 4.1, 2.1), bty="l")
plot(1,1,type="n",col="white",xlab="Mean Expression (Log2)",ylab="Differential Expression (Log2)",main=name,cex.main=2,cex.lab=2,ylim=c(-10,10),xlim=c(4,18),axes=FALSE)
axis(1)
axis(2,at=c(-10,-5,0,5,10))
a<-which((TT$logFC>=diff)&(TT$adj.P.Val<=co))
b<-which((TT$logFC<=(-diff))&(TT$adj.P.Val<=co))
c<-length(unique(TT$GeneName[a]))
d<-length(unique(TT$GeneName[b]))
e<-length(unique(TT$SystematicName[a]))
f<-length(unique(TT$SystematicName[b]))
points(TT$AveExpr[a],TT$logFC[a],pch=21, bg="#d95f02")
points(TT$AveExpr[b],TT$logFC[b],pch=21, bg="#7570b3")
xl<-sub("-.*","",x)
yl<-sub(".*-","",x)
legend("topright",legend=c(paste("Up in ",xl," = ",length(a)," (",e,")",sep=""),paste("Up in ",yl," = ",length(b)," (",f,")",sep="")),bty="n",cex=0.7)
dev.off()

}
   


###Tidy up objects and save R image
save.image(paste(rootdir,"/R_Images/1_Summary_Array_Data.RData",sep=""))
