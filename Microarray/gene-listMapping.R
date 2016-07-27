##Script to create Bed files of gene names, corrdinates and expression
##This is a supplement script to the analysis of Milliman_941 microarray experiment
##Gene - TSS mappings were from Kim's START-seq experiment and produced by the Integrative Bioinformatics Group (Andy Lavender)

library(reshape2)

RefSeq<-read.table("/Volumes/millimanej-1/My Documents/Project 1/Features/RefSeq_annotated_oTSS_Veh.bed", header = F,sep = "\t",quote = "")

geneList<-read.table("2015may06/siH1e2_24hrDEX_totalExpression.txt",header = F,quote="")

geneList<-merge(RefSeq,geneList, by.x = RefSeq$V4, by.y=geneList$V1)
