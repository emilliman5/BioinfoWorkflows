library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)

bams<-dir(recursive = T, pattern=".bam$")

features<-lapply(dir("GR_ChIP/GR_peaks_Valerie/", pattern="bed$",full.names = T), 
                 function(x) read.table(x, quote=""))
names(features)<-gsub(".bed", "", dir("GR_ChIP/GR_peaks_Valerie/", pattern="bed$"))

feats<-lapply(features, function(x){
    if(x[1,dim(x)[2]] %in% c("+","-")){
        GRanges(seqnames = x$V1, ranges = IRanges(start=x$V2, end=x$V3),strand = x[,dim(x)[2]])
    } else{
        GRanges(seqnames = x$V1, ranges = IRanges(start=x$V2, end=x$V3),strand = rep("*"))
    }
})    

names(feats)<-names(features)

cov<-lapply(bams, function(x) 
    aln<-readGAlignments(x)

