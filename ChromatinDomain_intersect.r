library(GenomicRanges)
library(gplots)
library(ggplot2)

setwd("~/Desktop/Project1/ATAC-seq/NS50193/Peaks/filtered_peaks/chromatin_domains/")
peaks<-lapply(dir(path="../*/",pattern="*fdr*.bed",full.names = T), read.table, header=F)
peak.names<-gsub("\\..*","", dir(path="../",pattern="fdr.*.bed"))


paths<-c("../siNT1_Dex/", "../siH1e_Dex/")
peaks<-lapply(paths, function(x) {
    z<-dir(path=x,pattern="*peaks.bed",full.names = T) 
    read.table(z, header=F) })
peak.names<-gsub("\\.bed", "", unlist(lapply(paths, function(x) dir(path=x,pattern="*peaks.bed"))))

peaks<-lapply(peaks, function(x) GRanges(x$V1, IRanges(start=x$V2, end=x$V3)))
names(peaks)<-peak.names

path<-c("~/workspace/Encode_Data/wgEncodeBroadHmm/",
        "~/workspace/Encode_Data/wgEncodeAwgSegmentation/Combined/", 
        "~/workspace/Encode_Data/wgEncodeAwgSegmentation/Segway/",
        "~/workspace/Encode_Data/wgEncodeAwgSegmentation/ChromHmm/")
for (i in seq_along(path)){
    exp<-gsub("~/.+wgEncode", "", path[i])
    exp<-gsub("/","",exp)
    BroadHmm<-lapply(dir(path=path[i],pattern=".bed",full.names = T), read.table, header=F)
    BroadHmm<-lapply(BroadHmm, function(x) GRanges(x$V1, IRanges(start=x$V2, end=x$V3), State=x$V4))
    
    overlaps<-lapply(peaks, function(z){
        lapply(BroadHmm, function(x) {
            x<-data.frame(State=mcols(x), Peaks=countOverlaps(subject=z, query=x,minoverlap = 90))
            x<-x[x$Peaks!=0,]
            x
        })
    })
    
    overlaps<-lapply(overlaps, function(y) do.call(cbind, lapply(y, function(x) tapply(x[,2],x[,1], sum))))
    
    overlaps.frac<-lapply(overlaps, 
                          function(x) data.frame(names=rownames(x), 
                                                 avg=rowMeans(t(t(x)/colSums(x))), 
                                                 stdev=apply(t(t(x)/colSums(x)), 1, sd)))
    overlaps.mean<-lapply(overlaps, 
                          function(x) data.frame(names=rownames(x), 
                                                 avg=rowMeans(x), 
                                                 stdev=apply(x, 1, sd)))
    
    for (x in seq_along(overlaps.mean))
    {
        makeBarplot(df1=overlaps.mean[[x]], exp=exp, peak=names(peaks)[x])
        makeBarplot(df1=overlaps.frac[[x]], exp=paste0(exp, "_fraction"), peak=names(peaks)[x])
    }
}

makeBarplot<-function(df1, limits, exp, peak)
{
    dodge<-position_dodge(width=0.9)
    p<-ggplot(data=df1, aes(x=names, y=avg, fill=names, ymin=avg-stdev, ymax=avg+stdev))
    p<-p+geom_bar(stat="identity", position=dodge) + 
        geom_errorbar(position=dodge, width=0.25) + 
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.x=element_blank()) +
        ggtitle(peak)
    ggsave(p, filename=paste0(exp,"_", peak, ".png"))
      
}

