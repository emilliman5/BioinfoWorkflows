#!/usr/bin/Rscript
suppressMessages(library(docopt))
suppressMessages(library(GenomicRanges))
suppressMessages(library(VennDiagram))
source("~/Documents/methods.R")

doc<-" This script will generate a Venn diagram for two sets of peaks in bed format as well as output a new bed file with the shared and unique peaks separated so that they can be fed into deeptools for making heatmaps.

Usage:  makePeakVenn.R --p1=<peak1> --p2=<peak2> [--n1=<name1>] [--n2=<name2>]

Options:
    --p1=<peak1>           Peak file 1
    --p2=<peak2>           Peak file 2
    --n1=<name1>           Name to use for Peak file 1 (Defaults to file name) 
    --n2=<name2>           Name to use for Peak file 2 (Defaults to file name)
    -h --help             This helpful message"

my_opts<-docopt(doc)
print(my_opts)
 
peaks<-lapply(c(my_opts$p1, my_opts$p2), read.table)
peaks<-lapply(peaks, bed2Grange)

if(is.null(my_opts$n1)){
    my_opts$n1<-gsub("*/","", my_opts$p1)
    my_opts$n1<-gsub("\\..*$","", my_opts$n1)
}

if(is.null(my_opts$n2)){
    my_opts$n2<-gsub("*/","", my_opts$p2)
    my_opts$n2<-gsub("\\..*$","", my_opts$n2)
}

print(my_opts)

names(peaks)<-c(my_opts$n1, my_opts$n2)

x<-peakOverlap(peaks[[1]],peaks[[2]])

png(paste0(my_opts$n1, "_",my_opts$n2, "_Venn.png"),height=900, width=900, units="px")
draw.pairwise.venn(x[1],x[2],x[3], fill = c("red4", "blue4"),cat.cex = c(2,2),cex = c(2,2,2),
                   cat.fontfamily = "sans serif",
                   category = c(my_opts$n1, my_opts$n2),
                   col = c("red","blue"),cat.pos=c(220,45),
                   print.mode = c("raw","percent"))
dev.off()
 
sampleNames<-c("Shared",paste0(my_opts$n1, " Unique"),paste0(my_opts$n2, " Unique"))
peakSets<-peakLists(peaks[[1]],peaks[[2]])
peakSets<-lapply(peakSets, function(x) data.frame(seqnames=seqnames(x),
                                            start=start(x),
                                            end=end(x),
                                            names=names(x),
                                            score=0,
                                            strand=gsub("\\*",".",strand(x)))
)

lapply(seq_along(peakSets), function(x) {
    write.table(peakSets[[x]], row.names = F,col.names = F,sep="\t", file=paste0(my_opts$n1, "_",my_opts$n2, "_venn.bed"),quote = F, append=T)
    write(paste0("#",sampleNames[x]), file=paste0(my_opts$n1, "_",my_opts$n2, "_venn.bed"), append = T)
})

# 
# grPeaks<-lapply(dir("../../../GR_ChIP_seq/Peaks_merged_runs/", 
#                     recursive = T,full.names = T,pattern="narrowPeak"),
#                 read.table)
# grPeaks2<-lapply(grPeaks, bed2Grange)
# names(grPeaks2)<-gsub(".narrowPeak","", dir("../../../GR_ChIP_seq/Peaks_merged_runs/", 
#                                             recursive = T,pattern="narrowPeak"))
# names(grPeaks2)<-gsub(".*/","",names(grPeaks2))
# 
# grvenn<-peakOverlap(grPeaks2[[3]],grPeaks2[[1]])
# 
# png("GR_mergerd_peaks_venn.png", height=800, width=800, units="px")
# draw.pairwise.venn(grvenn[1],grvenn[2],grvenn[3], 
#                    category = names(grPeaks2)[c(3,1)],
#                    col = c("red","blue"),
#                    print.mode = c("raw","percent"),
#                    fill = c("red2","blue2"))
# dev.off()
# 
# grvenn10<-peakOverlap(grPeaks2[[3]][grPeaks2[[3]]@elementMetadata$score>=quantile(grPeaks2[[3]]@elementMetadata$score,0.9)],
#                       grPeaks2[[1]][grPeaks2[[1]]@elementMetadata$score>=quantile(grPeaks2[[1]]@elementMetadata$score,0.9)])
# 
# png("GR_mergerd_Top10_peaks_venn.png", height=800, width=800, units="px")
# draw.pairwise.venn(grvenn10[1],grvenn10[2],grvenn10[3], 
#                    category = names(grPeaks2)[c(3,1)],
#                    col = c("red","blue"),
#                    print.mode = c("raw","percent"),
#                    fill = c("red2","blue2"))
# dev.off()
# 
# boxplot(peaks2[[1]]@elementMetadata$score,peaks2[[2]]@elementMetadata$score,peaks2[[3]]@elementMetadata$score, range = 0,names = names(peaks2))
# 
# grPeaks3<-peakLists(grPeaks2[[3]],grPeaks2[[1]])
# names(grPeaks3)<-c("Shared_GR_peaks","siCtrl_GR_unique_peaks", "siH1.4_GR_unique_peaks")
# 
# lapply(seq_along(grPeaks3), function(x) {
#     write.table(x=data.frame(seqnames=seqnames(grPeaks3[[x]]),
#                              start=start(grPeaks3[[x]]),
#                              end=end(grPeaks3[[x]]),
#                              names=names(grPeaks3[[x]]),
#                              strand=gsub("\\*",".",strand(grPeaks3[[x]]))), 
#                 file=paste0(names(grPeaks3)[x],".bed"), row.names=F, col.names=F)
#     
# })
# 
