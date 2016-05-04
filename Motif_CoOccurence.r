library(slam)
library(pheatmap)
library(amap)
library(RColorBrewer)

setwd("~/Desktop/A1-2_Project/Brg1_ChIP/NS50200/Peaks/Peaks_by_JH/")
brg1LD<-read.table("motifs/LabBrg1_DexPeaks_motifAnnotate.txt", header=T, 
                   quote="", allowEscapes=F, stringsAsFactors=F, sep="\t")

colnames(brg1LD)<-gsub("\\.Distance\\.From.+","",colnames(brg1LD), perl=T)
colnames(brg1LD)[1]<-"PeakID"
brg1LDmotifs<-brg1LD[,-c(1:21)]
brg1LDmotifs<-apply(brg1LDmotifs,1, gsub, 
                    pattern="^-*(\\d{1,3}).+", replacement="\\1", perl=T)
brg1LDmotifs<-t(brg1LDmotifs)
class(brg1LDmotifs)<-"numeric"
brg1LDmotifs<-ifelse(is.na(brg1LDmotifs), 0, ifelse(brg1LDmotifs>250,1,0))
brg1LDmotifs<-as.simple_triplet_matrix(brg1LDmotifs)

brg1LDcoOccur<-crossprod_simple_triplet_matrix(brg1LDmotifs)
brg1Hclust<-hcluster(brg1LDcoOccur, "euclid")

png("LabBrg1_Dexonly_Motif_coOccurence_transform.png", height=2000, width=2500, units="px")
pheatmap(sqrt(brg1LDcoOccur),color=colorRampPalette(brewer.pal(9,"Greens"))(100),
         cluster_rows = brg1Hclust, cluster_cols = brg1Hclust)
dev.off()

tree<-cutree(brg1Hclust,k = 5)
treeHeat<-tree %in% c(1,3,4)

png("LabBrg1_Dexonly_Motif_coOccurence_subset.png", height=1500, width=2500, units="px")
pheatmap(brg1LDcoOccur[treeHeat,treeHeat],border_color = NA,
         color=colorRampPalette(brewer.pal(9,"Greens"))(100))
dev.off()


brg1Peaks<-t(brg1LDmotifs)
peaksCoOccur<-crossprod_simple_triplet_matrix(brg1Peaks)
peakClust<-hcluster(peaksCoOccur, "eucli",nbproc = 16)

png("Brg1_Peak_adj_matrix.png",heigh=2500, width=3500, units="px")
pheatmap(peaksCoOccur, border_color=NA,color=colorRampPalette(brewer.pal(9,"Greens"))(100))


