library(ChIPQC)

setwd("~/Desktop/Project1/GR_ChIP_seq/")
files<-c(dir(path = "NS50038/Bams_merged/",pattern = ".bam$",full.names = T),
         dir(path = "NS50180/Bams/",pattern = ".bam$", full.names = T))

ids<-files
ids<-gsub("MES661_141114_M00421_ABFD6.", "",ids)
ids<-unlist(lapply(strsplit(ids, split = "\\."), function(x) x[1]))
samples<-data.frame(SampleID=ids ,rep=c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2), CellType=rep("UL3"), bamReads=files, 
                    Peaks="~/Desktop/Project1/GR_ChIP_seq/NS50038/Macs_output/UL3_GR_peaks_all.bed")

samples$SampleID<-makeUnique(samples$SampleID)
samples$Treatment<-rep("Veh")
samples$Treatment[grep("Dex", samples$SampleID)]<-"Dex"
samples$Condition<-rep("siNC1")
samples$Condition[grep("siH1e", samples$SampleID)]<-"siH1e"
samples$Condition[grep("siH1-4", samples$SampleID)]<-"siH1e"
samples$Factor<-rep("Input")
samples$Factor[grep("IP", samples$SampleID)]<-"GR"

qc<-ChIPQC(samples, annotation="hg19",chromosomes = c("chr1", "chr3","chr5","chr7","chr9","chr11","chr13"))
QCmetrics(qc)

ChIPQCreport(qc,facetBy = c("Treatment","Factor"))
write.table(QCmetrics(qc), "ChIP_qc.txt",sep="\t",F)

plotCoverageHist(qc, facetBy=c("Treatment","Factor"))
plotCC(qc, facetBy=c("Treatment","Factor"))
plotPeakProfile(qc, facetBy=c("Treatment","Factor"))
plotRegi(qc, facetBy=c("Treatment","Factor"))
plotCorHeatmap(qc, attributes=c("Treatment","Factor"))
plotPrincomp(qc, attributes=c("Treatment","Factor"))
