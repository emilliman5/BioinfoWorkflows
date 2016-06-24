path<-"~/Desktop/Tech_talk/ATAC_signal_correlations//NS50005//siH1e/"

setwd(path)
autoCorr<-read.table("tagAutocorrelation.txt", sep="\t",skip=1)
plot(autoCorr$V1, autoCorr$V2, type="l", col="blue")
lines(autoCorr$V1, autoCorr$V3, type="l", col="red")

png("AutoCorrelationPlot.png", height=600, width=1000, units="px")
par(cex=2)
plot(autoCorr$V1, autoCorr$V2, type="l", col="blue", xlim=c(-500,500), xlab="Distance (bp)", ylab="Number of Tags", ylim=c(0, max(c(autoCorr$V2, autoCorr$V3))))
lines(autoCorr$V1, autoCorr$V3, type="l", col="red", xlim=c(-500,500))
legend("topright", lty=1, lwd=c(2,2), col=c("blue","red"), legend = c("Same Strand", "Opposite Strand"), bty="n")
dev.off()

GC<-read.table("genomeGCcontent.txt", skip=1)
tagGC<-read.table("tagGCcontent.txt", skip=1)

png("GC_content.png", height=600, width=1000, units="px")
    par(cex=2)
    plot(GC$V1,GC$V3, type="l", col="blue", ylab="Normalized Fraction", xlab="GC content", main="GC content")
    lines(tagGC$V1,tagGC$V3, type="l", col="red")
    legend("topright", lty=c(1,1), lwd=c(2,2), col=c("blue","red"), legend=c("hg19","tags"))
dev.off()

tagClone<-read.table("tagCountDistribution.txt", skip=1)

png("TagClonalityPlot.png", height=600, width=1000, units="px")
    par(cex=2)
    plot(tagClone$V1, tagClone$V2, type="l", lwd=2, col="purple", ylab="Fraction of Reads", xlab="Tags per postion", main="Clonal Distribution")
dev.off()

tagOffset<-read.table("tagFreqUniq.txt", header=T)
matplot(tagOffset$Offset, tagOffset[,2:5], lty=1, type="l")
legend("topright", )

matplot(tagOffset$Offset, tagOffset[,6:9], lty=1, type="l")
legend("topright", )

matplot(tagOffset$Offset, tagOffset[,10:25], lty=1, type="l")
legend("topright", )