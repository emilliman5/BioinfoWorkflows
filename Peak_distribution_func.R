library(reshape2)
setwd("~/Desktop/Project1/NS50005_NS50031_analysis/Analysis/")
peaks<-lapply(list.files(full.names = T, path = "~/Desktop/Project1/NS50005_NS50031_analysis/Analysis/",pattern="*closestTSS.bed"), 
              function(x) read.table(x, sep = "\t",header = F, quote=""))

names(peaks)<-gsub(pattern = paste("UL3_","NS50005_NS50031_","_closestTSS\\.bed",sep = "|"),"", list.files(pattern = "*closestTSS.bed", path="~/Desktop/Project1/NS50005_NS50031_analysis/Analysis/"))

dist<-melt(peaks,measure.vars = "V10",id.vars = c("V1","V2","V3"))

png("ATAC_peak_distance_To_TSS.png", height=800, width=1400, units="px")
boxplot(log10(abs(dist$value)+1)~as.factor(dist$L1))
dev.off()

table.summary<-sapply(peaks, function(x) summary(abs(x$V10)))


binned_dist<-sapply(peaks, function (x) cbind(sum(abs(x$V10)<=50),
                                              sum(abs(x$V10)>50 & abs(x$V10)<=2500), 
                                              sum(abs(x$V10)>2500 & abs(x$V10)<=10000),
                                              sum(abs(x$V10)>10000 & abs(x$V10)<=50000),
                                              sum(abs(x$V10)>50000)))

rownames(binned_dist)<-c("TSS","0.05-2.5 kb","2.5-10 kb","10-50 kb",">50 kb")

binned<-apply(binned_dist,2, function (x) x/sum(x))

png("ATAC_peaks_bar_chart.png",height=10,width=8, units="in",res=150)
par(mar=c(8,10,2,8), mfrow=c(2,1))
barplot(as.matrix(binned[,c(4,1,3,5,2)]), names.arg = rev(c("Shared Peaks","siNC1 Unique Peaks","siH1.4 Unique Peaks","siNC1 Peaks","siH1.4 Peaks")),
        main="ATAC peak distance to TSSs",col=c("purple","green", "blue","red","orange"), 
        legend.text=rownames(binned),width=0.14, xlim=c(0,1), xlab="Percentage of Peaks",
        horiz=T, las=2,args.legend=list(x="topright",bty="n",inset=c(-0.35,0)))
par(mgp=c(5,1,0))
barplot(as.matrix(binned_dist[,c(4,1,3,5,2)]), names.arg = rev(c("Shared Peaks","siNC1 Unique Peaks","siH1.4 Unique Peaks","siNC1 Peaks","siH1.4 Peaks")),
        main="ATAC peak distance to TSSs",col=c("purple","green", "blue","red","orange"), 
        legend.text=rownames(binned_dist),width=0.14, xlim=c(0,max(colSums(binned_dist))), xlab="Number of Peaks",
        horiz=T, las=2,args.legend=list(x="topright",bty="n",inset=c(-0.35,0)))

dev.off()

