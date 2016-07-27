files<-dir("~/R_workspace/FactorBookMotifs/",full.names = T, ".txt")
motifs<-lapply(files, read.table, skip=1)
headers<-lapply(files, read.table, nrows=1, stringsAsFactors=F)
headers<-lapply(headers, function(x) if(ncol(x)>2){
    x<-x[,c(1,3)]
    colnames(x)<-c("V1","V2")
    x
} else{x})
headers<-do.call(rbind,headers)

lapply(sample(seq_along(motifs),size = 50), function(x) {
    png(paste0("motif",x,".png"), height=300, width=600,units="px")
    seqLogo(t(motifs[[x]]),xaxis = F,yaxis = F )
    dev.off()
})

hist(unlist(lapply(motifs, nrow)), main="Motif Size")
summary(unlist(lapply(motifs, nrow)))

oddsMax<-unlist(lapply(motifs, function(x) sum(log(apply(x,1,max)/0.25))))
hist(oddsMax, main="Motif Max score", xlab="log odds")

plot(jitter(unlist(lapply(motifs, nrow))), oddsMax, xlab="Motif Size",ylab="Max Motif Score")

odds<-lapply(motifs, function(x){
        unlist(lapply(1:nrow(x), function(y){
            z<-apply(x, 1, max)
            z[y]<-min(x[y,])+0.001
            sum(log(z/0.25))
        }))
})

headers[,"V3"]<-0.5076*oddsMax+2.2353
headers$V1<-gsub(">","",paste0(headers$V1,"_",seq_along(headers$V1)))
headers$V2<-paste0(">",headers$V2)

lapply(1:nrow(headers), function(x){
    write.table(headers[x,c(2,1,3)], file="FactorBookMotifs.homer",append = T, sep="\t", quote=FALSE, col.names=F, row.names=F)
    write.table(motifs[[x]], file="FactorBookMotifs.homer",append = T, sep="\t", quote=FALSE, col.names=F, row.names=F)
})

headers[,"V3"]<-unlist(lapply(odds, median))
lapply(1:nrow(headers), function(x){
    write.table(headers[x,c(2,1,3)], file="FactorBookMotifs_medianOdds.homer",append = T, sep="\t", quote=FALSE, col.names=F, row.names=F)
    write.table(motifs[[x]], file="FactorBookMotifs_medianOdds.homer",append = T, sep="\t", quote=FALSE, col.names=F, row.names=F)
})

headers[,"V3"]<-apply(cbind(unlist(lapply(odds, median)), 0.5076*oddsMax+2.2353), 1, min)
lapply(1:nrow(headers), function(x){
    write.table(headers[x,c(2,1,3)], file="FactorBookMotifs_minOdds.homer",append = T, sep="\t", quote=FALSE, col.names=F, row.names=F)
    write.table(motifs[[x]], file="FactorBookMotifs_minOdds.homer",append = T, sep="\t", quote=FALSE, col.names=F, row.names=F)
})

boxplot(odds[1:25], xlab="Motif",ylab="log odds score with one mismatch", ylim=c(0,max(oddsMax[1:25])))
points(oddsMax, col="blue",pch=19)
points(as.numeric(headers[1:25,3]), col="red",pch=20)
points(0.5076*oddsMax[1:25]+2.2353, col="green",pch=20)
legend(12,19,pch=c(19, 20, 20),col=c("blue","red", "green"), legend=c("Max log odds score", "Minimum Odds Threshold", "Homer threshold"), bty="n")
