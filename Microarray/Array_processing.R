library(limma)
library(lattice)
library(gplots)
library(RColorBrewer)
##functions

logratio2foldchange <- function (logratio, base=2)
{
    retval <- base^(logratio);
    retval <- ifelse(retval < 1, -1/retval, retval);
    retval;
}

## Read in Data and munge Array key into workable form.
targets<-readTargets("941_Archer/941_Archer_Array key.csv",sep = ",")
targets<-targets[1:15,1:2]
t(as.data.frame(strsplit(targets$observation,split = "_")))[,1]
targets$FileName<-paste("US22502532",t(as.data.frame(strsplit(targets$observation,split = "_")))[,1],"S01","GE1","1200","Jun14",
                        t(as.data.frame(strsplit(targets$observation,split = "_")))[,2],t(as.data.frame(strsplit(targets$observation,split = "_")))[,3], sep = "_")
targets$FileName<-paste(targets$FileName,"txt",sep = ".")
targets
x<-read.maimages(targets, source="agilent",path="941_Archer/", green.only = T)

##Diagnotics

png("Milliman_941_MAs.png", height=1500, width=900)
par(mfrow=c(5,3))
for(i in 1:dim(targets)[1]){
    smoothScatter(log2(rowMeans(x$E)), log2(x$E[,i])-log2(rowMeans(x$E[,-i])), main=targets[i,2], ylab="Expression log ratio (this sample vs. others)",
                  xlab="Average Expression", ylim=c(-6,6))
}
dev.off()

##Normalization
y<-backgroundCorrect(x,method="normexp")
y<-normalizeBetweenArrays(y,method="quantile")

##Post Normalization Diagnotics

png("Milliman_941_boxplots.png", height=1000, width=800)
par(mfrow=c(2,1))
boxplot(log2(x$E), names = targets[,2], ylab="log2 Intensity", main="Array Intensities PreNormalization")
boxplot(y$E, names = targets[,2], ylab="Quantile Normalized Intensity", main="Array Intensities Post Normalization")
dev.off()

png("Milliman_941_MAs_post_normal.png", height=1500, width=900)
par(mfrow=c(5,3))
for(i in 1:dim(targets)[1]){
    smoothScatter(log2(rowMeans(y$E)), log2(y$E[,i])-log2(rowMeans(y$E[,-i])), main=targets[i,2], ylab="Expression log ratio (this sample vs. others)",
                  xlab="Average Expression", ylim=c(-4,4))
}
dev.off()

##Setup Design and Contrasts Matrices for differential expression analysis
design<-data.frame("Filename"=targets$FileName, "Treatment"=c(0,1,1,1,1,1,1,0,0,0,1,1,1,0,0), "Condition"=c(1,1,1,1,1,1,1,0,0,0,0,0,0,1,1),
                   "Time"=c(4,4,4,4,24,24,24,4,4,4,4,4,4,4,4))
design[design$Treatment==0,2]<-"Veh"
design[design$Treatment==1,2]<-"Dex"
design[design$Condition==0,3]<-"siNT1"
design[design$Condition==1,3]<-"siH1e2"
design
TS<-paste(design$Condition,design$Treatment,design$Time, sep=".")
TS<-factor(TS)
TS

##My way for a model matrix
des<-model.matrix(~0+TS)
rownames(des)<-design$Filename

##James' way to create the model matrix
des2<-do.call(cbind,lapply(levels(TS), function(i){ (TS %in% i)*1}))
colnames(des2)<-levels(TS)
rownames(des2)<-design$Filename

des2

contrast.matrix<-makeContrasts(
    siNT1DexvsVeh=siNT1.Dex.4-siNT1.Veh.4,
    siH1eDexvsVeh=siH1e2.Dex.4-siH1e2.Veh.4,
    siH1e24hrDexvsVeh=siH1e2.Dex.24-siH1e2.Veh.4,
    siH1evssiNT1=siH1e2.Veh.4-siNT1.Veh.4,
    siH1e24vs4=siH1e2.Dex.24-siH1e2.Dex.4,
    siH1eDexvssiNT1Dex=(siH1e2.Dex.4-siH1e2.Veh.4)-(siNT1.Dex.4-siNT1.Veh.4),
    levels=des2
)
contrast.matrix

## Remove lowly expressed and control probes

# neg95<-apply(y$E[y$genes$ControlType==-1,],2, function(x) quantile(x,p=0.95))
# cutoff<-matrix(1.1*neg95,nrow(y), ncol(y), byrow=T)
# isexpr<-rowSums(y$E >cutoff) >=3
# table(isexpr)

y0<-y[y$genes$ControlType==0,]

##Fit model to expression data
fit.p<-lmFit(y0, des2)
fit.p2<-eBayes(fit.p,trend=T)

##Its a plot with points and lines... I think it is telling me that variation is higher at lower expression levels
png("SigmaVsAvgExp.png", height=500, width=1000, units="px")
plotSA(fit.p2,main="Probe-level S-A Plot")
dev.off()

summary(decideTests(fit.p2[,-1]))

##Gene-level averaging
# yAve<-avereps(y0,ID=y0$genes[,"SystematicName"])
# fit.g<-lmFit(yAve, des2)
# fit.g2<-eBayes(fit.g, trend=T)
# plotSA(fit.g2,main="Gene-level")

#probe-level contrasts
fit.p.c<-contrasts.fit(fit.p, contrast.matrix)
fit.p.c2<-eBayes(fit.p.c, trend=T)

results<-decideTests(fit.p.c2, lfc=.75, p.value = 0.01)
summary(results)
table(results[,1], results[,2])

png("Probe-level_vennDiagrams.png", height=1000, width=1000, units="px")
par(mfrow=c(2,2), oma=c(0,0,2,0))
vennDiagram(results[,c(1,2,6)])
vennDiagram(results[,c(2,3,5)])
vennDiagram(results[,c(1,2,4)])
mtext("adjPvalue ≤ 0.01 & logFC ≥0.75", outer=T, cex=2.5)
dev.off()

coefNames<-colnames(coef(fit.p.c2))
names(coefNames)<-coefNames

cutoffAveExpr <- 0;
cutoffAdjPVal <- 0.01;
cutoffLogFC <- log2(2);
cutoffPVal <- 0.01;

probe.hits<-lapply(coefNames, function(i){
    iTopTable <- topTable(fit.p.c2, coef=i, sort.by="none", number=nrow(results));
    iTopTable[,"hit"] <- (iTopTable[,"AveExpr"] >= cutoffAveExpr & 
                              iTopTable[,"adj.P.Val"] <= cutoffAdjPVal & 
                              abs(iTopTable[,"logFC"]) >= cutoffLogFC)*
        sign(iTopTable[,"logFC"]);
    iTopTable[,"hitP"] <- (iTopTable[,"AveExpr"] >= cutoffAveExpr & 
                               iTopTable[,"P.Value"] <= cutoffPVal & 
                               abs(iTopTable[,"logFC"]) >= cutoffLogFC)*
        sign(iTopTable[,"logFC"]);
    iTopTable[,"fold"] <- logratio2foldchange(iTopTable[,"logFC"]);
    colnames(iTopTable) <- paste(colnames(iTopTable), i);
    iTopTableDF <- data.frame(check.names=FALSE, "Region"=topTable(fit.p.c2, coef=i, sort.by="none", number=nrow(results))$SystematicName,
                              iTopTable[,grep(ignore.case = T, paste(c("hit", "p.val", "logFC", "fold"),collapse="|"), colnames(iTopTable)),drop=FALSE]);
});

for (i in 1:length(probe.hits)){
    png(paste0(names(probe.hits)[i],"volcanoPlot.png"),height=1250, width=3000, units="px")
    par(mfrow=c(1,2))
    smoothScatter(probe.hits[[i]][,2], -log10(probe.hits[[i]][,3]), main=names(probe.hits)[i],cex.main=2, cex.lab=1.5,cex.axis=1.5, xlab="Log2 Fold Change", ylab="-log10 pValue")
    smoothScatter(probe.hits[[i]][,2], -log10(probe.hits[[i]][,4]), main=names(probe.hits)[i], cex.main=2,cex.lab=1.5, cex.axis=1.5, xlab="Log2 Fold Change", ylab="-log10 Adj. pValue")
    dev.off()
}

png("LFC_scatterPlots.png", height=1500, width=1500, units="px")
splom(do.call(cbind, lapply(probe.hits, function(x) { x[,2]})), panel=panel.smoothScatter)
dev.off()

probeHits.df<-do.call(cbind, lapply(probe.hits, function(x) { x[,2]}))

isDiff<-rowSums(abs(probeHits.df)>=matrix(1, nrow(probeHits.df), ncol(probeHits.df)))>=1
n<-100
my_palette<-colorRampPalette(c("red","white", "blue"))(n=n-1)
col_breaks<-c(seq(-5,-1,length=n/3),seq(-1,1, length=n/3),seq(1,5,length=n/3))

png("siH1e_Dex_probeLevel_heatmap.png", height=1200, width=800, units="px")
heatmap.2(probeHits.df[isDiff,], trace="none", col=my_palette, margins =c(18,9),labRow = "", breaks=seq(-4,4, length.out = 100), lhei = c(.75,4))
dev.off()

png("siH1e_Dex_heatmap_4Poster.png", height=1200, width=800, units="px")
heatmap.2(probeHits.df[isDiff,c(1,2,4)], trace="none", 
          labCol = c("siNT1 Dex vs. Veh","siH1.4 Dex vs. Veh", "siH1.4 Veh vs.\n siNT1 Veh"), 
          col=my_palette, margins =c(22,9),labRow = "", breaks=seq(-4,4, length.out = 100), lhei = c(.75,4))
dev.off()

for (i in 1:length(probe.hits)){
    write.table(probe.hits[[i]][,c(1,2,3,4,7)] ,file=paste0(names(probe.hits)[i],"_expression.txt"),
              quote=F, sep="\t", row.names=F, col.names=T)
    
    write.table(probe.hits[[i]][probe.hits[[i]][,5]==1,c(1,2,3)],
                                file=paste0(names(probe.hits)[i],"_DEG_expression.txt"), 
                                quote=F,sep="\t", row.names=F, col.names=T)
}

##Gene-level contrasts
# fit.g.c<-contrasts.fit(fit.g, contrast.matrix)
# fit.g.c2<-eBayes(fit.g.c, trend=T)
# 
# summary(decideTests(fit.g.c2, lfc = 1,adjust.method = "none", p.value = 0.01))
# 
# results<-decideTests(fit.g.c2, lfc=1, adjust.method = "none", p.value = 0.01)
# vennDiagram(results[,c(1,2,6)])
# vennDiagram(results[,c(1,2)])

