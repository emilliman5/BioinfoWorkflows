library(reshape2)
library(limma)
library(lattice)
library(gplots)
library(RColorBrewer)
library(VennDiagram)
library(gtools)
library(fields)
library(made4)

##functions

extraFunFile <- "Array_processing_functions.R";
if (file.exists(extraFunFile)) {
   source(extraFunFile, keep.source=TRUE);
}

dir.create(paste0("../results/",getDate()))

## Read in Data and munge Array key into workable form.
targets<-readTargets("../data/972_Archer/972_Archer_Array key.csv",sep = ",")
targets<-targets[1:12,1:2]
t(as.data.frame(strsplit(targets$observation,split = "_")))[,1]

targets$FileName<-paste("US22502532",t(as.data.frame(strsplit(targets$observation,split = "_")))[,1],"S01","GE1","1200","Jun14",
                        t(as.data.frame(strsplit(targets$observation,split = "_")))[,2],t(as.data.frame(strsplit(targets$observation,split = "_")))[,3], sep = "_")
targets$FileName<-paste(targets$FileName,"txt",sep = ".")
targets
x<-read.maimages(targets, source="agilent",path="../data/972_Archer/", green.only = T)

## Add groupName and color by sample group
targets$groupName <- gsub(" Replicate.+", "", targets$BioRep);
targets$color <- brewer.pal(9, "Set1")[as.numeric(as.factor(targets$groupName))];

## I wrappered up the code above into a function, so you can re-run it after normalization
makeMAplots <- function
(x, targets, xlab="Average Expression", ylab="Expression log ratio (this sample vs. others)",
 ylim=c(-3,3), ablineH=c(-1,0,1), ncolX=NULL, nrowY=NULL,
 nbin=128, bandwidthN=300, useRaster=TRUE,
 ...)
{
   ## Purpose is to re-use a simple MA-plot wrapper
   if (max(x$E, na.rm=TRUE) > 30) {
      xE <- log2(x$E);
   } else {
      xE <- x$E;
   }
   if (is.null(ncolX)) {
      ncolX <- ceiling(sqrt(ncol(xE)));
   }
   if (is.null(nrowY)) {
      nrowY <- ceiling(ncol(xE) / ncolX);
   }
   parMfrow <- par("mfrow");
   par("mfrow"=c(nrowY, ncolX));
   for(i in 1:dim(targets)[1]){
       plotSmoothScatter(rowMeans(xE), xE[,i]-rowMeans(xE[,-i]), main=targets[i,2], ylab=ylab,
                     xlab=xlab, ylim=ylim, nbin=nbin, bandwidthN=bandwidthN, useRaster=useRaster, ...);
       abline(h=ablineH, lty="dashed", col="grey40");
   }
   par("mfrow"=parMfrow);
}

png(paste0("../results/",getDate(),"/Milliman_smoothMA_preNormal.png"), height=1500, width=900)
makeMAplots(x, targets, xlab="Avg Expression (raw data)");
dev.off();

## Alternate normalization without backgroundCorrect()
y<-normalizeBetweenArrays(x,method="quantile")
png(paste0("../results/",getDate(),"/Milliman_972_smoothMA_quantNormal.png"), height=1500, width=900)
makeMAplots(y, targets, xlab="Avg Expression (quantile-normalized)");
dev.off();

png(paste0("../results/",getDate(),"/Milliman_972_boxplot_quantNormal.png"), height=1500, width=900)
par(mfrow=c(2,1), oma=c(5,2,2,2))
boxplot(log2(x$E), names=gsub("Replicate", "Rep",targets$BioRep), las=2, ylab="Expression (log2)")
boxplot(y$E,names = gsub("Replicate", "Rep",targets$BioRep), las=2 ,ylab="Expression (quantile-normalized)");
dev.off();

## Sample correlation analysis
yCtr <- centerGeneData(y$E, returnGroupedValues=FALSE, mean=TRUE);
yCtrCor <- cor(yCtr);

png(paste0("../results/",getDate(),"/Milliman_972_sampleCorrelations_quantNormal.png"), height=800, width=800, units="px")
heatmap.2(yCtrCor, ColSideColors=targets$color, RowSideColors=targets$color, trace="none", tracecol="orange",
   main="Sample correlation\n(quantile-normalized)", #cellnote=round(yCtrCor, digits=2),
   labRow=targets$BioRep, labCol=targets$BioRep, margins=c(15,15), col=colorRampPalette(rev(brewer.pal(11, "RdBu"))));
dev.off()

png(paste0("../results/",getDate(),"/Milliman_972_sampleCorrelations_raw.png"), height=800, width=800, units="px")
xCtr <- centerGeneData(log2(x$E), returnGroupedValues=FALSE, mean=TRUE);
xCtrCor <- cor(xCtr);
heatmap.2(xCtrCor, ColSideColors=targets$color, RowSideColors=targets$color, trace="none", tracecol="orange",
   main="Sample correlation\n(raw data)", #cellnote=round(xCtrCor, digits=2),
   labRow=targets$BioRep, labCol=targets$BioRep, margins=c(15,15), col=colorRampPalette(rev(brewer.pal(11, "RdBu"))));
dev.off()

##Setup Design and Contrasts Matrices for differential expression analysis
design<-data.frame("Filename"=targets$FileName, "Treatment"=c(1,1,0,0,0,0,0,1,0,1,1,1), "Condition"=c(0,0,1,1,0,0,0,0,1,1,1,1),
                   "Time"=c(4,4,4,4,4,4,4,4,4,4,4,4))
design[design$Treatment==0,2]<-"Veh"
design[design$Treatment==1,2]<-"Dex"
design[design$Condition==0,3]<-"siNT1"
design[design$Condition==1,3]<-"siH1e2"
design$Treatment <- factor(design$Treatment, levels=c("Veh", "Dex"));
design
design$Time <- factor(design$Time);
TS<-paste(design$Condition,design$Treatment,design$Time, sep=".");
TS<-factor(TS);
TS;

##James' way to create the model matrix
des2<-do.call(cbind,lapply(levels(TS), function(i){ (TS %in% i)*1}));
colnames(des2)<-levels(TS);
rownames(des2)<-design$Filename;
des2

contrast.matrix<-makeContrasts(
    siNT1DexvsVeh=siNT1.Dex.4-siNT1.Veh.4,
    siH1eDexvsVeh=siH1e2.Dex.4-siH1e2.Veh.4,
    siH1evssiNT1=siH1e2.Veh.4-siNT1.Veh.4,
    siH1eDexvssiNT1Dex=(siH1e2.Dex.4-siH1e2.Veh.4)-(siNT1.Dex.4-siNT1.Veh.4),
    levels=des2
)
contrast.matrix

y1<-y[y$genes$ControlType==0,];
rownames(y1$genes) <- make.unique(y1$genes$ProbeName);

##Average Duplicate identical probes
dup.probes<-duplicated(y1$genes$ProbeName) | duplicated(y1$genes$ProbeName, fromLast=T)
tmp<-y1[dup.probes,]
tmp$E<-aggregate(tmp$E,list(as.factor(tmp$genes$ProbeName)),median)
tmp$Eb<-aggregate(tmp$Eb,list(as.factor(tmp$genes$ProbeName)),median)

tmp$genes<-tmp$genes[!duplicated(tmp$genes$ProbeName),]
tmp$genes<-tmp$genes[order(tmp$genes$ProbeName),]
# if(tmp$genes$ProbeName==as.character(tmp$E$Group.1)){
#     break;
# }
tmp$E<-tmp$E[,-1]
tmp$Eb<-tmp$Eb[,-1]
y0<-rbind(y1[!dup.probes,],tmp)
rm(tmp)
####

##Fit model to expression data
fit.p<-lmFit(y0, des2);
fit.p2<-eBayes(fit.p,trend=TRUE);

##Its a plot with points and lines... I think it is telling me that variation is higher at lower expression levels
png(paste0("../results/",getDate(),"/Milliman_972_SigmaVsAvgExp.png"), height=600, width=800, units="px")
plotSAsmooth(fit.p2,main="Probe-level S-A Plot")
dev.off()

#probe-level contrasts
fit.p.c<-contrasts.fit(fit.p, contrast.matrix);
fit.p.c2<-eBayes(fit.p.c, trend=TRUE);

cutoffAveExpr <- 0;
cutoffAdjPVal <- 0.05;
cutoffLogFC <- log2(2);
cutoffPVal <- 0.01;

results <- decideTests(fit.p.c2, lfc=cutoffLogFC, p.value=cutoffAdjPVal);
summary(results);

# tt1 <- eBayes2TopTables(fit.p.c2, cutoffPVal=0.01, cutoffFold=2^0.75, mergeDF=FALSE);
# tt <- eBayes2TopTables(fit.p.c2, cutoffPVal=0.01, cutoffFold=2^0.75, mergeDF=TRUE);

venn<-vennCounts(results[,c(1,2,4)])
venn.diagram(imagetype="png",cat.dist=c(0.1,0.1, 0.05) ,x = list(rownames(results[results[,1]!=0,]), rownames(results[results[,2]!=0,]),rownames(results[results[,4]!=0,])),
             category.names = c("siNT1 Dex vs Veh","siH1.4 Dex vs Veh", "(siNT1 Dex vs Veh) vs (siH1.4 Dex vs Veh)"),margin=c(0.2,0.2),
             filename=paste0("../results/",getDate(),"/Dex_probe_venn.png"), fill=c("firebrick4","darkblue", "darkorange"), 
             alpha=c(0.5,0.5, 0.5))

venn.diagram(imagetype="png",cat.dist=c(0.1,0.1, 0.05) ,x = list(rownames(results[results[,1]!=0,]), rownames(results[results[,2]!=0,]),rownames(results[results[,3]!=0,])),
             category.names = c("siNT1 Dex vs Veh","siH1.4 Dex vs Veh", "siNT1 vs siH1.4") ,margin=c(0.2,0.2),
             filename=paste0("../results/",getDate(),"/Dex_siH1_probe.png"), fill=c("firebrick4","darkblue", "darkorange"), 
             alpha=c(0.5,0.5, 0.5))

coefNames<-colnames(coef(fit.p.c2))
names(coefNames)<-coefNames

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
    png(paste0("../results/",getDate(),"/",names(probe.hits)[i],"_volcanoPlot.png"),height=1250, width=3000, units="px");
    par(mfrow=c(1,2));par("mar"=c(5,5,4,2));
    plotSmoothScatter(probe.hits[[i]][,2], -log10(probe.hits[[i]][,3]), main=names(probe.hits)[i],cex.main=2, cex.lab=1.5,
       cex.axis=1.5, xlab="Log2 Fold Change", ylab="-log10 pValue", ylim=c(0,6), xlim=c(-3,3), useRaster=TRUE)
    plotSmoothScatter(probe.hits[[i]][,2], -log10(probe.hits[[i]][,4]), main=names(probe.hits)[i], cex.main=2,cex.lab=1.5,
       cex.axis=1.5, xlab="Log2 Fold Change", ylab="-log10 Adj. pValue", ylim=c(0,6), xlim=c(-3,3), useRaster=TRUE)
    dev.off()
}

png(paste0("../results/",getDate(),"/Milliman_972_LFC_scatterPlots.png"), height=1500, width=1500, units="px")
pairs(do.call(cbind, lapply(probe.hits, function(x) { x[,2]})), diag.panel=panel.density, panel=function(...) plotSmoothScatter(...,add=T))
dev.off()

pids<-rownames(probe.hits[[3]][probe.hits[[3]]$`hitP siH1evssiNT1`!=0,])
prIntensity<-y0$E[y0$genes$ProbeName %in% pids, grep("Veh", y0$targets$BioRep)]
colnames(prIntensity)<-gsub("H1e1","H1.4", grep("Veh", y0$targets$BioRep, value=T))
colnames(prIntensity)<-gsub("_", " ", colnames(prIntensity))
colnames(prIntensity)<-gsub("Veh Replicate", "Rep", colnames(prIntensity))
prFC<-probe.hits[[3]][probe.hits[[3]]$`hitP siH1evssiNT1`!=0,2]
rowSideColors<-colorRamps::matlab.like2(99)[cut(prFC,100)]

png(paste0("../results/",getDate(),"/sih1e_vssiNT1_rep_heatmap2.png"), height=1500, width=800, units="px")
heatmap.2(prIntensity, trace="none", labRow="", margins =c(24,9), 
          ColSideColors = targets$color[grep("Veh", y0$targets$BioRep)],
          col=grDevices::heat.colors(32),density.info = "none",keysize = 1,
          RowSideColors = rowSideColors,
          key.xlab = "Probe Intensity",cexCol = 3, key.par=list(mar=c(5,0.25,6,0.25), cex=1.5))
image.plot(as.matrix(prFC),horizontal=T,col=colorRamps::matlab.like2(99),
           legend.only = T,legend.lab = "log2 Fold Change", legend.shrink=0.5)
dev.off()

probeHits.df<-do.call(cbind, lapply(probe.hits, function(x) { x[,2]}))
probeHits.df<-cbind(probe.hits[[2]][,"Region"], probeHits.df)

isDiff<-rowSums(abs(probeHits.df)>=matrix(1, nrow(probeHits.df), ncol(probeHits.df)))>1
n<-100
my_palette<-colorRampPalette(rev(brewer.pal(12, "RdBu")))(n=n-1)
col_breaks<-c(seq(-5,-1,length=n/3),seq(-1,1, length=n/3),seq(1,5,length=n/3))

png(paste0("../results/",getDate(),"/siH1e_AllDex_probeLevel_heatmap.png"), height=1200, width=800, units="px")
heatmap.2(probeHits.df[isDiff,-1], trace="none", col=my_palette, margins =c(18,9),labRow = "", breaks=seq(-4,4, length.out = 100), lhei = c(.75,4))
dev.off()

png(paste0("../results/",getDate(),"/siH1e_Dex_probeLevel_heatmap.png"), height=1200, width=800, units="px")
heatmap.2(probeHits.df[isDiff,c(-1,-5)], 
          trace="none", col=my_palette, 
          margins =c(18,9),labRow = "", 
          breaks=seq(-4,4, length.out = n), 
          lhei = c(.75,4),
          labCol=c("siNT1 Dex vs. Veh","siH1.4 Dex vs. Veh","siH1.4 vs. siNT1"))
dev.off()

################################################
##Duplicate gene probe QC
################################################
# rownames(probe.hits)<-probe.hits$Region

dup.genes<-lapply(probe.hits, function(x) duplicated(x[,1]) | duplicated(x[,1], fromLast = T))
gene.dups<-mapply(x=probe.hits, y=dup.genes, function(x,y) x[y,], SIMPLIFY = FALSE)
gene.dups.median<-lapply(gene.dups, function(x) 
    tapply(x[,2], droplevels(x[,"Region"]), median))

gene.dups.median<-lapply(gene.dups.median, function(x) data.frame(Region=names(x), Median=x))
gene.dups<-mapply(x=gene.dups,y=gene.dups.median, function(x,y) 
    merge(x,y,by="Region", all.x=T), SIMPLIFY =F)

hist(tapply(droplevels(gene.dups[[1]][,1]),droplevels(gene.dups[[1]][,1]), length))
summary(tapply(droplevels(gene.dups[[1]][,1]),droplevels(gene.dups[[1]][,1]), length))

png(paste0("../results/",getDate(),"/Replicate_Gene_Probe_vs_Median_Gene_scatter.png"), height=1200, width=2100, units="px")
par(mfrow=c(2,3))
lapply(gene.dups, function(x) plotSmoothScatter(cex=1.5, x[,2],x[,"Median"], main=colnames(x)[2], ylab="Individual Probe log2", xlab="Median log2 Fold Change per Gene"))
dev.off()

rm(tmp, gene.dups, dup.genes, dup.probes)

##############################################################################################
##Gene Venn Diagrams
##############################################################################################

gene.results<-as.data.frame(cbind(probe.hits[[1]][,"Region"],do.call(cbind, lapply(probe.hits, function(x) as.numeric(x[,5])))))
colnames(gene.results)<-c("Region", names(probe.hits))

head(gene.results)

threeWayVenn(gene.results[,c(2,3,4)], filename = "siH1_Dex_gene_venn.png", categoryNames = c("siNT1 Dex vs. Veh","siH1.4 Dex vs. Veh", "siH1.4 Veh vs. siNT1 Veh"))
threeWayVenn(gene.results[,c(2,3,5)], filename = "Dex_dd_gene_venn.png", categoryNames = c("siNT1 Dex vs Veh","siH1.4 Dex vs Veh", "(siNT1 Dex vs Veh) vs (siH1.4 Dex vs Veh)"))

png(paste0("../results/",getDate(),"/Pairwise_venn.png"), height=800, width=800, units="px")
draw.pairwise.venn(xx[2,1]+xx[2,2], xx[2,2]+xx[1,2],xx[2,2],
                   fill=c("firebrick4","darkblue"), 
                   category = c("siNT1 Dex vs Veh", "siH1.4 Dex vs. Veh"), 
                   alpha=c(0.5,0.5), cex = 2, cat.cex=2, cat.dist = 0.1, margin=c(0.2,0.2))
dev.off()

###########################################################
## Gene info Retreival
###########################################################

library(org.Hs.eg.db)
library(biomaRt)
##Create a list of all unique Probes in the experiment for unified look-up
##We will iterativly search the look up tables for the two sets of Ids in the df probes.
##1. Agilent IDs from biomaRt
##2. Refseq IDs
##3. Ensembl IDs
##4. Genbank Accessions

probes<-unique(do.call(rbind, lapply(probe.hits, function(x) data.frame(probeID=rownames(x),Region=as.character(x[,"Region"])))))
rownames(probes)<-probes$probeID

##Use biomaRt to get Agilent probe ID to Entrez Gene ID mappings
##Other mappings can be added here to query biomaRt, however the 
##values parameter inthe getBM will need to be updated to include the necessary id data
ensembl<-useMart("ensembl", dataset="hsapiens_gene_ensembl")
attrs<-c("efg_agilent_wholegenome_4x44k_v1","entrezgene")
patterns<-list(c("A_\\d\\d_","efg_agilent_wholegenome_4x44k_v1"))
agilent<-do.call("rbind",lapply(patterns, function(x) getBM(attributes = attrs, filters=x[2],
                                                               values=rownames(probe.hits[[1]]), 
                                                               mart=ensembl)))

genbank<-as.data.frame(org.Hs.egACCNUM2EG)
ensembl<-as.data.frame(org.Hs.egENSEMBLTRANS2EG)
refseq<-as.data.frame(org.Hs.egREFSEQ2EG)

probes<-merge(x=probes, y=agilent, by.x="probeID", by.y="efg_agilent_wholegenome_4x44k_v1", all.x=T)
colnames(probes)[3]<-"Agilent2EG"
probes<-merge(x=probes, y=refseq, by.x="Region", by.y="accession", all.x=T)
colnames(probes)[4]<-"Refseq2EG"
probes$Refseq2EG<-as.integer(probes$Refseq2EG)
probes<-merge(x=probes, y=ensembl, by.x="Region", by.y="trans_id", all.x=T)
colnames(probes)[5]<-"Ensembl2EG"
probes$Ensembl2EG<-as.integer(probes$Ensembl2EG)
probes<-merge(x=probes, y=genbank, by.x="Region", by.y="accession", all.x=T)
colnames(probes)[6]<-"Accession2EG"
probes$Accession2EG<-as.integer(probes$Accession2EG)

##How many probes do not have any gene information
probeNA<-apply(probes[,3:6], 1, function(x) sum(!is.na(x)) )
##How many probes have more than one entrez gene id associated with it.
probeMultis<-apply(probes[,3:6], 1, function(x) length(unique(x[!is.na(x)]))<=1 )

symbol<-as.data.frame(org.Hs.egSYMBOL2EG)
alias<-as.data.frame(org.Hs.egALIAS2EG)
colnames(alias)[2]<-"symbol"
symbol<-rbind(symbol, alias)
symbol<-tapply(symbol[,2], INDEX = symbol[,1], function(x) paste(x, collapse = ","))

probes2<-melt(probes)
probes2<-probes2[!is.na(probes2$value),]
eg2eg<-melt(probes[,3:6],id.vars = 1)
eg2eg<-eg2eg[!is.na(eg2eg$Agilent2EG) & !is.na(eg2eg$value),c(1,3)]
eg2eg<-eg2eg[eg2eg$Agilent2EG!=eg2eg$value,]
colnames(eg2eg)<-c("parent","child")
eg2eg<-data.frame(parent=c(eg2eg[,1],eg2eg[,2]), child=c(eg2eg[,2],eg2eg[,1]))


geneTable<-


scRNA<-read.table("../data/RefSeq_annotated_oTSS_Veh.bed", header=F, sep="\t", stringsAsFactors = F)
##Merge
scRNA<-merge(scRNA, eg.id, by.x="V4", by.y="trans_id", all.x=T)

##Merge
allgenes<-lapply(allgenes, function(x) merge(x,scRNA, all.x=T))
hits<-lapply(allgenes, function(x) x[(!is.na(x[,"gene_id"]) & x[,7]!=0 & x[,"Range"]<2),])
####I proceed under protest, gene annotation should not be this difficult


hits.table<-table(gene.results[,2:3])
q.hits<-hits.table[1,1]+hits.table[3,3]-1
m.hits<-sum(rowSums(hits.table[c(1,3),]))
n.hits<-sum(rowSums(hits.table))-m.hits
k.hits<-sum(colSums(hits.table[,c(1,3)]), hits.table[2,1],hits.table[2,3])
    
phyper(q=q.hits,m=m.hits, n=n.hits,k=k.hits)

x<-gene.results[,1:3]
x<-cbind(gene.results[,1],as.factor(abs(as.numeric(as.character(x[,2])))),as.factor(abs(as.numeric(as.character(x[,3])))))
table(x[,3],x[,2])
fisher.test(matrix(c(180,118,269,31812), nrow=2), alternative="greater")$p.value

agree<-hits.table[1,1]+hits.table[3,3]
disagree<-sum(hits.table[1,2],
              hits.table[1,3],
              hits.table[2,1],
              hits.table[2,3],
              hits.table[3,2],
              hits.table[3,1])
concordance<-(agree-disagree)/sum(agree,disagree)

################################################
### Gene List table output
################################################

for (i in 1:length(hits)){
    
    ##Write excel-stype tables of all expression data and DEGs-only data.
    write.table(hits[[i]][,c(1,2,3,4,6,9,7,8)] ,file=paste0("../results/",getDate(),"/",names(hits)[i],"_DEGs.txt"),
                quote=F, sep="\t", col.names=colnames(hits[[i]])[c(1,2,3,4,6,9,7,8)], row.names=F)
    
    write.table(allgenes[[i]][,c(1,2,3,4,6,9,7,8)] ,file=paste0("../results/",getDate(),"/",names(hits)[i],"_expression.txt"),
                quote=F, sep="\t", col.names=colnames(allgenes[[i]])[c(1,2,3,4,6,9,7,8)], row.names=F)
    
    ## Write BED format entries for all data and DEGs only. the Score field contains the log2 FC.
    write.table(hits[[i]][!is.na(hits[[i]][,"V4"]),c(14,15,16,13,4,18)] ,file=paste0("../results/",getDate(),"/",names(hits)[i],"_expression.bed"),
                quote=F, sep="\t", row.names=F, col.names=F)
    
#     write.table(hits[[i]][abs(hits[[i]][,5]==1),c(12,13,14,11,3,16)],
#                 file=paste0("../results/",getDate(),"/",names(hits)[i],"_DEG_expression.txt"), 
#                 quote=F,sep="\t", row.names=F, col.names=F)    
}

write.table(cbind(y0$genes$SystematicName, rowMeans(y0$E)), file=paste0("../results/",getDate(),"/Milliman_972_AvgExp.txt"), quote=F,sep="\t",row.names=F,col.names=F)

for (i in levels(as.factor(targets$groupName))){ 
    write.table(cbind(y0$genes$SystematicName,rowMeans(y0$E[,colnames(y0$E)%in%gsub(".txt","",targets[targets$groupName==i,3])])), file=paste0("../results/",getDate(),"/",i,"_totalExpression.txt"), row.names=F, col.names=F,sep="\t",quote=F)
    write.table(cbind(y0$genes$SystematicName,rowMeans(y0$E[,colnames(y0$E)%in%gsub(".txt","",targets[targets$groupName==i,3])])), file=paste0("../results/",getDate(),"/",i,"_totalExpression.bed"), row.names=F, col.names=F,sep="\t",quote=F)
    }





####Testing area####
