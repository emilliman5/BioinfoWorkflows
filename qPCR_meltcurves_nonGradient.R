#Script to process qPCR d(RFU) dissociation data

library(reshape2)
library(ggplot2)
library(xlsx)
library(grid)
library(gridExtra)

## Navigate to directory with qPCR files...
## setwd("~/Documents/qPCR_QC/2014-08-20/")


qpcr.table<-read.xlsx2(dir()[grep("Melt Curve Derivative Results.xlsx", dir())], 
                  sheetIndex = 1,as.data.frame = TRUE, colClasses = rep("numeric", 500))

qpcr.amp<-read.xlsx2(dir()[grep("Quantification Amplification Results", dir())], 
                     sheetIndex = 1,as.data.frame = TRUE, colClasses = rep("numeric", 500))


qpcr.amp<-qpcr.amp[,2:length(colnames(qpcr.amp))]
qpcr.table<-qpcr.table[,2:length(colnames(qpcr.table))]

q.amp<-melt(qpcr.amp, id = "Cycle", variable.name="Well", value.name="RFU")
qt<-melt(qpcr.table, id = "Temperature", variable.name="Well", value.name="dRFU")

qt$row<-as.factor(substring(as.character(qt$Well), 0, 1))
qt$column<-as.factor(substring(as.character(qt$Well), 2))

q.amp$row<-as.factor(substring(as.character(q.amp$Well), 0, 1))
q.amp$column<-as.factor(substring(as.character(q.amp$Well), 2))

##Adding well annotations
q.amp[q.amp$row %in% c("A","B","E","F","I","J"),6]<-"NTC"
q.amp[is.na(q.amp$V6),6]<-"DNA"

qt[qt$row %in% c("A","B","E","F","I","J"),6]<-"NTC"
qt[is.na(qt$V6),6]<-"DNA"

q.amp[q.amp$row %in% c("A","B","C","D"),7]<-"Platinum.Taq"
q.amp[q.amp$row %in% c("E","F","G","H"),7]<-"HSTaq"
q.amp[q.amp$row %in% c("I","J","K","L"),7]<-"Taq"

qt[qt$row %in% c("A","B","C","D"),7]<-"Platinum.Taq"
qt[qt$row %in% c("E","F","G","H"),7]<-"HSTaq"
qt[qt$row %in% c("I","J","K","L"),7]<-"Taq"

q.amp[q.amp$column %in% c("1","2"),8]<-"DSB-70and71"
qt[qt$column %in% c("1","2"),8]<-"DSB-70and71"
q.amp[q.amp$column %in% c("3","4"),8]<-"Fam105A"
qt[qt$column %in% c("3","4"),8]<-"Fam105A"
q.amp[q.amp$column %in% c("5","6"),8]<-"GPDgen"
qt[qt$column %in% c("5","6"),8]<-"GPDgen"
q.amp[q.amp$column %in% c("7","8"),8]<-"H1c"
qt[qt$column %in% c("7","8"),8]<-"H1c"
q.amp[q.amp$column %in% c("9","10"),8]<-"H1d"
qt[qt$column %in% c("9","10"),8]<-"H1d"
q.amp[q.amp$column %in% c("11","12"),8]<-"H10"
qt[qt$column %in% c("11","12"),8]<-"H10"
q.amp[q.amp$column %in% c("13","14"),8]<-"NucA"
qt[qt$column %in% c("13","14"),8]<-"NucA"
q.amp[q.amp$column %in% c("15","16"),8]<-"NucB"
qt[qt$column %in% c("15","16"),8]<-"NucB"
q.amp[q.amp$column %in% c("17","18"),8]<-"NucF"
qt[qt$column %in% c("17","18"),8]<-"NucF"
q.amp[q.amp$column %in% c("19","20"),8]<-"NeoR"
qt[qt$column %in% c("19","20"),8]<-"NeoR"
q.amp[q.amp$column %in% c("21","22"),8]<-"PNMT-300"
qt[qt$column %in% c("21","22"),8]<-"PNMT-300"
q.amp[q.amp$column %in% c("23","24"),8]<-"SNAI1_GBS"
qt[qt$column %in% c("23","24"),8]<-"SNAI1_GBS"



##Make Plots
for (id in 1:length(primer.list)){  
  
  if(length(grep(pattern = "\\s$", x = primer.list[id]))>0){
    primer.list[id]<-substr(primer.list[id], 1, nchar(primer.list[id])-1) 
  }
}

ggplot(data = qt, aes(Temperature, dRFU, col=c(V6, V7, V8)))+geom_line(size=1)


for (i in levels(qt$V8)) {
  
  png(paste("meltCurve_AmpCurve_", i, ".png", sep =""), height=1800, width=1200, units="px")
  p1<-ggplot(data = qt[qt$V8==i,], aes(Temperature, dRFU, group=Well, col=V9))+geom_line(size=1) + ggtitle(as.character(i)) + theme(legend.text=element_text(size=16), plot.title=element_text(size=16))
  p2<-ggplot(data = q.amp[q.amp$V8==i,], aes(Cycle, RFU, group=Well, col=V9)) + geom_line(size=1) + ggtitle(as.character(i)) + theme(legend.text=element_text(size=16),plot.title=element_text(size=16))
  plot_list=list(p1, p2)
  do.call(grid.arrange, c(plot_list, list(nrow=2)))
  dev.off()
  
}

rm(p1, p2)


qPCR[[l]]<-list(substr(paste(getwd()), nchar(paste(getwd()))-9, nchar(paste(getwd()))), qt, q.amp, primer.list)
l<-l+1
