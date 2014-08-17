#Script to process qPCR d(RFU) dissociation data

library(reshape2)
library(ggplot2)
library(xlsx)

setwd("~/Documents/qPCR_QC/2014-08-15/")
qpcr.table<-read.xlsx2("20140815_174558_CT009851_PRIMER_VALID -  Melt Curve Derivative Results.xlsx", 
                  sheetIndex = 1,as.data.frame = TRUE, colClasses = rep("numeric", 500))

qpcr.table<-qpcr.table[,2:length(colnames(qpcr.table))]

qt<-melt(qpcr.table, id = "Temperature", variable.name="Well", value.name="dRFU")

qt$row<-as.factor(substring(as.character(qt$Well), 0, 1))
qt$column<-as.factor(substring(as.character(qt$Well), 2))

primer<-c("NotI 70-71","NotI 84-85", "NotI 88-89", "NotI 90-91","Not-4","H10", "GPD gen","FAM1095A", "SNAI1 GBS", "PNMT -300", "NFK1a 1800", "GILZ plus150", "NucA-Luc", "NucA","NucA-2", "NucB", "NucF", "", "", "", "", "", "", "")

for (i in levels(qt$column)) {
  
  png(paste("meltCurve_2step_", primer[as.numeric(i)], ".png", sep =""), height=800, width=1200, units="px")
  p<-ggplot(qt[qt$column==i,], aes(Temperature, dRFU, col=row)) + geom_line() + labs(title=primer[as.numeric(i)])
  print(p)
  dev.off()
  
}

