#Script to process qPCR d(RFU) dissociation data

library(reshape2)
library(ggplot2)
library(xlsx)

## Navigate to directory with qPCR files...
## setwd("~/Documents/qPCR_QC/2014-08-20/")


qpcr.table<-read.xlsx2(dir()[grep("Melt Curve Derivative Results.xlsx", dir())], 
                  sheetIndex = 1,as.data.frame = TRUE, colClasses = rep("numeric", 500))

qpcr.table<-qpcr.table[,2:length(colnames(qpcr.table))]

qt<-melt(qpcr.table, id = "Temperature", variable.name="Well", value.name="dRFU")

qt$row<-as.factor(substring(as.character(qt$Well), 0, 1))
qt$column<-as.factor(substring(as.character(qt$Well), 2))

primer<-read.xlsx2(dir(pattern="plate"), sheetIndex = 2,as.data.frame = TRUE, header=FALSE)

primer$X2<-gsub("/", "-", primer$X2)

for (i in levels(qt$column)) {
  
  png(paste("meltCurve_2step_MMTV_locus_", primer[as.numeric(i),2], ".png", sep =""), height=800, width=1200, units="px")
  p<-ggplot(qt[qt$column==i,], aes(Temperature, dRFU, col=row)) + geom_line() + labs(title=primer[as.numeric(i),2])
  print(p)
  dev.off()
  
}
rm(p)
