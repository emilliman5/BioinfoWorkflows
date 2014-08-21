#Script to process qPCR d(RFU) dissociation data

library(reshape2)
library(ggplot2)
library(xlsx)

setwd("~/Documents/qPCR_QC/2014-08-20/")
qpcr.table<-read.xlsx2("20140820_133749_CT009851_PRIMER_VALID -  Melt Curve Derivative Results.xlsx", 
                  sheetIndex = 1,as.data.frame = TRUE, colClasses = rep("numeric", 500))

qpcr.table<-qpcr.table[,2:length(colnames(qpcr.table))]

qt<-melt(qpcr.table, id = "Temperature", variable.name="Well", value.name="dRFU")

qt$row<-as.factor(substring(as.character(qt$Well), 0, 1))
qt$column<-as.factor(substring(as.character(qt$Well), 2))

primer<-c("MMTV 142-143","MMTV 144-145", "MMTV 146-147", "MMTV 148-149","MMTV 150-151", 
          "MMTV 152-153","MMTV 154-155", "MMTV 156-157", "MMTV 158-159", "MMTV 160-161", "TPA-25 Dilute ChIP", "TPA-25 ChIP", "","", "", "", "", "", "", "", "", "", "")

for (i in levels(qt$column)) {
  
  png(paste("meltCurve_2step_MMTV_locus_", primer[as.numeric(i)], ".png", sep =""), height=800, width=1200, units="px")
  p<-ggplot(qt[qt$column==i,], aes(Temperature, dRFU, col=row)) + geom_line() + labs(title=primer[as.numeric(i)])
  print(p)
  dev.off()
  
}

