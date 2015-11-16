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

primer<-read.xlsx2(dir(pattern="plate"), sheetIndex = 1, as.data.frame = TRUE, header=FALSE)

primer.list<-unname(unlist(primer[2,2:25]))[1:24]
primer.list<-droplevels(primer.list)
primer.list<-gsub("/", "-", primer.list)

for (id in 1:length(primer.list)){  
  
  if(length(grep(pattern = "\\s$", x = primer.list[id]))>0){
    primer.list[id]<-substr(primer.list[id], 1, nchar(primer.list[id])-1) 
  }
}

for (i in levels(qt$column)) {
  
  png(paste("meltCurve_AmpCurve_", primer.list[as.numeric(i)], ".png", sep =""), height=1800, width=1200, units="px")
  p1<-ggplot(qt[qt$column==i,], aes(Temperature, dRFU, col=row)) + geom_line(size=1) + ggtitle(as.character(primer.list[as.numeric(i)]))
  p2<-ggplot(q.amp[q.amp$column==i,], aes(Cycle, RFU, col=row)) + geom_line(size=1) + ggtitle(as.character(primer.list[as.numeric(i)]))
  plot_list=list(p1, p2)
  do.call(grid.arrange, c(plot_list, list(nrow=2)))
  dev.off()
  
}

rm(p1, p2)


qPCR[[l]]<-list(substr(paste(getwd()), nchar(paste(getwd()))-9, nchar(paste(getwd()))), qt, q.amp, primer.list)
l<-l+1
