library(gdata)
library(docopt)

"Usage: qPCRFormat.R <files>..." -> doc

myopts<-docopt(doc)
print(myopts)


outfiles<-gsub(".xls","", myopts$files)
samples<-paste0("S", 1:36)
targets<-c("RPL13A","GAPDH","B2M","FAM90A1","H1e")

data<-lapply(myopts$files, read.xls, sheet=1)
data1<-lapply(data, function(x) x[c(2,8)])

data2<-lapply(data1, function(x){
    dat1<-data.frame(rep1=numeric(),rep2=numeric(),rep3=numeric(),rep4=numeric())
    for(row in seq(0,383,by=48)){
        for(column in seq(0,23,by=2)){
            dat1<-rbind(dat1, x[c(row+column+1,row+column+2,row+column+25, row+column+26),2])      
        }
    }
    dat1
})

data2<-lapply(data2, setNames, nm=c("Rep1","Rep2","Rep3","Rep4"))
data2<-mapply(x=data2, function(x,y) write.csv(x, file=paste0(y,".txt")), y=outfiles)
