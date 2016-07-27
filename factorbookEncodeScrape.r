library(XML)
library(RCurl)
library(parallel)

dir.create("FactorBookMotifs")
fbFwdURL<-"http://compbio.mit.edu/encode-motifs/logos/table/logos/mat/fwd/"
fbRevURL<-"http://compbio.mit.edu/encode-motifs/logos/table/logos/mat/rev/"

fbFwd<-getURL("http://compbio.mit.edu/encode-motifs/logos/table/logos/mat/fwd/")
fbRev<-getURL("http://compbio.mit.edu/encode-motifs/logos/table/logos/mat/rev/")

fbFwd<-htmlParse(fbFwd)
fbRev<-htmlParse(fbRev)

matsFwd<-xpathSApply(fbFwd ,"//a/@href") 
matsRev<-xpathSApply(fbFwd ,"//a/@href") 

matsRev<-grep(".txt",matsRev,value=T)
matsFwd<-grep(".txt",matsFwd,value=T)

mclapply(matsFwd, mc.cores = 4, 
                    function(x) download.file(paste0(fbFwdURL,x),quiet = T,
                                              destfile =paste0("FactorBookMotifs/fwd_",x)))

mclapply(matsFwd, mc.cores = 4, 
         function(x) download.file(paste0(fbRevURL,x),quiet = T,
                                   destfile =paste0("FactorBookMotifs/rev_",x)))
