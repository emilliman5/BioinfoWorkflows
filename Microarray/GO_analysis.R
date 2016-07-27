library(GO.db)
library(annotate)
library(topGO)
library(RColorBrewer)
library(xtable)
library(Rgraphviz)
library(org.Hs.eg.db)
library(biomaRt)

##remove genes without entrezIDs

# GoDegList<-lapply(hits, function(x) x[(!is.na(x[,"gene_id"]) & x[,7]!=0 & x[,"Range"]<2),c(1,7)])

# GoDegs<-lapply(GoDegList, function(x) x[,2])
# GoDegs<-lapply(GoDegList, function(x) names(x)<-x[,1]) 
# 
# GoDegs<-GoDegList[[1]][,2]
# eg.names<-GoDegList[[1]][,1]
# names(GoDegs)<-eg.names
##remove genes with no GO mapping

# GoDegs<-lapply(GoDegs, function(x) x[sapply(mget(names(x), org.Hs.egGO), 
#                function(y) {
#                    if (length(y) ==1 && is.na(y))
#                         FALSE
#                    else TRUE
#                 }),])

# GoDegList<-lapply(GoDegList, function(x) x[!duplicated(x[,1]),])
# GoDegList<-lapply(GoDegList, function(x) names(x)<-x[,1])

agilentUniverse<-hits[[1]][!is.na(hits[[1]][,"gene_id"]),c(1,7)]
entrezUniverse<-agilentUniverse[,2]
names(entrezUniverse)<-agilentUniverse[,1]

GOparams<- new("topGOdata",
               geneSel=function(p) p!=0,
               allGenes=entrezUniverse,
               ontology="BP",
               annot=annFUN.org,
               mapping= "org.Hs.eg.db",
               ID="Entrez")

resultsFisher<-runTest(GOparams, algorithm = "classic",statistic = "fisher")
GenTable(GOparams, classicFisher = resultsFisher, topNodes = 10)

showSigOfNodes(GOparams, score(resultsFisher), firstSigNodes = 10, useInfo = 'all')
