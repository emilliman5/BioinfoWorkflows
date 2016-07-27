library(reshape2)
library(igraph)
library(slam)
library(pheatmap)
library(amap)
library(RColorBrewer)

files<-dir(pattern=".txt")

annot<-lapply(files, function(x) read.table(x, header=T,quote="", 
                   allowEscapes=F, stringsAsFactors=F, sep="\t", fill=T))


names(annot)<-gsub("_motifAnnotations.txt", "", files)

for(i in seq_along(annot)){
    colnames(annot[[i]])<-gsub("\\.Distance\\.From.+","",colnames(annot[[i]]), perl=T)
    colnames(annot[[i]])[1]<-"PeakID"
}

annotmotifs<-lapply(annot, function(x) x[,-c(1:21)])
annotmotifs<-lapply(annotmotifs, function(x) {
                                apply(x,1, gsub, 
                                        pattern="\\([ACTG]+,[\\+-],\\d\\.\\d\\d\\)", replacement="", perl=T)
                })
annotmotifs<-lapply(annotmotifs, t)
annotmotifs<-lapply(annotmotifs, function(x) ifelse(x=="", 0,ifelse(is.na(x)==T, 0, 1)))

for(i in seq_along(annotmotifs)){
    class(annotmotifs[[i]])<-"numeric"
    annotmotifs[[i]]<-annotmotifs[[i]][,which(!is.na(col_sums(annotmotifs[[i]])))]
}
matrixD<-lapply( annotmotifs, function(x) {
    crossprod_simple_triplet_matrix(as.simple_triplet_matrix(!x))
})
# annotmotifs<-as.simple_triplet_matrix(annotmotifs)

annotcoOccur<-lapply( annotmotifs, function(x){
    crossprod_simple_triplet_matrix(as.simple_triplet_matrix(x)) })

annotJaccard<-lapply(seq_along(annot), function(x) {
    annotcoOccur[[x]]/(nrow(annot[[x]])-matrixD[[x]])
})

motifJaccardclust<-lapply(annotJaccard, hcluster,  method="euclid",nbproc = 8)
motifAdjMatClust<-lapply( annotcoOccur,  hcluster, method="euclid",nbproc = 8)

lapply(seq_along(annotcoOccur), function(x) {
    png(paste0(names(annot)[x],"_minOdds_Motif_coOccurence.png"), height=2000, width=2500, units="px")
        pheatmap(sqrt(annotcoOccur[[x]]),color=colorRampPalette(brewer.pal(9,"Greens"))(100),
            cluster_rows = motifAdjMatClust[[x]], cluster_cols = motifAdjMatClust[[x]])
    dev.off()

    png(paste0(names(annot)[x],"_minOdds_Motif_Jaccard.png"), height=2000, width=2500, units="px")
        pheatmap(annotJaccard[[x]],color=colorRampPalette(brewer.pal(9,"Greens"))(100),
            cluster_rows = motifJaccardclust[[x]], cluster_cols = motifJaccardclust[[x]])
    dev.off()
})
###Remove Motifs that have perfect concordance...

edges<-lapply(annotJaccard, function(x){ 
    e<-melt(x)
    e[as.integer(e$Var2)>as.integer(e$Var1),]
    })

nodesToRemove<-lapply(seq_along(edges), function(x)
    rownames(annotJaccard[[x]]) %in% edges[[x]][edges[[x]]$value==1,2]
)
annotJaccard2<-mapply(x=annotJaccard, y=nodesToRemove,
                      function(x,y) x[!y,!y]
)
motifJaccardclust<-lapply(annotJaccard2, hcluster,method="euclid",nbproc = 8)
motifAdjMatClust<-lapply(annotJaccard2, hcluster, method="euclid",nbproc = 8)

lapply(seq_along(annotJaccard2), function(x) {
    png(paste0(names(annot)[x],"_minOdds_Motif_Jaccard2.png"), height=2000, width=2500, units="px")
    pheatmap(annotJaccard2[[x]],color=colorRampPalette(brewer.pal(9,"Greens"))(100),
             cluster_rows = motifJaccardclust[[x]], cluster_cols = motifJaccardclust[[x]])
    dev.off()
})


edgelist1<-lapply(annotJaccard2, function(x){
    e<-melt(x)
    e[as.integer(e$Var1)<as.integer(e$Var2),]
})

edgesIDs<-lapply(edgelist1, function(x)
    paste0(as.character(x$Var1),"_",as.character(x$Var2)))

edgesUnion<-(edgesIDs[[1]] %in% edgesIDs[[2]]) & (edgesIDs[[1]] %in% edgesIDs[[3]])

edgeDF<-lapply(edgesIDs, function(x) x %in% edgesIDs[[1]][edgesUnion])
edgesDF<-mapply(x=edgelist1, y=edgeDF, function(x,y) x[y,], SIMPLIFY = F)

edgeDF<-do.call(cbind, edgesDF)
edgeDF<-edgeDF[,-c(4,5,7,8)]
colnames(edgeDF)<-c("to", "from",names(annot))

kclust<-kmeans(edgeDF[,3:5],centers = 8)
ids<-lapply(1:8, function(x) sample(which(kclust$cluster==x),100))
ids<-unlist(ids)

png("Motif_agg1.png", height=3000, width=1000, units="px")
pheatmap(edgeDF[ids,3:5], annotation_row = edgeDF[ids,1:2],border_color = NA,
         cluster_rows = F, cluster_cols = F,
         color=colorRampPalette(brewer.pal(9,"Greens"))(100))
dev.off()

png("Motif_agg2.png", height=3000, width=1000, units="px")
    pheatmap(edgeDF[,3:5], annotation_row = edgeDF[,1:2],k_means=kclust,
             cluster_rows = F, cluster_cols = F,border_color=NA,
             color=colorRampPalette(brewer.pal(9,"Greens"))(100))
dev.off()

edgelist<-lapply(edgelist1, function(x)
    x[x$value>quantile(x$value,0.9),])

for(i in seq_along(edgelist)){
    colnames(edgelist[[i]])<-c("from","to","Weight")
}

graph<-lapply(edgelist, graph_from_data_frame ,directed=FALSE)

nodeSize<-lapply(seq_along(annotJaccard2), function(x){
    col_sums(annotmotifs[[x]])[colnames(annotmotifs[[x]]) %in% V(graph[[x]])$name]
})

for(i in seq_along(graph)){
    V(graph[[i]])$Size<-nodeSize[[i]]
}

lapply(seq_along(graph), function(x) {
    write_graph(graph[[x]], file=paste0(names(annotmotifs)[x],"_network.graphml")
                ,format="graphml")
    })
graphLay<-layout_with_fr(graph,weights = 1-edgelist$Weight)
graphLay<-layout_with_kk(graph,weights = -log(edgelist$Weight))
graphLay<-layout_with_dh(graph)

png("motif_network.png", height=2500, width=2500, units="px")
plot(graph, vertex.size=2)
dev.off()

png("motif_network_v7.png", height=2500, width=2500, units="px")
plot(graph, vertex.size=2, layout=graphLay)
dev.off()
