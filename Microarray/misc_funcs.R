getExp<-function(gene="AdamTS8",expt=1)
{
    id<-unique(eg.id[grep(gene, tmp[,3], ignore.case = T),1])
    head(hits[[expt]][hits[[expt]][,1]==as.character(id),])[1,]
}


exp<-as.data.frame(cbind(as.character(hits[[1]][,8]), hits[[1]][,7],do.call(cbind, lapply(hits, function(x) x[,10]))))

exp[,-c(1,2)]<-as.numeric(exp[,-c(1,2)])

exp<-exp[(abs(exp[,3])>1|abs(exp[,4])>1|abs(exp[,5])>1|abs(exp[,6])>1|abs(exp[,7])>1|abs(exp[,8])>1),]
    
    

write.table(as.data.frame(cbind(as.character(hits[[1]][,8]), hits[[1]][,7],
                                do.call(cbind, lapply(hits, function(x) x[,10])))),
            "../results/2015jun08/geneExpTable.txt",quote=F,sep = "\t")
