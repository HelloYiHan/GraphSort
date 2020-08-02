#GraphSort preprocessing script
#Author: Yi Han et al.

#R packages requirements: "KEGGgraph","KEGG.db","testit","SparseM","graph","funr"

#Input expression file requirements: Tab-delimited with no quotations and no missing entries.
#                                    Ensembl gene IDs in column 1. Mixture labels in row 1.
#                                    Microarray data should should be quantile normalized and in non-log space.

rm(list=ls())

#install and library required packages
install.packages("BiocManager")
BiocManager::install(c("KEGGgraph","KEGG.db","SparseM","graph"))
list.of.packages <- c("testit", "funr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
sapply(c("KEGGgraph","KEGG.db","testit","SparseM","graph","funr"), require, character.only=TRUE)



Args <- commandArgs(trailingOnly = TRUE)

arrange_format<-function(pfiles,graph_genes, bulk_rna_seq, gene_name_type){
  
  graph_genes<-graph_genes[order(graph_genes$Entrez),]
  
  graphs<- sapply(pfiles, parseKGML2Graph)
  
  panGraph <- mergeGraphs(graphs)
  
  nodes(panGraph)<-translateKEGGID2GeneID(nodes(panGraph))
  
  panSubGraph<-subGraph(as.character(graph_genes$Entrez), panGraph)
  
  pan_sparse<-graph2SparseM(panSubGraph)
  pan_coo<-as.matrix.coo(pan_sparse)
  
  numNode<-numNodes(panSubGraph)
  
  assert("num nodes is not 2504", numNode==2504)
  
  assert("row names of bulk_rna_seq should not have duplicated items!",sum(duplicated(row.names(bulk_rna_seq)))==0)
  
  assert("row names of bulk_rna_seq should not have NA!",sum(is.na(row.names(bulk_rna_seq)))==0)
  
  if(gene_name_type=='Entrez'){
    bulk_rna_seq<-merge(graph_genes,bulk_rna_seq,by.x=1, by.y='row.names', all.x=T)
  } else if(gene_name_type=='Ensembl'){
    bulk_rna_seq<-merge(graph_genes,bulk_rna_seq,by.x=2, by.y='row.names', all.x=T)
  } else if(gene_name_type=='Symbol'){
    bulk_rna_seq<-merge(graph_genes,bulk_rna_seq,by.x=3, by.y='row.names', all.x=T)
  } else{
    stop("gene_name_type must be one of 'Entrez', 'Ensembl', 'Symbol'!")
  }
  
  bulk_rna_seq[is.na(bulk_rna_seq)]<-0
  
  bulk_rna_seq<-bulk_rna_seq[order(bulk_rna_seq$Entrez),]
  
  assert("rows of bulk_rna_seq are different from those of panSubGraph!",sum(bulk_rna_seq$Entrez==nodes(panSubGraph))==2504)
  
  row.names(bulk_rna_seq)<-bulk_rna_seq[,1]
  
  bulk_rna_seq<-bulk_rna_seq[,-c(1,2,3)]
  
  assert("row num of raw file should be 2504!",dim(bulk_rna_seq)[1]==2504)
  
  if(gene_name_type=='Entrez'){
    assert("row name of raw_file should be same with that of graph_genes",sum(row.names(bulk_rna_seq)==graph_genes[,1])==2504)
  } else if(gene_name_type=='Ensembl'){
    assert("row name of raw_file should be same with that of graph_genes",sum(row.names(bulk_rna_seq)==graph_genes[,2])==2504)
  } else if(gene_name_type=='Symbol'){
    assert("row name of raw_file should be same with that of graph_genes",sum(row.names(bulk_rna_seq)==graph_genes[,3])==2504)
  } else{
    stop("gene_name_type must be one of 'Entrez', 'Ensembl', 'Symbol'!")
  }
  
  simBulk<-as.data.frame(t(t(bulk_rna_seq)/colSums(bulk_rna_seq)*1e6))
  
  
  numGraph<-(dim(simBulk)[2])
  
  return.list<-list("simBulk"=simBulk, "pan_coo"=pan_coo, "numGraph"=numGraph, "numNode"=numNode)
  
  return(return.list)
}

write_file_no_pro<-function(bulk, pan_coo, numGraph, numNode, folder_path, file_name){
  
  #graph indicator
  graphIndicator<-c()
  for (Ngraph in c(1:numGraph)) {
    graphIndicator<-c(graphIndicator,rep(Ngraph,numNode))
  }
  write.table(graphIndicator,paste(folder_path,"/",file_name,"/",file_name,"_graph_indicator.txt",sep = ''),quote = F,sep = ',',row.names = F,col.names = F)
  
  #node attributes
  nodeAttribute<-c()
  for (Ngraph in c(1:(numGraph))) {
    nodeAttribute<-c(nodeAttribute,bulk[,Ngraph])
  }
  write.table(nodeAttribute,paste(folder_path,"/",file_name,"/",file_name,"_node_attributes.txt",sep = ''),quote = F,sep = ',',row.names = F,col.names = F)
  
  #A (edge adjacent matrix)
  AEdge<-cbind(pan_coo@ia,pan_coo@ja)
  for (Ngraph in c(1:(numGraph-1))) {
    AEdge<-rbind(AEdge,cbind(pan_coo@ia,pan_coo@ja)+Ngraph*numNode)
  }
  write.table(AEdge,paste(folder_path,"/",file_name,"/",file_name,"_A.txt",sep = ''),quote = F,sep = ',',row.names = F,col.names = F)
  
}




#parameters
if(!is.na(Args[2])){stop("More than one argument!")}

input_name<-Args[1]

path<-sys.script()

path<-paste(head(strsplit(x = path,split = "/")[[1]],-1),collapse = "/")

#read input
input<-read.table(file = paste(path,input_name,sep = "/"),sep = "\t",row.names = 1,header = T)

#used to remove batch effect of input and training files
train<-read.table(file = paste(path,"/data_for_rem_bat_eff.txt",sep=""),sep="\t",row.names = 1,header = F)

#check if input file has required genes
assert("The input file has too little required genes.", sum(rownames(input) %in% rownames(train))>0)

merge_input_train<-merge(train,input,by.x='row.names',by.y='row.names')

rownames(merge_input_train)<-merge_input_train[,1]

merge_input_train<-merge_input_train[,-1]

batch_num<-c(rep(1,2000),rep(2,dim(input)[2]))

#Reference of ComBat_seq
#Zhang, Y., Parmigiani, G., & Johnson, W. E. (2020). ComBat-Seq: batch effect adjustment for RNA-Seq count data. bioRxiv, 904730.
source(paste(path,"/ComBat_seq.R",sep=""))
source(paste(path,"/helper_seq.R",sep=""))

data_no_batch<-ComBat_seq(counts = as.matrix(merge_input_train),batch = batch_num)

#transform input file into required format for GraphSort
file<- c("hsa04010.xml","hsa04620.xml","hsa04621.xml" ,"hsa05133.xml","hsa04218.xml",
          "hsa04060.xml","hsa04061.xml","hsa04062.xml","hsa04064.xml","hsa05340.xml",
          "hsa04110.xml","hsa04115.xml","hsa04141.xml","hsa04142.xml","hsa05332.xml",
          "hsa04145.xml","hsa04151.xml","hsa04210.xml","hsa04217.xml","hsa05330.xml",
          "hsa04380.xml","hsa04514.xml","hsa04612.xml","hsa05323.xml","hsa05321.xml",
          "hsa04630.xml","hsa04640.xml","hsa05132.xml","hsa05320.xml","hsa05310.xml",
          "hsa04650.xml","hsa04657.xml","hsa04658.xml","hsa05152.xml","hsa04672.xml",
          "hsa04659.xml","hsa04660.xml","hsa04668.xml","hsa04940.xml","hsa04625.xml",
          "hsa04662.xml","hsa04664.xml","hsa04666.xml","hsa05140.xml","hsa05144.xml",
          "hsa05146.xml","hsa05160.xml","hsa05161.xml","hsa05163.xml","hsa05202.xml",
          "hsa05166.xml","hsa05167.xml","hsa05169.xml","hsa05170.xml","hsa05142.xml",
          "hsa05205.xml","hsa05219.xml","hsa05221.xml","hsa05235.xml","hsa05164.xml",
          "hsa05416.xml","hsa05134.xml","hsa05145.xml","hsa05150.xml","hsa05162.xml")

pfile<- paste(path,"kgml", file, sep="/")

graph_gene<-read.csv(paste(path,"GraphGenes.csv",sep = '/'),stringsAsFactors = F)

#arrange data to required format
arrange_results<-arrange_format(pfile,graph_gene,data_no_batch[,2001:dim(merge_input_train)[2]],'Ensembl')

assert("The input file has too little required genes.", sum(colSums(arrange_results$simBulk)>0)==dim(arrange_results$simBulk))

write_file_no_pro(arrange_results$simBulk, arrange_results$pan_coo, arrange_results$numGraph, arrange_results$numNode,path, 'InputFile')
