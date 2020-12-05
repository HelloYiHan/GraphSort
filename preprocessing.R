#! /usr/bin/Rscript
#Preprocess expression data for GraphSort
#Author: Han, Y. et al.
#R packages requirements: "KEGGgraph","KEGG.db","testit","SparseM","graph","funr","edgeR"
#Input expression file requirements: Tab-delimited with no quotations and no missing entries.
#                                    Ensembl gene IDs in column 1. Mixture labels in row 1.
#                                    Microarray data should should be quantile normalized and in non-log space.
#Parameters: input file, training file, input data type(rnaseq or microarray)

rm(list=ls())

#1.install and library packages---------------------------------------------
if(!"XML" %in% installed.packages()[,"Package"]) install.packages("XML", repos = "http://www.omegahat.net/R", quiet = T)

list.of.packages <- c("BiocManager","XML","testit", "funr")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages,quiet = T)

list.of.packages.bioconductor <- c("KEGGgraph","KEGG.db","SparseM","graph","edgeR","sva","preprocessCore")

new.packages <- list.of.packages.bioconductor[!(list.of.packages.bioconductor %in% installed.packages()[,"Package"])]

if(length(new.packages)) BiocManager::install(new.packages, quiet = T, update = F)

suppressPackageStartupMessages({
  sapply(c("KEGGgraph","KEGG.db","testit","SparseM","graph","funr","edgeR","sva","preprocessCore"), require, character.only=TRUE)
})


#2.parameters---------------------------------------------------------------------------
Args <- commandArgs(trailingOnly = TRUE)

if(is.na(Args[3])){stop("Not enough arguments!")}

path<-sys.script()

path<-paste(head(strsplit(x = path,split = "/")[[1]],-1),collapse = "/")


#3.functions----------------------------------------------------------------
if (Args[3]=="pancreatic") {
  arrange_format<-function(pfiles,graph_genes, bulk_rna_seq, gene_name_type){
    
    graph_genes<-graph_genes[order(graph_genes$Entrez),]
    
    graphs<- sapply(pfiles, parseKGML2Graph)
    
    panGraph <- mergeGraphs(graphs)
    
    nodes(panGraph)<-translateKEGGID2GeneID(nodes(panGraph))
    
    panSubGraph<-subGraph(as.character(graph_genes$Entrez), panGraph)
    
    pan_sparse<-graph2SparseM(panSubGraph)
    pan_coo<-as.matrix.coo(pan_sparse)
    
    numNode<-numNodes(panSubGraph)
    
    assert("num nodes is not 1818", numNode==1818)
    
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
    
    assert("rows of bulk_rna_seq are different from those of panSubGraph!",sum(bulk_rna_seq$Entrez==nodes(panSubGraph))==1818)
    
    row.names(bulk_rna_seq)<-bulk_rna_seq[,1]
    
    bulk_rna_seq<-bulk_rna_seq[,-c(1,2,3)]
    
    assert("row num of raw file should be 1818!",dim(bulk_rna_seq)[1]==1818)
    
    if(gene_name_type=='Entrez'){
      assert("row name of raw_file should be same with that of graph_genes",sum(row.names(bulk_rna_seq)==graph_genes[,1])==1818)
    } else if(gene_name_type=='Ensembl'){
      assert("row name of raw_file should be same with that of graph_genes",sum(row.names(bulk_rna_seq)==graph_genes[,2])==1818)
    } else if(gene_name_type=='Symbol'){
      assert("row name of raw_file should be same with that of graph_genes",sum(row.names(bulk_rna_seq)==graph_genes[,3])==1818)
    } else{
      stop("gene_name_type must be one of 'Entrez', 'Ensembl', 'Symbol'!")
    }
    
    simBulk<-as.data.frame(t(t(bulk_rna_seq)/colSums(bulk_rna_seq)*1e6))
    
    
    numGraph<-(dim(simBulk)[2])
    
    return.list<-list("simBulk"=simBulk, "pan_coo"=pan_coo, "numGraph"=numGraph, "numNode"=numNode)
    
    return(return.list)
  }
} else {
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
  
  #edge adjacent matrix
  AEdge<-cbind(pan_coo@ia,pan_coo@ja)
  for (Ngraph in c(1:(numGraph-1))) {
    AEdge<-rbind(AEdge,cbind(pan_coo@ia,pan_coo@ja)+Ngraph*numNode)
  }
  write.table(AEdge,paste(folder_path,"/",file_name,"/",file_name,"_A.txt",sep = ''),quote = F,sep = ',',row.names = F,col.names = F)
  
}


#4.pre-defined variables and data----------------------------------------------------------
#transform input file into required format for GraphSort
if(Args[3]=='pancreatic') {
  files<- c('hsa04950','hsa04972','hsa04974','hsa04610','hsa04512',
            'hsa04976','hsa04151','hsa05205','hsa04360','hsa04911',
            'hsa05215','hsa04668','hsa05230','hsa04510','hsa05323',
            'hsa04964','hsa04068','hsa01521','hsa04010','hsa05218',
            'hsa00430','hsa00051','hsa00140','hsa04971','hsa04024',
            'hsa05146','hsa00480','hsa05144','hsa04913','hsa04910',
            'hsa04922','hsa04931','hsa01522','hsa04930','hsa04940')
  pfile<- paste(path,"p_kgml", file, sep="/")
  
  graph_gene<-read.csv(paste(path,"pancreas_graph_genes.csv",sep = '/'),stringsAsFactors = F)
  
} else {
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
  
}


#5.read input----------------------------------------------------------------------------
input<-read.table(file = Args[1],sep = "\t",row.names = 1,header = T,stringsAsFactors = F)

#used to remove batch effect of input and training files
train<-read.table(file = Args[2],sep="\t",row.names = 1,header = T,stringsAsFactors = F)

n_train<-dim(train)[2]

batch_num<-c(rep(1,n_train),rep(2,dim(input)[2]))

#check if input file has required genes
assert("The input file has too little required genes.", sum(rownames(input) %in% rownames(train))>0)


#6.preprocess------------------------------------------------------------------------------------------
if (Args[3]=="rnaseq") {
  merge_input_train<-merge(train,input,by.x='row.names',by.y='row.names')

  rownames(merge_input_train)<-merge_input_train[,1]

  merge_input_train<-merge_input_train[,-1]
  
  #Reference of ComBat_seq, Zhang, Y., Parmigiani, G., & Johnson, W. E. (2020). ComBat-Seq: batch effect adjustment for RNA-Seq count data. bioRxiv, 904730.
  source(paste(path,"/ComBat_seq.R",sep=""))
  source(paste(path,"/helper_seq.R",sep=""))
  
  data_no_batch<-ComBat_seq(counts = as.matrix(merge_input_train),batch = batch_num)
  
  #arrange data to required format
  arrange_results<-arrange_format(pfile,graph_gene,data_no_batch[,(n_train+1):dim(merge_input_train)[2]],'Ensembl')
  
} else if (Args[3]=="microarray") {
  #input should be non-log, rownames should be gene symbol
  input_row_names<-rownames(input)
  input_column_names<-colnames(input)
  input<-as.data.frame(preprocessCore::normalize.quantiles(as.matrix(input)))
  colnames(input)<-input_column_names
  rownames(input)<-input_row_names
  
  merge_input_train<-merge(train,input,by.x='row.names',by.y='row.names')
  
  rownames(merge_input_train)<-merge_input_train[,1]
  
  merge_input_train<-merge_input_train[,-1]
  
  
  data_no_batch<-ComBat(dat=as.matrix(merge_input_train), batch = batch_num)
  
  arrange_results<-arrange_format(pfile,graph_gene,data_no_batch[,(n_train+1):dim(merge_input_train)[2]],'Symbol')
  
} else if (Args[3]=="pancreatic") {
  
  merge_input_train<-merge(train,input,by.x=0,by.y=0)
  
  rownames(merge_input_train)<-merge_input_train$Row.names
  
  merge_input_train<-merge_input_train[,-1]
  
  #Reference of ComBat_seq, Zhang, Y., Parmigiani, G., & Johnson, W. E. (2020). ComBat-Seq: batch effect adjustment for RNA-Seq count data. bioRxiv, 904730.
  source(paste(path,"/ComBat_seq.R",sep=""))
  source(paste(path,"/helper_seq.R",sep=""))
  
  data_no_batch<-ComBat_seq(counts = as.matrix(merge_input_train),batch = batch_num)
  
  #arrange data to required format
  arrange_results<-arrange_format(pfile,graph_gene,data_no_batch[,(n_train+1):dim(merge_input_train)[2]],'Symbol')
  
}


#7.save file
assert("The input file has too little required genes.", sum(colSums(arrange_results$simBulk)>0)==dim(arrange_results$simBulk)[2])

input_name<-tail(strsplit(x = Args[1],split = "/")[[1]],1)

dir_name<-paste("InputFile",input_name,sep = "_")

dir.create(dir_name)

write_file_no_pro(arrange_results$simBulk, arrange_results$pan_coo, arrange_results$numGraph, arrange_results$numNode,".", dir_name)

write.table(colnames(input),paste(".",dir_name,"input_file_sample_info.txt",sep = "/"),quote = F,sep = "\t",row.names = F,col.names = F)

system(paste("zip -r",paste0(dir_name,".zip"),dir_name,sep = " "))
