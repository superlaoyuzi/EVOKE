# ---- library packages ----
#devtools::install_github("Japrin/STARTRAC")
#devtools::install_github("neurorestore/Augur")
lapply(c("Startrac","ggplot2","tictoc","ggpubr","ComplexHeatmap","RColorBrewer","ggsci",
         "circlize","reshape2","clusterProfiler","tidyverse","sscVis","Seurat","readr","qs","BiocParallel","Augur"), 
       library, character.only = T)
setwd("/home/data/t070404/tumor_progression")
# ---- Load Seurat object ----
file<-list.files("./stage1")
file<-file[grep("rds",file)]
# ---- Process the data circularly ----
for(i in file){
  patientid<-str_extract(i, "[A-Z]+_P\\d+")
  SO<-readRDS(paste("./stage1/",i,sep = ""))
  ##Proportion of cell types
  Idents(SO) <- SO$celltype
  Cellratio <- prop.table(table(SO$celltype, SO$sample_type), margin = 2)#计算各组样本不同细胞群比例
  write.table(Cellratio,paste("./patient_cell_prop/",patientid,"_cell_prop.txt",sep = ""),quote = F,row.names = T,sep = "\t")
  ##Ro/e Cell tissue preference analysis(Based on Startrac)
  metadata <- SO@meta.data
  metadata <- metadata[,c('patient','sample_type','celltype')]
  Roe <- calTissueDist(metadata,
                       byPatient = F,
                       colname.cluster = "celltype",# Different subpopulations of cells
                       colname.patient = "patient", # Different samples
                       colname.tissue = "sample_type",# different organizations
                       method = "chisq", # "chisq", "fisher", and "freq"
                       min.rowSum = 0) 
  if(length(table(SO$sample_type))==1){
    Roe<-data.frame(Roe,row.names = unique(SO$celltype))
    colnames(Roe)<-unique(SO$sample_type)
  }
  write.table(Roe,paste("./patient_Roe/",patientid,"__Roe_celltype_sampletype.txt",sep = ""),quote = F,row.names = T,sep = "\t")
##Determine cellular differential disturbances(Based on Augur)
  if(length(table(SO$sample_type))>1){
    augur <- calculate_auc(SO,
                           cell_type_col = "celltype", # cell subpopulation data
                           label_col = "sample_type", # Experimental grouping data
                           n_threads = 8)
    qsave(augur,paste("./patient_augur/",patientid,"_augur.qs",sep = ""))
    write.table(augur$AUC,paste("./patient_augur/",patientid,"_augur_auc.txt",sep = ""),quote = F,row.names = F,sep = "\t")
    
  }
  ##Differentially expressed genes in cells
  if("Health" %in% unique(SO$sample_type)){
    markers_in_health<-FindAllMarkers(SO[,which(SO$sample_type == "Health")],only.pos = T) 
    markers_in_health$sample_type<-"Health"
  }else{markers_in_health<-c()}
  if("Inflamed" %in% unique(SO$sample_type)){
    markers_in_precancerous<-FindAllMarkers(SO[,which(SO$sample_type == "Inflamed")],only.pos = T) 
    markers_in_precancerous$sample_type<-"Inflamed"
  }else{markers_in_precancerous<-c()}
  if("Margin" %in% unique(SO$sample_type)){
    markers_in_precancerous<-FindAllMarkers(SO[,which(SO$sample_type == "Margin")],only.pos = T) 
    markers_in_precancerous$sample_type<-"Margin"
  }else{markers_in_precancerous<-c()}
  if("Primary" %in% unique(SO$sample_type)){
    markers_in_primary<-FindAllMarkers(SO[,which(SO$sample_type == "Primary")],only.pos = T) 
    markers_in_primary$sample_type<-"Primary"
  }else{markers_in_primary<-c()}
  if("Metastasis" %in% unique(SO$sample_type)){
    markers_in_metastasis<-FindAllMarkers(SO[,which(SO$sample_type == "Metastasis")],only.pos = T) 
    markers_in_metastasis$sample_type<-"Metastasis"
  }else{markers_in_metastasis<-c()}

  
  markers<-rbind(markers_in_health,markers_in_precancerous,markers_in_primary,markers_in_metastasis)
  markers$patient<-patientid
  write.table(markers,paste("./patient_cellmarker/",patientid,"_markers.txt",sep = ""),quote = F,row.names = F,sep = "\t")
  
  top5<-markers %>% group_by(sample_type,cluster) %>% top_n(10,avg_log2FC)
  
  cell_gene<-aggregate(top5$gene,list(top5$sample_type,top5$cluster),paste,collapse = ";")
  cell_gene<-spread(cell_gene,Group.1,x)
  cn<-paste(colnames(cell_gene),"_gene_signature",sep  = "")
  colnames(cell_gene)<-c("celltype",cn[2:length(cn)])
  ##GO enrichment analysis
  genetable<-aggregate(markers$gene,list("sample_type" = markers$sample_type,"celltype" = markers$cluster),paste,collapse = ";")
  GOenrich<-apply(genetable, 1, function(x){
    gene <- unlist(strsplit(x[length(x)],";"))
    GO <- enrichGO(gene = gene,OrgDb = 'org.Hs.eg.db', ont = "BP",keyType = "SYMBOL")
    GO<-subset(GO@result,p.adjust<0.01)
    GO<-cbind(GO,rep(x[1],nrow(GO)),rep(x[2],nrow(GO)))
    return(GO)
  })
  GOenrich<-do.call(rbind,GOenrich)
  colnames(GOenrich)<-c(colnames(GOenrich)[1:12],"celltype","sample_type")
  write.table(GOenrich,paste("./patient_cellfunc/",patientid,"_GO.txt",sep = ""),quote = F,row.names = F,sep = "\t")
  ##Organize cellular data
  ##cell proportion
  prop<-as.data.frame.matrix(Cellratio)
  colnames(prop)<-paste(colnames(prop),"_prop",sep = "")
  prop$celltype<-rownames(prop)
  ##roe
  Roe<-as.data.frame.matrix(Roe)
  colnames(Roe)<-paste(colnames(Roe),"_Ro/e",sep = "")
  Roe$celltype<-rownames(Roe)
  ##Cellular differential disturbance
  if(length(table(SO$sample_type))>1){
    augur_AUC<-as.data.frame(augur$AUC)
    colnames(augur_AUC)<-c("celltype","AUC")
  }else{
    augur_AUC<-cbind(prop$celltype,NA)
    colnames(augur_AUC)<-c("celltype","AUC")
  }
  
  ##Differentially expressed genes in cells
  DEG<-cell_gene
  ##top5 pathway
  top5<-GOenrich %>% group_by(celltype,sample_type) %>% top_n(5, -log(p.adjust))
  aa<-aggregate(top5$ID,list("celltype" = top5$celltype,"sample_type" = top5$sample_type),paste,collapse = ";")
  top5<-spread(aa,celltype,x)
  colnames(top5)<-c("celltype",paste(colnames(top5)[-1],"_GO",sep = ""))
  ##combined
  # Put all data boxes into the list
  df_list <- list(prop,Roe,augur_AUC,DEG,top5)
  # Merge all data boxes using the Reduce function
  cell_table<-Reduce(function(x,y){merge(x,y,all= TRUE,by = 'celltype')},df_list)
  cell_table$patient<-patientid
  write.table(cell_table,paste("./patient_celltable/",patientid,"_celltable.txt",sep = ""),quote = F,row.names = F,sep = "\t")
  
  
  }


