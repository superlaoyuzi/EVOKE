# ---- library packages ----
library(Seurat)
# ---- Seurat  ----
file<-list.files("./stage1")
file<-file[grep("rds",file)]
# ---- Genetic differences ----
for(i in file){
  patientid<-substring(i,1,8)
  SO<-readRDS(paste("./stage1/",i,sep = ""))
  if(length(unique(SO$sample_type))<2){next}
  ##For patients with two stages, calculate the differences of the same cell type in different stages
  celltype<-unique(SO$celltype)
  allmarkers<-c()
  for(j in celltype){
    if(j == "Unknown"){next}
    SO_cell<-subset(SO,celltype == j)
    if(length(unique(SO_cell$sample_type))<2){next}
    sampletype<-factor(unique(SO_cell$sample_type),levels = c("Health","Inflamed","Margin","Primary","Metastasis"))
    if(length(sampletype) == 2){
      cell1<-colnames(subset(SO_cell,sample_type == sampletype[1]))
      cell2<-colnames(subset(SO_cell,sample_type == sampletype[2]))
      markers<-rownames_to_column(FindMarkers(SO_cell,ident.1  = cell2,ident.2  = cell1),var = "gene")
      markers$from_stage<-sampletype[1]
      markers$to_stage<-sampletype[2]
    }else if(length(sampletype) == 3){
      cell1<-colnames(subset(SO_cell,sample_type == sampletype[1]))
      cell2<-colnames(subset(SO_cell,sample_type == sampletype[2]))
      cell3<-colnames(subset(SO_cell,sample_type == sampletype[3]))
      marker1<-rownames_to_column(FindMarkers(SO_cell,ident.1  = cell2,ident.2  = cell1),var = "gene")
      marker1$from_stage<-sampletype[1]
      marker1$to_stage<-sampletype[2]
      marker2<-rownames_to_column(FindMarkers(SO_cell,ident.3  = cell3,ident.2  = cell2),var = "gene")
      marker2$from_stage<-sampletype[2]
      marker2$to_stage<-sampletype[3]
      markers<-rbind(marker1,marker2)
    }else if(length(sampletype) == 4){
      cell1<-colnames(subset(SO_cell,sample_type == sampletype[1]))
      cell2<-colnames(subset(SO_cell,sample_type == sampletype[2]))
      cell3<-colnames(subset(SO_cell,sample_type == sampletype[3]))
      cell4<-colnames(subset(SO_cell,sample_type == sampletype[4]))
      marker1<-rownames_to_column(FindMarkers(SO_cell,ident.1  = cell2,ident.2  = cell1),var = "gene")
      marker1$from_stage<-sampletype[1]
      marker1$to_stage<-sampletype[2]
      marker2<-rownames_to_column(FindMarkers(SO_cell,ident.3  = cell3,ident.2  = cell2),var = "gene")
      marker2$from_stage<-sampletype[2]
      marker2$to_stage<-sampletype[3]
      marker3<-rownames_to_column(FindMarkers(SO_cell,ident.3  = cell4,ident.2  = cell3),var = "gene")
      marker3$from_stage<-sampletype[3]
      marker3$to_stage<-sampletype[4]
      marker3<-rbind(marker1,marker2,marker3)
    }
    markers$celltype<-j
    markers<-subset(markers,p_val_adj<0.05)
    allmarkers<-rbind(allmarkers,markers)
  }
  allmarkers$patient<-patientid
  write.table(allmarkers,paste("./patient_key_gene/",patientid,"_key_gene.txt",sep = ""),quote = F,row.names = F,sep = "\t")
}

