##library packages
library(tidyverse)
library(Seurat)
##load patient data
dataset<-data.table::fread("GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt")
dataset<-column_to_rownames(dataset,colnames(dataset)[1])
head(colnames(dataset))
##Extract the number of each patient
aa<-unlist(lapply(colnames(dataset),function(x){
  a<-strsplit(x,"_")[[1]][1]
  return(strsplit(a,"-")[[1]][1])
}))
table(aa)

##Extract same patient data
patient1<-dataset[,grep("^SMC01",x = colnames(dataset))]
for(i in unique(aa)){
  patientdata<-dataset[,grep(i,x = colnames(dataset))]
  write.table(patientdata,paste(""))
}


##Single cell data preprocessing
##create Seurat object
SO<-CreateSeuratObject(patient1,min.cells = 3, min.features = 200)
##QC
SO[["percent.mt"]] <- PercentageFeatureSet(SO, pattern = "^MT-")
VlnPlot(SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
SO <- subset(SO, subset  = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
##preprocessing
SO <- NormalizeData(object = SO,normalization.method = "LogNormalize", scale.factor = 10000)
SO <- FindVariableFeatures(object = SO,selection.method = "vst", nfeatures = 2000)
SO <- ScaleData(object = SO)
SO <- RunPCA(object = SO, features = VariableFeatures(object = SO))
SO <- FindNeighbors(object = SO,dims = 1:15)
SO <- FindClusters(object = SO,resolution = 0.5)
SO <- RunUMAP(object = SO,dims = 1:15)
DimPlot(object = SO, reduction = "umap",group.by = "orig.ident")
##Calculate differentially expressed genes
markers<-FindAllMarkers(SO)
##Annotate cell types
source("./sctype/sctype_wrapper.R"); 
  ##known_tissue_type = Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
  SO <- run_sctype(SO,known_tissue_type="Intestine",custom_marker_file="./sctype/ScTypeDB_full.xlsx",name="sctype_classification",plot=TRUE)
  ##Further classify immune cells：
  SO <- run_sctype(SO,known_tissue_type="Immune system",custom_marker_file="./sctype/ScTypeDB_full.xlsx",name="sctype_immune_classification",plot=TRUE)
  ##Summarize and generalize cell types
  view(unique(cbind(SO$seurat_clusters,SO$sctype_classification,SO$sctype_immune_classification)))
  SO$celltype<-ifelse(SO$sctype_classification %in% c("Lymphoid cells","Myeloid cells","Unknown"),SO$sctype_immune_classification,SO$sctype_classification)
  ##The epithelial cells in cancer samples are defined as cancer
  SO$celltype[intersect(grep("epithelial",SO$celltype),which(SO$sample_type == "Primary"))]<-"Cancer cells"
  celltype<-unique(cbind(SO$seurat_clusters,SO$sctype_classification,SO$sctype_immune_classification))
  write.table(celltype,paste("./stage1/",patientID[i],"_celltype.txt",sep = ""),quote = F,row.names = F,sep = "\t")
  ##Export data
  saveRDS(SO,paste("./stage1/",patientID[i],"_after_cell_annotation.rds",sep = ""))
  ##UMap Export dimension reduction data for visualization
  umap = SO@reductions$umap@cell.embeddings %>%
    as.data.frame() %>% 
    cbind("celltype" = SO@meta.data$celltype, "sample_type" = SO@meta.data$sample_type)
  write.table(umap,paste("./patient_umap/",patientID[i],"_umap_celltype_sampletype.txt",sep = ""),quote = F,row.names = F,sep = "\t")
