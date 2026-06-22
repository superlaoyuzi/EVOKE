# ---- library packages ----
library(BiocParallel)
library(tradeSeq)
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
library(monocle)
library(slingshot)
library(RColorBrewer)
library(zoo)
setwd("/home/data/t070404/tumor_progression")
##Prepare documents
# ---- Seurat object ----
file<-list.files("./stage1")
file<-file[grep("rds",file)]
# ---- Loop processing data ----
set.seed(10086)
for(i in file[29:length(file)]){
  patientid<-substring(i,1,8)
  SO<-readRDS(paste("./stage1/",i,sep = ""))
  #load slingshot sds
  lineages<-readRDS(paste("./patient_slingshot/",patientid,"_slingshot_data.rds",sep = ""))
  curves <- as.SlingshotDataSet(getCurves(
    data          = lineages,
    thresh        = 1e-1,
    stretch       = 2,
    allow.breaks  = F,
    approx_points = 150
  ))
  pseudotime <- slingPseudotime(curves, na = F)
  pseudotime_na <- slingPseudotime(curves, na = T)
  cellWeights <- slingCurveWeights(curves)
  #load sc counts
  count<-GetAssayData(SO,layer = "counts",assay = "RNA")
  # icMat <- evaluateK(counts = count, 
  #                    sds = crv, nGenes = 200,
  #                    k = 3:10,verbose = T, plot = TRUE)
  ##Select cells and genes
  sel_cells <- split(colnames(SO@assays$RNA@data), SO$celltype)
  sel_cells <- unlist(lapply(sel_cells, function(x) {
    if(length(x)>50){
      return(sample(x, 50))
    }else{
      return(x)
    }
    
  }))

  gv <- as.data.frame(na.omit(scran::modelGeneVar(SO@assays$RNA@data[, sel_cells])))
  gv <- gv[order(gv$bio, decreasing = T), ]
  sel_genes <- sort(rownames(gv)[1:500])

    sceGAM <- fitGAM(
      counts = drop0(SO@assays$RNA@data[sel_genes, sel_cells]),
      pseudotime = pseudotime[sel_cells, ],
      cellWeights = cellWeights[sel_cells, ],
      nknots = 5, verbose = T, parallel = T, sce = TRUE,
      BPPARAM = BiocParallel::MulticoreParam())
  saveRDS(sceGAM,paste("./sceGAM/",patientid,"_GAM.rds",sep = ""))

  
  ##key gene
  startRes <- startVsEndTest(sceGAM,lineages=T)
  associatRes<-associationTest(sceGAM,lineages=T)
  for(n in 1:ncol(pseudotime)){
    startRes_line<-startRes[,grep(paste("lineage",n,"$",sep = ""),colnames(startRes))]
    associatRes_line<-na.omit(associatRes[,grep(paste("_",n,"$",sep = ""),colnames(associatRes))])
    colnames(startRes_line)<-c("S.VS.E_waldStat", "S.VS.E_df","S.VS.E_pvalue","S.VS.E_logFC"  )
    colnames(associatRes_line)<-c("association_waldStat", "association_df","association_pvalue"  )
    
    startRes_line<-subset(startRes_line,S.VS.E_pvalue<0.05)
    associatRes_line<-subset(associatRes_line,association_pvalue<0.05)
    
    startRes_line <- startRes_line[order(startRes_line$S.VS.E_waldStat, decreasing = T), ]
    associatRes_line <- associatRes_line[order(associatRes_line$association_waldStat, decreasing = T), ]
    
    startRes_line<-rownames_to_column(startRes_line,var = "gene")
    associatRes_line<-rownames_to_column(associatRes_line,var = "gene")
    combine_Res<-merge(associatRes_line,startRes_line,by = "gene")
    
    combine_Res$trajectory<-colnames(pseudotime)[n]
    write.csv(combine_Res,paste("./trajectory_gene/",patientid,"_",colnames(pseudotime)[n],"_trajectory_gene.csv",sep = ""),quote = F,row.names = F)
    umap = SO@reductions$umap@cell.embeddings %>%
      as.data.frame() %>% 
      cbind("celltype" = SO@meta.data$celltype, "sample_type" = SO@meta.data$sample_type)
    for(x in combine_Res$gene){
      umap_2<-cbind(cbind(umap,"curve_pseudotime" = pseudotime_na[, n]),"gene_express" = count[x,])
      umap_3<-umap_2[order(umap_2$curve_pseudotime,na.last = NA),]
      line<-as.data.frame(curves@curves[[n]][1])
      colnames(line)<-c("UMAP_1","UMAP_2")
      ggplot(umap_3,aes(x = UMAP_1,y = UMAP_2))+geom_point(aes(color = log(gene_express+1)),size =.5)+scale_color_gradient2(mid = "yellow",high = "red")+
        theme_bw()+ geom_path(data = line,arrow = arrow())+labs(title=x)
      ggsave(paste("./trajectory_gene_express/",patientid,"_",colnames(pseudotime)[n],"_",x,"_express.pdf",sep = ""))
    }
    umap_n<-cbind(cbind(umap,"curve_pseudotime" = pseudotime_na[, n]),t(count[combine_Res$gene,]))
    write.csv(umap_n,paste("./trajectory_gene_umap/",patientid,"_",colnames(pseudotime)[n],"_trajectory_gene_umap.csv",sep = ""),quote = F,row.names = T)
    
    }
  
  

}



