# ---- library packages ----
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
library(monocle)
library(slingshot)
library(RColorBrewer)
library(zoo)
setwd("/home/data/t070404/tumor_progression")
# Define some color palette
pal <- c(scales::hue_pal()(8), RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))
set.seed(10086)
pal <- rep(sample(pal, length(pal)), 200)


# ---- Seurat object ----
file<-list.files("./stage1")
file<-file[grep("rds",file)]
# ---- Pseudotime analysis ----
for(i in file[29:59]){
  patientid<-substring(i,1,8)
  SO <- readRDS(paste("./stage1/",i,sep = ""))
  umap = SO@reductions$umap@cell.embeddings %>%
    as.data.frame() %>% 
    cbind("celltype" = SO@meta.data$celltype, "sample_type" = SO@meta.data$sample_type)
  mm <- sparse.model.matrix(~ 0 + factor(SO$celltype))
  colnames(mm) <- levels(factor(SO$celltype))
  centroids2d <- as.matrix(t(t(SO@reductions$umap@cell.embeddings) %*% mm) / Matrix::colSums(mm))
  #slingshot
  # sce <- as.SingleCellExperiment(SO, assay = "RNA")
  # sim <- slingshot(sce,      #输入单细胞对象
  #                             reducedDim = 'UMAP',  #降维方式
  #                             clusterLabels = sce$celltype,  #cell类型
  #                             approx_points = 150,
  #                             reweight = FALSE)
  lineages <- as.SlingshotDataSet(getLineages(
    data           = SO@reductions$umap@cell.embeddings,
    clusterLabels  = SO$celltype,
    dist.method    = "slingshot"# It can be: "simple", "scaled.full", "scaled.diag", "slingshot" or "mnn"
  ))
  
  plot(SO@reductions$umap@cell.embeddings, col =  pal[SO$seurat_clusters],cex = .5, pch = 16)
  lines(lineages, lwd = 1, col = "black", cex = 2)
  text(centroids2d, labels = rownames(centroids2d), cex = 0.8, font = 2, col = "lightgrey")
  #View Path
  # Define curves
  curves <- as.SlingshotDataSet(getCurves(
    data          = lineages,
    thresh        = 1e-1,
    stretch       = 2,
    allow.breaks  = F,
    approx_points = 150
  ))
  # Plots
  {
    plot(SO@reductions$umap@cell.embeddings, col =  pal[SO$seurat_clusters],pch = 16)
    lines(curves, lwd = 2, col = "black")
    text(centroids2d, labels = rownames(centroids2d), cex = 1, font = 2)
  }
  ##compute the differentiation pseudotime
  pseudotime <- slingPseudotime(curves, na = F)
  cellWeights <- slingCurveWeights(curves)
  
  x <- rowMeans(pseudotime)
  x <- x / max(x)
  o <- order(x)
  
  {
    plot(SO@reductions$umap@cell.embeddings[o, ],
         main = paste0("pseudotime"), pch = 16, cex = 0.4, axes = T, xlab = "", ylab = "",
         col = colorRampPalette(c("grey70", "orange3", "firebrick", "purple4"))(99)[x[o] * 98 + 1]
    )
    points(centroids2d, cex = 2.5, pch = 16, col = "#FFFFFF99")
    text(centroids2d, labels = rownames(centroids2d), cex = 1, font = 2)
  }

  ##Extract the number of paths
  pseudotime_na <- slingPseudotime(curves, na = T)
  Lineage<-colnames(pseudotime_na)
  ##visualization
  for(j in Lineage){
    # Sort by pseudo time
    sorted_cells <- colnames(SO)[order(pseudotime_na[, j], na.last = NA)]
    sample<-factor(unique(SO$sample_type),levels = c("Health","Inflamed","Margin","Primary","Metastasis"))
    if(length(sample)<2){next}
    if(length(sample) == 2){
      df <- data.frame(
        time = 1:length(SO@meta.data[sorted_cells,"sample_type"]),
        char = SO@meta.data[sorted_cells,"sample_type"]
      ) %>%
        mutate(
          cum_1 = cumsum(char == sample[1]),  # sample1 
          cum_2 = cumsum(char == sample[2])   # sample2 
        )
    }else if(length(sample) == 3){
      df <- data.frame(
        time = 1:length(SO@meta.data[sorted_cells,"sample_type"]),
        char = SO@meta.data[sorted_cells,"sample_type"]
      ) %>%
        mutate(
          cum_1 = cumsum(char == sample[1]),  # sample1 
          cum_2 = cumsum(char == sample[2]),   # sample2 
          cum_3 = cumsum(char == sample[3])
        )
    }
    
    # Sliding window calculation frequency (window size=100)
    window_size <- ceiling(ncol(SO)/5000)*100
    df$prop_1 <- rollapply(
      data = (df$char == sample[1]),
      width = window_size,
      FUN = function(x) mean(x),
      fill = NA,
      align = "right"
    )
    
    df$prop_2 <- rollapply(
      data = (df$char == sample[2]),
      width = window_size,
      FUN = function(x) mean(x),
      fill = NA,
      align = "right"
    ) 
    if(length(sample) == 3){
      df$prop_3 <- rollapply(
        data = (df$char == sample[3]),
        width = window_size,
        FUN = function(x) mean(x),
        fill = NA,
        align = "right"
      ) 
    }
    if(length(sample) == 2){    
      df_long_prop <- df %>%
      pivot_longer(cols = c(prop_1, prop_2), names_to = "Character", values_to = "Proportion")
    }else if(length(sample) == 3){
      df_long_prop <- df %>%
        pivot_longer(cols = c(prop_1, prop_2, prop_3), names_to = "Character", values_to = "Proportion")
      
    }
    
    if(length(sample) == 2){    
      p<-ggplot(df_long_prop, aes(x = time, y = Proportion, color = Character)) +
        geom_line(linewidth = 1) +
        scale_color_manual(values = c("prop_1" = "blue4", "prop_2" = "red4"), labels = sample[1:2]) +
        labs(title = paste("sliding window (size", window_size, ") frequency"), x = "pseudotime", y = "frequency") +
        ylim(0, 1) +
        theme_minimal()
    }else if(length(sample) == 3){
      p<-ggplot(df_long_prop, aes(x = time, y = Proportion, color = Character)) +
        geom_line(linewidth = 1) +
        scale_color_manual(values = c("prop_1" = "blue4", "prop_2" = "yellow3","prop_3" = "red4"), labels = sample[1:3]) +
        labs(title = paste("sliding window (size", window_size, ") frequency"), x = "pseudotime", y = "frequency") +
        ylim(0, 1) +
        theme_minimal()
      
    }


    
    ##Export data
    write.csv(df_long_prop,paste("./lineage_slingshot_table/",patientid,"_",j,"frenquce.csv",sep = ""),quote = F,row.names = F)
    
    ggsave(paste("./lineage_slingshot_plot/",patientid,"_",j,"frenquce.pdf",sep = ""),p)
    ##Single path pseudo time sequence
    pseudotime_na <- slingPseudotime(curves, na = T)
    umap_2<-cbind(umap,"curve_pseudotime" = pseudotime_na[, j])
    pdf(paste("./curve_slingshot_plot/",patientid,"_",j,"_curve.pdf",sep = ""))  
    plot(umap_2[order(umap_2$curve_pseudotime,na.last = NA),1:2],
         main = gsub("Lineage","Curve",j), pch = 16, cex = 0.4, axes = T, xlab = "", ylab = "",
         col = colorRampPalette(c("grey70", "orange3", "firebrick", "purple4"))(99)[x[o] * 98 + 1]
    )
    lines(curves@curves[[j]], lwd = 2, col = "black")
    dev.off()
    umap_out<-umap_2[order(umap_2$curve_pseudotime,na.last = T),]
    write.csv(rownames_to_column(umap_out,var = "cellnames"),paste("./curve_umap_table/",patientid,"_",j,"curve_umap_table.csv",sep = ""),quote = F,row.names = F)
    write.csv(curves@curves[[j]]$s,paste("./curve_umap_table/",patientid,"_",j,"curve_line_table.csv",sep = ""),quote = F,row.names = F)
    }
  ##Patient's proposed time sequence
  
  write.csv(cbind(umap,"pseudotime" = x)[o,],paste("./patient_pseudotime_umap/",patientid,"_pseudotime_umap_celltype_sampletype.csv",sep = ""),quote = F,row.names = T)
  saveRDS(lineages,paste("./patient_slingshot/",patientid,"_slingshot_data.rds",sep = ""))
  aa<-t(as.data.frame(lapply(curves@lineages,paste,collapse = "->")))
  las<-lapply(grep("Curve",capture.output(curves),value = T),function(x){
    grep("^\\d",unlist(lapply(strsplit(x,"[\t]"),strsplit,"[: ]")),value = T)
  })
  TRAN<-as.data.frame(cbind(aa,do.call(rbind,las)))
  colnames(TRAN)<-c("cell_trajectory","Curve_length","Curve_samples")
  TRAN$patient<-patientid
  TRAN$StartStage<-sample[1]
  TRAN$EndStage<-sample[length(sample)]
  write.table(TRAN,paste("./patient_trajectory_table/",patientid,"_trajectory_table.txt",sep = ""),quote = F,row.names = T,sep = "\t")
  }



