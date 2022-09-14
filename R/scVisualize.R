#' Creating basic Seurat visualizations as wel as a heatmap if your choice is so
#'
#' 
#' @param Seurat - Seurat object
#' @param features - genes of interest
#' @param create_heatmap - choose Seurat assay
#' @param heatmap_grouping - metadata
#' @keywords single cell, Seurat, correlation
#' @export


scVisualize <- function(Seurat_Object, features, create_heatmap, heatmap_grouping) {
  
  DefaultAssay(Seurat_Object) <- "SCT"
  normalized_counts_rna <- as.data.frame(t(as.matrix(Seurat_Object@assays[['SCT']]@data)))
  
  
  plot_list <- list()
  
  plot_list[[1]] <- DimPlot(Seurat_Object)
  
  plot_list[[2]] <- VlnPlot(Seurat_Object, features = features, pt.size = 0, combine = TRUE, ncol = 2)
  
  plot_list[[3]] <- FeaturePlot(Seurat_Object, features = features, pt.size = 0.5, 
                                label = TRUE)
  plot_list[[4]] <- RidgePlot(Seurat_Object, features = features, ncol = 2)
  
  plot_list[[5]] <- DoHeatmap(subset(Seurat_Object, downsample = 100), features = features, size = 3)
  
  if (create_heatmap == TRUE) {
    normalized_counts_rna <- normalized_counts_rna[colnames(normalized_counts_rna) %in% features]
    heatmap_grouping <- as.factor(heatmap_grouping)
    normalized_counts_rna <- cbind(normalized_counts_rna, heatmap_grouping)
    
    heatmapt_input <- melt(as.data.table(normalized_counts_rna), id.vars = c('heatmap_grouping'))
    heatmapt_input <- heatmapt_input %>%
      group_by(heatmap_grouping, variable) %>%
      summarize(Mean = mean(value)) %>%
      as.data.frame()
    heatmapt_input$Mean <- as.numeric(heatmapt_input$Mean)
    
    heatmapt_input$variable <- factor(heatmapt_input$variable, levels = rev(features[features %in% heatmapt_input$variable]))
    
    heatmap <- ggplot(heatmapt_input,aes(x=heatmap_grouping,y=variable,fill=Mean)) + geom_tile(colour="white",size=0.05) + 
      labs(x=paste(heatmap_grouping),y="Genes") + labs(title = "Expression heatmap") 
    plot_list[[6]] <- heatmap + scale_fill_viridis(option = "A", discrete = FALSE) 
  }
  
  plot_list
  
}

