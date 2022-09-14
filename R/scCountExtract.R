#' Extracting counts from the Seurat object to use other types of analysis 
#'
#' 
#' @param Seurat - Seurat object
#' @param use_metadata - TRUE/FALSE whether you want to extract specific metadata
#' @param features - metadata you want to use
#' @keywords single cell, Seurat, correlation
#' @export

scCountExtract <- function(Seurat, use_metadata, features) {
  
  if (use_metadata == TRUE) {
    
    normalized_counts_rna <- as.data.frame(t(as.matrix(Seurat@assays[['SCT']]@data)))
    meta <- data.frame(matrix(0, ncol = length(features), nrow = nrow(normalized_counts_rna)))

    for (i in 1:length(features)) {
      
      meta[i] <- Seurat@meta.data[colnames(Seurat@meta.data) == features[i]]

    }
    colnames(meta) <- features
    
    
    normalized_counts_rna <- data.frame(normalized_counts_rna, meta)
    normalized_counts_rna
    
  } else if (use_metadata == FALSE) {
    
    normalized_counts_rna <- as.data.frame(t(as.matrix(Seurat@assays[['SCT']]@data)))
    normalized_counts_rna
    
  } else {
    print("Error, wrong input data. If you think you have a bug, contact the author
        at mateusz.jakiel98@gmail.com")
  }
}