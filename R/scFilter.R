#' Use the results of IntroQC for quality control application
#'
#' 
#' @param Seurat - Seurat object
#' @param QC_results - results of IntroQC
#' @keywords single cell, quality control
#' @export

scFilter <- function(QC_results, Seurat_Object) {
  
  type <- class(data)
  x <- QC_results[[1]]
  
  if (type[1] == "Seurat") {
    
    DefaultAssay(data) <- "RNA"
    Seurat_Objet <- subset(x = Seurat_Object, subset = nFeature_RNA > 300 & nFeature_RNA < x$temp[1] & 
                         percent_mito < x$temp[2] & percent_ribo > x$temp[3] & percent_hb < x$temp[4])
  
  } else if (type[1] == "list") {
  
    Seurat_Object <- lapply(seq_along(Seurat_Object), function(i) {
      Seurat_Object <- subset(x = Seurat_Object[[i]], subset = nFeature_RNA > 300 & nFeature_RNA < x$temp[1] &
                             percent_mito < x$temp[2] & percent_ribo > x$temp[3] & percent_hb < x$temp[4])
    
  }
)
  } else {
    print("Wrong data type. If you think you have a bug, contact the author
        at mateusz.jakiel98@gmail.com")
  }
}