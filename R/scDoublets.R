#' Finding doublets in Seurat Objects
#'
#' 
#' @param Seurat_Object - Seurat object
#' @export

scDoublets <- function(Seurat_Object) {
  
  sweep.res.Seurat_Object <- paramSweep_v3(Seurat_Object, PCs = 1:10, sct = TRUE, num.cores = detectCores() - 1)
  sweep.stats_Seurat_Object <- summarizeSweep(sweep.res.Seurat_Object, GT = FALSE)
  bcmvn_Seurat_Object <- find.pK(sweep.stats_Seurat_Object)

  homotypic.prop <- modelHomotypic(Seurat_Object@meta.data$seurat_clusters)
  nExp_poi <- round(0.075*nrow(Seurat_Object@meta.data))  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  Seurat_Object <- doubletFinder_v3(Seurat_Object, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

  metadata <- Seurat_Object@meta.data
  colnames(metadata)[15] <- "doublet_finder"
  Seurat_Object@meta.data <- metadata 

  Seurat_Object <- subset(Seurat_Object, doublet_finder == "Singlet")
  Seurat_Object

}