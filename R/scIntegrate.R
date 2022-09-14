#' Integrating scRNA-seq samples into a single Seurat Object
#'
#' 
#' @param Seurat - Seurat object
#' @export


scIntegrate <- function(Seurat_Object) {
  
  type <- class(Seurat_Object)
  
  if (type[1] == "Seurat") {
    
    Seurat_Object <- SCTransform(Seurat_Object)
    Seurat_Object <- FindVariableFeatures(Seurat_Object, selection.method = "vst", nfeatures = 2000)
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
    Seurat_Object <- CellCycleScoring(Seurat_Object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    Seurat_Object <- ScaleData(Seurat_Object, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Seurat_Object))
    Seurat_Object <- RunPCA(Seurat_Object, npcs = 30, verbose = FALSE)
    Seurat_Object <- RunUMAP(Seurat_Object, reduction = "pca", dims = 1:30)
    Seurat_Object <- FindNeighbors(Seurat_Object, reduction = "pca", dims = 1:30)
    Seurat_Object <- FindClusters(Seurat_Object, resolution = 0.5)
    
    Seurat_Object
  
  
  } else if (type[1] == "list") {
  
    Seurat_Object <- lapply(Seurat_Object, SCTransform)
    features <- SelectIntegrationFeatures(object.list = Seurat_Object, nfeatures = 3000)
    
    Seurat_Object <- lapply(X = Seurat_Object, FUN = function(x) {
      x <- ScaleData(x, features = features, verbose = FALSE)
      x <- RunPCA(x, features = features, verbose = FALSE)
    })
    Seurat_Object <- PrepSCTIntegration(object.list = Seurat_Object, anchor.features = features)
    Seurat_Object <- FindIntegrationAnchors(object.list = Seurat_Object, normalization.method = "SCT",
                                        anchor.features = features, reduction = "rpca")
    Seurat_Object <- IntegrateData(anchorset = Seurat_Object, normalization.method = "SCT")
    
    Seurat_Object <- FindVariableFeatures(Seurat_Object, selection.method = "vst", nfeatures = 2000)
    
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
    
    Seurat_Object <- CellCycleScoring(Seurat_Object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    Seurat_Object <- ScaleData(Seurat_Object, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Seurat_Object))
    Seurat_Object <- RunPCA(Seurat_Object, npcs = 30, verbose = FALSE)
    Seurat_Object <- RunUMAP(Seurat_Object, reduction = "pca", dims = 1:30)
    Seurat_Object <- FindNeighbors(Seurat_Object, reduction = "pca", dims = 1:30)
    Seurat_Object <- FindClusters(Seurat_Object, resolution = 0.5)
    
    Seurat_Object
}
  else {
    print('Something went wrong. If you think you have a bug, contact the author
        at mateusz.jakiel98@gmail.com')
  }
}