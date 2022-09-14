#' Feature selection using the Monte Carlo Feature Selection algorithm
#'
#' 
#' @keywords single cell, Seurat, correlation
#' @export


scFeatureSelection <- function(Seurat, use_metadata, features) {
  
  normalized_counts_rna <- as.data.frame(t(as.matrix(Seurat@assays[['SCT']]@data)))
  colnames(normalized_counts_rna) <- toupper(normalized_counts_rna)
  genes <- toupper(features)
  genes <- genes[-!genes %in% colnames(normalized_counts_rna)]
  
  if(use_metadata == TRUE) {
  
    meta <- Seurat@meta.data[colnames(Seurat@meta.data) == features]
    
    colnames(meta) <- "feature"
    
    normalized_counts_rna <- data.frame(normalized_counts_rna, meta)
    
    MCFS <- mcfs(formula = feature~., normalized_counts_rna, 
                      threadsNumber = 12, mode = 2, 
                      finalRuleset = FALSE, splitSetSize	= 1000, finalCV = TRUE)
    
    topFeatures <- MCFS$RI$attribute[1:MCFS$cutoff_value]
    plotRI <- plot(MCFS, type = "ri", size = 100, gg = TRUE)
    
    plot(MCFS, cv_measure = "acc", size = 20)
    plot(MCFS, type = "cv", cv_measure = "wacc")
    
    gid <- build.idgraph(MCFS, size = 25)
    plot(gid)
  }
  else if (use_metadata == FALSE) {
    if (length(features) > 1) {
  
                        mcfs_polyfeat <- list()
                        
                        for (i in 1:length(features)) {
                        
                          mcfs_polyfeat[[i]][1] <- mcfs(formula = features[i]~., normalized_counts_rna, 
                                                        threadsNumber = 12, mode = 2, 
                                                        finalRuleset = FALSE, splitSetSize	= 1000, finalCV = TRUE)
                        
                          mcfs_polyfeat[[i]][2] <- mcfs_polyfeat[[i]]$RI$attribute[1:mcfs_polyfeat$cutoff_value]
                          mcfs_polyfeat[[i]][3] <- plot(mcfs_polyfeat[[i]], cv_measure = "acc", size = 20)

                          mcfs_polyfeat[[i]][4] <- plot(mcfs_polyfeat[[i]], type = "cv", cv_measure = "wacc")
                        
                          MCFS[[i]][5] <- MCFS[[i]]$ID[MCFS[[i]]$ID$weight > MCFS[[i]][["cutoff"]][["minID"]][3],]
                        
                        }
  names(MCFS) <- c("MCFS experiment", "SelectedFeatures", "Top25Features", "ClassifierComparison", "Interdepencies")
  mcfs_polyfeat
  
    } 
    else {

        MCFS <- list()
        
        MCFS[[1]] <- mcfs(formula = features~., normalized_counts_rna, 
                                      threadsNumber = 12, mode = 2, 
                                      finalRuleset = FALSE, splitSetSize	= 1000, finalCV = TRUE)
        
        MCFS[[2]] <- MCFS[[1]]$RI$attribute[1:MCFS[[1]]$cutoff_value]
        
        MCFS[[3]] <- plot(MCFS[[1]], cv_measure = "acc", size = 20)
        
        MCFS[[4]] <- plot(MCFS[[1]], type = "cv", cv_measure = "wacc")
        
        MCFS[[5]] <- MCFS[[1]]$ID[MCFS[[1]]$ID$weight > MCFS[[1]][["cutoff"]][["minID"]][3],]
        
    }
    names(MCFS) <- c("MCFS experiment", "SelectedFeatures", "Top25Features", "ClassifierComparison", "Interdepencies")
    MCFS
} 
  else {
print("Error, define data correctly. If you think you have a bug, contact the author
        at mateusz.jakiel98@gmail.com")
}
}