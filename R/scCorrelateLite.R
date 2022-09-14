#' Creating a dataframe of correlated genes to our genes of interest
#'
#' 
#' @param Seurat - Seurat object
#' @param genes - genes of interest
#' @param assay - choose Seurat assay
#' @param corr_test - which type of correlation test 
#' @param p_val - p-value cutoff for FDR
#' @keywords single cell, Seurat, correlation
#' @export

scCorrelateLite <- function(Seurat, gene_symbols) {
  
normalized_counts_rna <- as.data.frame(t(as.matrix(Seurat@assays[['SCT']]@data)))

smp_size <- floor(0.05 * nrow(normalized_counts_rna))
train_ind <- sample(seq_len(nrow(normalized_counts_rna)), size = smp_size)
normalized_counts_rna <- normalized_counts_rna[train_ind, ]

normalized_counts_rna <- normalized_counts_rna[vapply(normalized_counts_rna, function(x) length(unique(x)) > 1, logical(1L))] 

gene_symbols <- toupper(gene_symbols)

colnames(normalized_counts_rna) <- toupper(colnames(normalized_counts_rna))

gene_symbols <- gene_symbols[gene_symbols %in% colnames(normalized_counts_rna)]

corelation_results <- list()
res <- list()
for(k in 1:length(gene_symbols)){
  
  tempData <- normalized_counts_rna[,!colnames(normalized_counts_rna) %in% gene_symbols[k]]
  
  for(i in 1:length(colnames(tempData))){
    corr_results <- cor.test(normalized_counts_rna[, colnames(normalized_counts_rna)==gene_symbols[k]], tempData[, i], method = "pearson")
    names(corr_results$statistic) <- NULL
    corelation_results[[i]] <- data.frame(first_var=gene_symbols[k],
                                          second_var=colnames(tempData)[i],
                                          cor_stat=corr_results$estimate[[1]],
                                          pvalue=corr_results$p.value)
    
  }
 res[[k]] <- corelation_results
 next
}

results_corelation <- list()

for (i in 1:length(gene_symbols)) {
  
  results_corelation[[i]] <- rbindlist(res[[i]])

  results_corelation[[i]]$adj_pval <- p.adjust(results_corelation[[i]]$pvalue, method = "fdr", length(results_corelation[[i]]$pvalue))

  results_corelation[[i]] <- na.omit(results_corelation[[i]])
  
  results_corelation[[i]] <- results_corelation[[i]][!results_corelation[[i]]$adj_pval > 0.05,]
}







results_corelation
}