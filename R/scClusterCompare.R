#' Comparing FindMarkers results
#'
#' This function allows you to compare the results of the FindMarkers function from the Seurat package.
#' This function is for people who have one set of markers that they feel are validated for a specific cell type
#' and want to reuse them for annotation in a new dataset. If you do not know, which marker genes were crucial for the annotation
#' this can be difficult without the help of a biologist. That is why this package allows to automize this process, by calculating
#' a score that can to some degree show which clusters are similar to which.
#' 
#' The function returns a dataframe with n number of columns, corresponding to the number of original dataset clusters plus one with the score.
#' The number of rows corresponds to the number of clusters in the new dataset.
#' @param x - the validated markers
#' @param y - your newly generated markers from a new dataset
#' @keywords Seurat, cell type, markers
#' @export

scClusterCompare <- function(known_markers, novel_markers) {
old_markers <- known_markers
old_markers <- old_markers[7:8]
old_markers$cluster <- as.factor(old_markers$cluster)

new_markers <- novel_markers
new_markers <- new_markers[7:8]
new_markers$cluster <- as.numeric(new_markers$cluster)
new_markers$cluster <- new_markers$cluster + 1
new_markers$cluster <- as.factor(new_markers$cluster)

factors <- names(summary(old_markers$cluster))
nr <- data.frame()


for (y in 1:length(summary(new_markers$cluster))) {
  for (x in 1:length(factors)) {
    nr[y, x] <- c((sum(new_markers$gene[new_markers$cluster == y] %in%  
                         old_markers$gene[old_markers$cluster == factors[x]])/
                     length(new_markers$gene[new_markers$cluster == y])) * 100)
    
  }
}

colnames(nr) <- c(factors)
rownames(nr) <- 1:length(summary(new_markers$cluster))

nr <- as.data.frame(t(nr))

plot_data <- list()

for (i in 1:ncol(nr)) {
  plot_data[[i]] <- as.data.frame(c(nr[order(nr[i], decreasing = TRUE), ][i], as.data.frame(rownames(nr[order(nr[i], decreasing = TRUE), ][i]))))
  plot_data[[i]] <- plot_data[[i]][1:3,]
  colnames(plot_data[[i]]) <- c("New cluster", "Top comparison clusters")
}




plot <- list()

for (i in 1:length(plot_data)) {
  plot[[i]] <- ggplot(plot_data[[i]], aes(x = `Top comparison clusters`, y = `New cluster`)) + 
    geom_col(aes(fill = `Top comparison clusters`)) +
    geom_text(mapping = aes(label = round(`New cluster`)), position = position_stack(vjust = .5)) + 
    labs(y = paste("New cluster", i-1)) +
    theme(legend.position="none")
}

do.call('grid.arrange', c(plot))

nr
}