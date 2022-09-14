#' scRNA-seq quality control checks
#'
#'This function calculates a couple of basic quality control metrics and returns a plot containing
#'multiple violin plots with the mean and median of:
#'- nFeature_RNA
#'- percentage of mitochondrial genes in cell
#'- percentage of ribosomal genes in cell
#'- percentage of red blood cell genes
#' 
#' The function returns the Seurat file (or list of files) with added quality control metrics and a 
#' plot containing the quality control metrics as violin plots.
#' To select the best object in a pre-integrated list of cells, the file with the highest number of cells will be
#' used as a benchmark. 
#' 
#' @param x - the Seurat file (or list of files)
#' @keywords Seurat, quality control
#' @export

scIntroQC <- function(x) {
  
data <- x
type <- class(data)

if (type[1] == "Seurat"){

DefaultAssay(data) <- "RNA"

feats <- c("nFeature_RNA", "percent_mito", "percent_ribo", "percent_hb")

plot <- VlnPlot(data, features = feats, pt.size = 0, ncol = 2) + NoLegend()

temp <- numeric()
for (i in 1:length(feats)) {
  x <- mean(eval(parse(text = paste0("data@meta.data$", feats[i]))))
  temp <- c(temp, x)
}
h <- data.frame(temp, feats)

for (i in 1:length(feats)) {
  plot[[i]] <- plot[[i]] + geom_hline(yintercept = h[i, 1], linetype="dashed", color = "black") +
    geom_text(x = 1, y = h[i, 1], label = paste("mean =", round(h[i, 1], digits = 3)), size = 5, hjust = 0.5, vjust = -1)
  }

y <- list()
y[[1]] <- h
y[[2]] <- plot
return(y)

} 

else if (type[1] == "list"){
  
choice <- data.frame(choice = c(1:length(data)))

for (i in 1:length(data)){
  choice[i, ] <- ncol(data[[i]])
}

num <- which.max(choice[,1])

chosen_data <- data[[num]]
Idents(chosen_data) <- as.factor(paste("Sample", num))

feats <- c("nFeature_RNA", "percent_mito", "percent_ribo", "percent_hb")

plot <- VlnPlot(chosen_data, features = feats, pt.size = 0, ncol = 2) + NoLegend()

temp <- numeric()
for (i in 1:length(feats)) {
  x <- mean(eval(parse(text = paste0("chosen_data@meta.data$", feats[i]))))
  temp <- c(temp, x)
}
h <- data.frame(temp, feats)

for (i in 1:length(feats)) {
    plot[[i]] <- plot[[i]] + geom_hline(yintercept = h[i, 1], linetype="dashed", color = "black") +
    geom_text(x = 1, y = h[i, 1], label = paste("mean =", round(h[i, 1], digits = 3)), size = 5, hjust = 0.5, vjust = -1)
  }
y <- list()
y[[1]] <- h
y[[2]] <- plot
return(y)

} 
else {
  
  print("You have input invalid data. Please make sure you are either inputting a single sample of 3 files (barcodes, feature names and counts)
  or a list of samples. If you wish to create a list please use the CreateSampleList() function. If you think you have a bug, contact the author
        at mateusz.jakiel98@gmail.com")
  
  }
}
        






