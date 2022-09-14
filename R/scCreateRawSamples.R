#' Creating a list of scRNA-seq objects for further analysis
#'
#' This function allows you to combine raw scRNA-seq files into a list of Seurat objects.
#' 
#' @param x - The directory of your files
#' @keywords Seurat objecy
#' @export

scCreateRawSamples <- function(x = getwd(), min.cells, min.features) {

filenames <- dir(path=x)

if (TRUE %in% grepl('.tsv', filenames)) {
  if (length(filenames)%%3 == 0 & length(filenames)>3) {
    num <- grepl("barcodes", filenames)
    val <- which(num)
    
    matrix_mash <- function(x) {
      barcode_path <- filenames[x]
      features_path <- filenames[x+1]
      matrix_path <- filenames[x+2]
      barcode_names = read.delim(barcode_path, 
                                 header = FALSE,
                                 stringsAsFactors = FALSE)
      feature_names = read.delim(features_path, 
                                 header = FALSE,
                                 stringsAsFactors = FALSE)
      mat <- readMM(file = matrix_path)
      colnames(mat) = barcode_names$V1
      rownames(mat) = feature_names$V2
      mat
    }
    
    mash_it <- lapply(val, matrix_mash)
    
    data <- lapply(seq_along(mash_it), function(i) {
      seurat_object <- CreateSeuratObject(counts = mash_it[[i]], 
                                          min.cells = 5, min.features = 500)
    }
  )
    
    mito <- function(i) {PercentageFeatureSet(data[[i]], "^MT-", col.name = "percent_mito")}
    ribo <- function(i) {PercentageFeatureSet(data[[i]], "^RP[SL]", col.name = "percent_ribo")}
    hb <- function(i) {PercentageFeatureSet(data[[i]], '^HB[^(P)]', col.name = "percent_hb")}
    
    data <- lapply(seq_along(data), mito)
    data <- lapply(seq_along(data), ribo)
    data <- lapply(seq_along(data), hb)
    
  data
}
  

  else if (length(filenames)%%3 == 0 & length(filenames)==3) {

    num <- grepl("barcodes", filenames)
    val <- which(num)
    
    barcode_path <- filenames[val]
    features_path <- filenames[val+1]
    matrix_path <- filenames[val+2]
    barcode_names = read.delim(barcode_path, 
                               header = FALSE,
                               stringsAsFactors = FALSE)
    feature_names = read.delim(features_path, 
                               header = FALSE,
                               stringsAsFactors = FALSE)
    mat <- readMM(file = matrix_path)
    colnames(mat) = barcode_names$V1
    rownames(mat) = feature_names$V2
    data <- CreateSeuratObject(counts = mat, 
                               min.cells = 5, min.features = 500)
    
    data <- PercentageFeatureSet(data, "^MT-", col.name = "percent_mito")
    data <- PercentageFeatureSet(data, "^RP[SL]", col.name = "percent_ribo")
    data <- PercentageFeatureSet(data, '^HB[^(P)]', col.name = "percent_hb")
    
    data

  }
  
} else if (TRUE %in% grepl('.txt', filenames)) {
  
    if (length(filenames) > 1) {
      
      val <- 1:length(filenames)
      
    table_mash <- function(x) {
      table <- read.table(filenames[x], row.names = 1, header = TRUE)
      table
    }
    
    data <- lapply(val, table_mash)
    
    data <- lapply(seq_along(mash_it), function(i) {
      seurat_object <- CreateSeuratObject(counts = mash_it[[i]], 
                                          min.cells = 5, min.features = 500)
    }
    )
    
    mito <- function(i) {PercentageFeatureSet(data[[i]], "^MT-", col.name = "percent_mito")}
    ribo <- function(i) {PercentageFeatureSet(data[[i]], "^RP[SL]", col.name = "percent_ribo")}
    hb <- function(i) {PercentageFeatureSet(data[[i]], '^HB[^(P)]', col.name = "percent_hb")}
    
    data <- lapply(seq_along(data), mito)
    data <- lapply(seq_along(data), ribo)
    data <- lapply(seq_along(data), hb)
    
    data

    } else {
      data <- read.table(filenames, row.names = 1, header = TRUE)
      data <- CreateSeuratObject(counts = data, 
                                          min.cells = 5, min.features = 500)
      data <- PercentageFeatureSet(data, "^MT-", col.name = "percent_mito")
      data <- PercentageFeatureSet(data, "^RP[SL]", col.name = "percent_ribo")
      data <- PercentageFeatureSet(data, '^HB[^(P)]', col.name = "percent_hb")
      
      data
      
      
    }
}
  else {"Data is not in barcode/features/matrix format. If you wish to use a different format that is not supported, please
    contact the author at mateusz.jakiel98@gmail.com"}
}