#' Feature selection using gradient boosting 
#'
#' 
#' @keywords single cell, Seurat, correlation
#' @export


scFeatureSelectionLite <- function(Seurat, use_metadata, features) {
  
  normalized_counts_rna <- as.data.frame(t(as.matrix(Seurat@assays[['SCT']]@data)))
  colnames(normalized_counts_rna) <- toupper(colnames(normalized_counts_rna))

  
  if(use_metadata == TRUE) {
    
    features <- as.factor(features)
    meta <- Seurat@meta.data[colnames(Seurat@meta.data) == features]
    
    smp_size <- floor(0.15 * nrow(normalized_counts_rna))
    train_ind <- sample(seq_len(nrow(normalized_counts_rna)), size = smp_size)
    
    train <- normalized_counts_rna[train_ind, ]
    meta <- meta[train_ind, ]
    
    smp_size <- floor(0.8 * nrow(train))
    train_ind <- sample(seq_len(nrow(train)), size = smp_size)
    
    meta <- meta[train_ind]
    
    train <- na.omit(train[train_ind, ])
    test <- na.omit(train[-train_ind, ])
    meta_train <- meta
    meta_test <- meta[-train_ind]
    
    train_x = data.matrix(train)
    test_x = data.matrix(test)

    xgb_train = xgb.DMatrix(data = train_x, label = meta_train)
    xgb_test = xgb.DMatrix(data = test_x, label = meta_test)
    
    watchlist <- list(train=xgb_train, test=xgb_test)
    
    model = xgb.train(data = xgb_train, max.depth = 3, watchlist = watchlist, nrounds = 100)
    
    importance <- xgb.importance(colnames(xgb_train), model = model)
    
    importance
    
  }
  else if (use_metadata == FALSE) {
    
    smp_size <- floor(0.15 * nrow(normalized_counts_rna))
    train_ind <- sample(seq_len(nrow(normalized_counts_rna)), size = smp_size)
    
    train <- normalized_counts_rna[train_ind, ]

    smp_size <- floor(0.8 * nrow(train))
    train_ind <- sample(seq_len(nrow(train)), size = smp_size)
    
    train <- na.omit(train[train_ind, ])
    test <- na.omit(train[-train_ind, ])
    
    if (length(features) > 1) {
      
      importance <- list()
      for (i in 1:length(genes)) {
        
        train_x = data.matrix(train[,!colnames(train) %in% genes[i]])
        train_y = train[,colnames(train) == genes[i]]
        test_x = data.matrix(test[,!colnames(test) %in% genes[i]])
        test_y = test[,colnames(test) == genes[i]]
        
        xgb_train = xgb.DMatrix(data = train_x, label = train_y)
        xgb_test = xgb.DMatrix(data = test_x, label = test_y)
        
        watchlist <- list(train=xgb_train, test=xgb_test)
        
        model = xgb.train(data = xgb_train, max.depth = 3, watchlist = watchlist, nrounds = 100)

        importance[i][[1]] <- xgb.importance(colnames(xgb_train), model = model)
      }
      names(importance) <- genes
      
      importance
  }    
    else {
      
      train_x = data.matrix(train[,!colnames(train) %in% genes])
      train_y = train[,colnames(train) == genes]
      test_x = data.matrix(test[,!colnames(test) %in% genes])
      test_y = test[,colnames(test) == genes]
      
      xgb_train = xgb.DMatrix(data = train_x, label = train_y)
      xgb_test = xgb.DMatrix(data = test_x, label = test_y)
      
      watchlist <- list(train=xgb_train, test=xgb_test)
      
      model = xgb.train(data = xgb_train, max.depth = 3, watchlist = watchlist, nrounds = 100)
      
      importance <- xgb.importance(colnames(xgb_train), model = model)
      
      importance
      
    }
  } 
  else {
    print("Error, define data correctly. If you think you have a bug, contact the author
        at mateusz.jakiel98@gmail.com")
  }
}
