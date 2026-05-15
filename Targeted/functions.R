convert2int <- function(dat){
  
  rows <- rownames(dat)
  cols <- colnames(dat)
  tryCatch({
    output <- apply(dat, MARGIN = 2, function(x) as.numeric(x))
  })
  rownames(output) <- rows
  colnames(output) <- cols
  
  stopifnot(all(rownames(dat) == rownames(output)))
  stopifnot(all(colnames(dat) == colnames(output)))
  
  return(output)
}
row_min <- function(row){
  
  output <- min(row[!is.na(row)])
  
  if(is.na(output))stop("Output should not be NA, there was a problem with input!")
  return(output)
}

row_median <- function(row){
  
  output <- median(row[!is.na(row)])
  
  if(is.na(output))stop("Output should not be NA, there was a problem with input!")
  return(output)
}
imputeMissing <- function(dat,
                          fraction = 0.2,
                          method = c("minimum", "median")){
  method = match.arg(method)
  if(method == "minimum"){
    impute_values <- apply(dat, MARGIN = 1, row_min) * fraction
  }else if(data_impute == "median"){
    impute_values <- apply(dat, MARGIN = 1, row_median)
  }
  for(i in seq(nrow(dat))){
    
    metab = rownames(dat)[i]
    dat[metab,][is.na(dat[metab,])] <- impute_values[metab]
    rm(metab)
  }
  #clear space
  rm(i)
  return(dat)
}

metabDE <- function(dat,
                    method=c("ttest", "lm"),
                    paired=FALSE,
                    feature_dat,
                    sample_dat,
                    group_label_col,
                    comparison,
                    ref){
  method<-match.arg(method)
  
  stopifnot(!(group_label_col %in% colnames(sample_dat)))
  stopifnot(!(c(ref, comparison) %in% sample_dat[["group_label_col"]]))
  
  ref_samples <- rownames(sample_dat)[sample_dat[["group_label_col"]] == ref]
  comparison_samples <- rownames(sample_dat)[sample_dat[["group_label_col"]] == comparison]
  stopifnot(all(c(ref_samples, comparison_samples) %in% rownames(dat)))
  
  ref_group <- dat[ref_samples, ]
  comparison_group <- dat[comparison_group, ]
  
  if(!all(colnames(ref_group) == colnames(comparison_group))){
    stop("metabolites across group1 and group2 are not in same order!")
  }
  
  if(paired){
    if(!all(rownames(ref_group) == rownames(comparison_group))){
      stop("samples across group1 and group2 are not in same order!")
  }
  }
  p <- ncol(group1)
  stopifnot(all(rownames(metabolite_metadata) == colnames(group1)))
  stopifnot(all(rownames(metabolite_metadata) == colnames(group2)))
  stopifnot(all(colnames(group1) == colnames(group2)))
  
  metab_info <- data.frame(metabolite_metadata)
  metab_info$comparison <- paste0(z, "vsPre")
  metab_info$foldchange <- colMeans(group1) - colMeans(group2)
  metab_info$fcdirection <- sapply(1:p, function(i) ifelse(metab_info$foldchange[i] > 0, "Up", "Down"))
  metab_info$statistic <- sapply(1:p, function(i) t.test(group1[,i], group2[,i], paired = TRUE, var.equal = FALSE)$statistic)
  metab_info$pvalue <- sapply(1:p, function(i) t.test(group1[,i], group2[,i], paired = TRUE, var.equal = FALSE)$p.value)
  metab_info$qvalue <- p.adjust(metab_info$pvalue, "BH")
  metab_info$DEstatus <- sapply(1:p, function(i) ifelse(abs(metab_info$qvalue[i]) >= 0.05, FALSE, TRUE))
  
  ##order metabolites by padj
  metab_info <- metab_info[order(metab_info$qvalue),]
}
  

