library(stringr)
library(janitor)
library(readxl)
source("./functions.R")
################################################################################
# static inputs #
################################################################################
#folder paths
parent_folder <- "../"
data_folder <- "01-raw_data/06-Sato/"
table_output_folder <- "02-processed_data/00-tables/"
output_folder <- "02-processed_data/06-Sato/"
file_names <- list(rppos="6-Sato_Plasma_Metabolomics Data_ALLFACTORSINCLUDED.xlsx")
sheet_name <- list(rppos="Sheet1")
#flags
missingness = 0.3
adj_mode = "unadjusted"
#can be minimum or median
data_impute = "minimum"
impute_min_scaler = 0.2
log_transform = TRUE

################################################################################
# pre-process #
################################################################################
##CREATE OUTPUT FOLDER IF IT DOESN'T EXIST
tryCatch({
  dir.create(paste0(parent_folder, output_folder))
  message(paste0("Creating ",paste0(parent_folder, output_folder," folder!")))
},
warning = function(w){
  message("Could not create Folder!")
  message(conditionMessage(w))})

################################################################################
# Process study #
################################################################################
processed_dat <- vector("list", length(names(file_names)))
names(processed_dat) <- names(file_names)
for(filename in names(file_names)){
  
  dat <- as.data.frame(read_xlsx(paste0(parent_folder, data_folder, file_names[[filename]]),
                                 sheet = sheet_name[[filename]], col_names = FALSE))
  rownames(dat) <- dat[, 6]
  
  sample_metadata <- dat[, c(5, 6, 3, 2, 4, 13)]
  sample_metadata$age <- NA
  # sample_metadata <- sample_metadata[, c(1, 2, 7, 3, 5, 4, 6)]
  sample_metadata <- sample_metadata[, c(1, 2, 7, 3, 5, 4)]
  colnames(sample_metadata) <- c("subject", "sample", "age", "sex", 
                                 "phase", "group")
  sample_metadata <- sample_metadata[-1, ]
  rownames(sample_metadata) <- sample_metadata$sample
  
  sample_metadata$timepoint <- ifelse(sample_metadata$group == "Exercise", "IPE", "Pre")
  sample_metadata$timepoint <- factor(sample_metadata$timepoint, levels = c("Pre", "IPE"))
  sample_metadata$timepoint2 <- sample_metadata$timepoint
  
  sample_metadata$age <- NA
  sample_metadata$weight <- NA
  sample_metadata$BMI <- NA
  sample_metadata$study <- "Sato"
  sample_metadata$platform <- filename
  sample_metadata$species <- "rat"
  sample_metadata$dataset <- paste0(unique(sample_metadata$study),
                                    ":", filename)
  sample_metadata$exhaustion <- "submax"
  sample_metadata$has_metadata <- 0
  rownames(sample_metadata) <- sample_metadata$sample
  
  dat <- dat[, -seq(16)]
  metab_names <- unlist(dat[1,])
  refmet_names <- metab_names
  clean_metabs <- make_clean_names(metab_names)
  
  colnames(dat) <- clean_metabs
  dat <- dat[-1, ]
  
  dat <- t(dat)
  dat <- convert2int(dat)
  unkowns <- rownames(dat)[grepl("^[Xx]_[0-9]{5}", rownames(dat))]
  dat <- dat[-match(unkowns, rownames(dat)), ]
  
  feature_metadat <- data.frame(original = metab_names,
                                clean_metabolites = clean_metabs,
                                refmet_name = refmet_names,
                                row.names = clean_metabs)
  feature_metadat <- feature_metadat[rownames(dat),]
  
  
  
  #missingness filter
  if(!is.na(missingness)){
    ##apply missingness filter at 70%
    filter_70p <- apply(dat, 1, function(x) ifelse(sum(is.na(x)) < ceiling(ncol(dat)*missingness), TRUE, FALSE))
    dat <- dat[filter_70p,]
    feature_metadat <- feature_metadat[filter_70p,]
    #clear space
    rm(filter_70p)
  }
  
  #data imputation
  if(!is.na(data_impute)){
    dat <- imputeMissing(dat=dat,
                         method = "minimum",
                         fraction = impute_min_scaler)
  }
  
  #log transform
  if(log_transform){
    dat <- log2(dat)
  }
  
  write.csv(dat, paste0(parent_folder, output_folder, filename, "_FUEL_fitered_imputed_log_dat.csv"))
  write.csv(feature_metadat, paste0(parent_folder, output_folder, filename, "_FUEL_feature_metadata.csv"))
  
  refmet_metabolites <- feature_metadat[!is.na(feature_metadat$refmet), ]
  refmet_metabolites <- refmet_metabolites[refmet_metabolites$refmet != "", ]
  if(nrow(refmet_metabolites) > 0){
    refmet_dat <- cbind(refmet_metabolites[, c("original", "refmet_name")],
                        dat[refmet_metabolites$clean_metabolites, ])
    rownames(refmet_dat) <- make_clean_names(refmet_dat$refmet_name)
    refmet_dat <- t(refmet_dat[,-seq(2)])
    # refmet_dat <- t(dat[refmet_metabolites$clean_metabolites,])
    refmet_dat <- merge(sample_metadata, refmet_dat, by.x = "sample", by.y="row.names", all=TRUE)
    # rownames(refmet_dat) <- refmet_dat$sample
  }else{
    refmet_dat <- NULL
  }
  #add to output list
  processed_dat[[filename]] <- list(data=dat,
                                    refmet_dat=refmet_dat,
                                    samples=sample_metadata,
                                    features=feature_metadat)
  
  #clear
  rm(dat, sample_metadata, feature_metadat)
}

rm(list=ls()[-match(c("parent_folder", "table_output_folder","processed_dat"), ls())])
save(processed_dat, file = paste0(parent_folder, table_output_folder, "06-Sato.Rdata"))