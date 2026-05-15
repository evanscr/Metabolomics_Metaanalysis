library(stringr)
library(janitor)
library(readxl)
source("./functions.R")
################################################################################
# static inputs #
################################################################################
#folder paths
parent_folder <- "../"
data_folder <- "01-raw_data/07-ST001749/"
table_output_folder <- "02-processed_data/00-tables/"
output_folder <- "02-processed_data/07-ST001749/"
factors_file <- "ST001749-REACH_metabolomics_Metadata_ToBeCompleted.csv"
file_names <- list(rppos="REACH Metabolomics Study Reversed phase POSITIVE ION MODE_1.csv",
                   rpneg="REACH Metabolomics Study Reversed phase NEGATIVE ION MODE.csv",
                   hilicneg="REACH Metabolomics Study HILIC NEGATIVE ION MODE.csv")

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

sample_factors <- read.csv(paste0(parent_folder, data_folder, factors_file), header = TRUE)[, seq(3)]
sample_factors <- data.frame(subject = sample_factors$RefMet_name,
                             sample = sample_factors$RefMet_name,
                             age = NA,
                             sex = NA,
                             weight = NA,
                             height = NA,
                             BMI = NA,
                             group = sample_factors$Treatment,
                             timepoint = sample_factors$Timepoint,
                             row.names = sample_factors$RefMet_name)
################################################################################
# Process study #
################################################################################
processed_dat <- vector("list", length(names(file_names)))
names(processed_dat) <- names(file_names)
for(filename in names(file_names)){
  
  dat <- read.csv(paste0(parent_folder, data_folder, file_names[[filename]]), header = FALSE)
  metab_names <- unlist(dat[seq(3, nrow(dat)), 1])
  clean_metabs <- make_clean_names(metab_names)
  refmet_names <- unlist(dat[seq(3, nrow(dat)), 2])
  sample_names <- unlist(dat[1, seq(3, ncol(dat))])

  sample_metadata <- sample_factors[sample_names, ]
  
  sample_metadata$timepoint[sample_metadata$timepoint == "Pre"] <- "Pre"
  sample_metadata$timepoint[sample_metadata$timepoint == "Post"] <- "IPE"
  sample_metadata$timepoint <- factor(sample_metadata$timepoint, levels = c("Pre", "IPE"))
  sample_metadata$timepoint2 <- sample_metadata$timepoint
  
  sample_metadata$study <- "ST001749"
  sample_metadata$platform <- filename
  sample_metadata$species <- "human"
  sample_metadata$dataset <- paste0(unique(sample_metadata$study),
                                    ":", filename)
  sample_metadata$exhaustion <- "submax"
  sample_metadata$has_metadata <- 0
  rownames(sample_metadata) <- sample_metadata$sample
  
  colnames(dat) <- dat[1,]
  rownames(dat) <- make_clean_names(dat[,1])
  dat <- dat[-seq(2), -seq(2)]
  dat <- convert2int(dat)
  
  # unkowns <- rownames(dat)[grepl("^[Xx]_[0-9]{5}", rownames(dat))]
  # dat <- dat[-match(unkowns, rownames(dat)), ]
  
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
save(processed_dat, file = paste0(parent_folder, table_output_folder, "07-ST001749.Rdata"))