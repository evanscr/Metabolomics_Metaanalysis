library(stringr)
library(janitor)
source("./functions.R")
################################################################################
# flags #
################################################################################
#folder paths
output_folder <- "02-processed_data/08-TrainSed_LCRHCR_rats/"
parent_folder <- "../"
data_folder <- "../01-raw_data/08-TrainSed_LCRHCR_rats/"
table_output_folder <- "02-processed_data/00-tables/"
file_names <- c("LCRHCR_TrainSed_Plasma_RPLCpos_raw_data_with_manifest_and_factor_info_QCaNVaS.csv",
                "LCRHCR_TrainSed_Plasma_RPLCneg_raw_data_with_manifest_and_factor_info_QCaNVaS.csv",
                "LCRHCR_TrainSed_Plasma_IPCneg_raw_data_with_manifest_and_factor_info_QCaNVaS.csv")
names(file_names) <- c("rppos", "rpneg", "ionpneg")
#flags
missingness = 0.3
adj_mode = "unadjusted"
#can be minimum or median
data_impute = "minimum"
impute_min_scaler = 0.2
log_transform = TRUE
################################################################################

##CREATE OUTPUT FOLDER IF IT DOESN'T EXIST
tryCatch({
  dir.create(paste0(parent_folder, output_folder))
  message(paste0("Creating ",paste0(parent_folder, output_folder," folder!")))
},
warning = function(w){
  message("Could not create Folder!")
  message(conditionMessage(w))})



processed_dat <- vector("list", length(names(file_names)))
names(processed_dat) <- names(file_names)
for(i in names(file_names)){
  
  dat <- read.csv(paste0(data_folder, file_names[[i]]),header = FALSE)
  
  sample_names <- unlist(dat[-1, 2])
  sample_names <- sample_names[!grepl("QC", sample_names)]
  sample_names <- sample_names[!grepl("blank", sample_names)]
  sample_names <- sample_names[!grepl("CHEAR", sample_names)]
  sample_names <- sample_names[sample_names != ""]
  sample_metadata <- dat[,seq(12)]
  colnames(sample_metadata) <- sample_metadata[1, ]
  rownames(sample_metadata) <- sample_metadata[, 2]
  sample_metadata <- sample_metadata[-1, ]
  
  sample_metadata <- sample_metadata[sample_names, ]
  sample_metadata$age <- NA
  sample_metadata$weight <- NA
  sample_metadata$height <- NA
  sample_metadata$BMI <- NA
  sample_metadata$subject <- sample_metadata$SampleName
  sample_metadata <- sample_metadata[, c(17, 2, 13, 7, 14, 15, 16, 9, 8)]
  colnames(sample_metadata) <- c("subject", "sample", "age", "sex", "weight",
                                 "height", "BMI", "group", "timepoint")
  sample_metadata$timepoint[sample_metadata$timepoint == "REST"] <- "Pre"
  sample_metadata$timepoint[sample_metadata$timepoint == "MAX"] <- "IPE"
  sample_metadata$timepoint[sample_metadata$timepoint == "HALF MAX"] <- "half_max_during"
  sample_metadata$timepoint[sample_metadata$timepoint == "REC"] <- "60min"
  sample_metadata$timepoint <- factor(sample_metadata$timepoint, levels = c("Pre", "half_max_during",
                                                                            "IPE", "60min"))
  sample_metadata$timepoint2 <- sample_metadata$timepoint
  sample_metadata$timepoint2[!(sample_metadata$timepoint2 %in% c("Pre", "IPE"))] <- NA
  sample_metadata$timepoint2 <- factor(sample_metadata$timepoint2, levels = c("Pre", "IPE"))
  sample_metadata$study <- "LCR_HCR"
  sample_metadata$platform <- i
  sample_metadata$species <- "rat"
  sample_metadata$dataset <- paste0(unique(sample_metadata$study),
                                    ":", i)
  sample_metadata$exhaustion <- "max"
  sample_metadata$has_metadata <- 0
  rownames(sample_metadata) <- sample_metadata$sample
  
  metab_names <- unlist(dat[1, seq(13, ncol(dat))])
  clean_metabs <- make_clean_names(metab_names)
  refmet_names <- metab_names
  colnames(dat) <- make_clean_names(dat[1,])
  rownames(dat) <- dat[,2]
  dat <- dat[-1, -seq(12)]

  dat <- dat[sample_names, ]
  
  feature_metadat <- data.frame(original = metab_names,
                                clean_metabolites = clean_metabs,
                                refmet_name = refmet_names,
                                row.names = clean_metabs)

  dat <- t(dat)
  dat <- convert2int(dat)
  
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
  
  write.csv(dat, paste0(parent_folder, output_folder, "LCRHCR_", i, "_log_dat.csv"))
  write.csv(sample_metadata, paste0(parent_folder, output_folder, "LCRHCR_", i, "_sample_metadata.csv"))
  write.csv(feature_metadat, paste0(parent_folder, output_folder, "LCRHCR_", i, "_feature_metadata.csv"))
  
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
  processed_dat[[i]] <- list(data=dat,
                                    refmet_dat=refmet_dat,
                                    samples=sample_metadata,
                                    features=feature_metadat)
  
  #clear
  rm(dat, sample_metadata, feature_metadat)
 
}

rm(list=ls()[-match(c("parent_folder", "table_output_folder","processed_dat"), ls())])
save(processed_dat, file = paste0(parent_folder, table_output_folder, "08-LCR_HCR.Rdata"))
