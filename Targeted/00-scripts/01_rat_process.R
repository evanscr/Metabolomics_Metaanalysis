library(stringr)
library(janitor)
source("./functions.R")
################################################################################
# static inputs #
################################################################################
#folder paths
parent_folder <- "../"
data_folder <- "01-raw_data/01-MotrPAC_rat/"
table_output_folder <- "02-processed_data/00-tables/"
output_folder <- "02-processed_data/01-MotrPAC_rat/"
named_alignment_file <- "MotrPAC_refmet_names.csv"
pid_key_file <- "QC_EXPORT_DMAQC_Transfer_PASS_1A.6M_0.01_3-Data_Sets_QC_EXPORT_TXFR0.01_1A.6M_DS_MoTrPAC.PASS_Animal.BICLabelData_20190613.csv"
animal_key_file <- "QC_EXPORT_DMAQC_Transfer_PASS_1A.6M_0.01_3-Data_Sets_QC_EXPORT_TXFR0.01_1A.6M_DS_MoTrPAC.PASS_Animal.Key.csv"
factors_file <- "01-PASS1A-FACTORS.csv"
file_names <- list(rppos=list(results="metab_u_rppos/motrpac_results_metab_u_rppos_T31_plasma_named.txt",
                              metab="metab_u_rppos/motrpac_metadata_metabolites_metab_u_rppos_T31_plasma_named.txt"),
                   rpneg=list(results="metab_u_rpneg/motrpac_results_metab_u_rpneg_T31_plasma_named.txt",
                              metab="metab_u_rpneg/motrpac_metadata_metabolites_metab_u_rpneg_T31_plasma_named.txt"),
                   hilicpos=list(results="metab_u_hilicpos/motrpac_results_metab_u_hilicpos_T31_plasma_named.txt",
                                 metab="metab_u_hilicpos/motrpac_metadata_metabolites_metab_u_hilicpos_T31_plasma_named.txt"),
                   ionpneg=list(results="metab_u_ionpneg/motrpac_results_metab_u_ionpneg_T31_plasma_named.txt",
                                metab="metab_u_ionpneg/motrpac_metadata_metabolites_metab_u_ionpneg_T31_plasma_named.txt"))
# file_names <- list(rppos=list(results="metab_u_rppos/motrpac_results_metab_u_rppos_T31_plasma_named.txt",
#                               metab="metab_u_rppos/motrpac_metadata_metabolites_metab_u_rppos_T31_plasma_named.txt"),
#                    rpneg=list(results="metab_u_rpneg/motrpac_results_metab_u_rpneg_T31_plasma_named.txt",
#                               metab="metab_u_rpneg/motrpac_metadata_metabolites_metab_u_rpneg_T31_plasma_named.txt"),
#                    hilicpos=list(results="metab_u_hilicpos/motrpac_results_metab_u_hilicpos_T31_plasma_named.txt",
#                                  metab="metab_u_hilicpos/motrpac_metadata_metabolites_metab_u_hilicpos_T31_plasma_named.txt"),
#                    ionpneg=list(results="metab_u_ionpneg/motrpac_results_metab_u_ionpneg_T31_plasma_named.txt",
#                                 metab="metab_u_ionpneg/motrpac_metadata_metabolites_metab_u_ionpneg_T31_plasma_named.txt"),
#                    lrppos=list(results="metab_u_lrppos/motrpac_results_metab_u_lrppos_T31_plasma_named.txt",
#                                metab="metab_u_lrppos/motrpac_metadata_metabolites_metab_u_lrppos_T31_plasma_named.txt"),
#                    lrpneg=list(results="metab_u_lrpneg/motrpac_results_metab_u_lrpneg_T31_plasma_named.txt",
#                                metab="metab_u_lrpneg/motrpac_metadata_metabolites_metab_u_lrpneg_T31_plasma_named.txt"))

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


tmp <- read.csv(paste0(parent_folder, data_folder, named_alignment_file))
named_alignment <- vector(mode="list", length = length(file_names))
names(named_alignment) <- names(file_names)
counter=1
tryCatch({
  for(i in names(named_alignment)){
    tmp2 <- tmp[, c(counter, counter+1)]
    colnames(tmp2) <- c("original", "refmet_name")
    tmp3 <- cbind(data.frame(original = tmp2$original),
                  data.frame(clean_metabolites = make_clean_names(tmp2$original)),
                  data.frame(refmet_name = tmp2$refmet))
    tmp3 <- tmp3[tmp3$original != "", ]
    named_alignment[[i]] <- tmp3
    counter = counter+2
  }
},
error=function(e){
  message("column does not exist")
  message(conditionMessage(e))})
rm(tmp, counter, i, tmp2, tmp3)
################################################################################
# Process study #
################################################################################

processed_dat <- vector("list", length(names(file_names)))
names(processed_dat) <- names(file_names)
for(filename in names(file_names)){
  
  dat <- read.table(paste0(parent_folder, data_folder, file_names[[filename]][["results"]]),header = FALSE,
                    sep = "\t",quote = "",fill = FALSE)
  colnames(dat) <- dat[1,]
  rownames(dat) <- make_clean_names(dat$metabolite_name)
  dat <- dat[-1,-1]
  dat <- dat[,!grepl("QC", colnames(dat))]
  dat <- dat[,!grepl("CS", colnames(dat))]
  dat <- dat[,!grepl("SED", colnames(dat))]
  dat <- dat[,!grepl("Sample", colnames(dat))]
  dat <- dat[,!grepl("IPE", colnames(dat))]
  feature_metadat <- read.table(paste0(parent_folder, data_folder, file_names[[filename]][["metab"]]),header = FALSE,
                                sep = "\t",quote = "",fill = FALSE)
  colnames(feature_metadat) <- as.character(feature_metadat[1,])
  feature_metadat <- feature_metadat[-1,]
  feature_metadat$clean_metabolites <- make_clean_names(feature_metadat$metabolite_name)
  feature_metadat <- feature_metadat[, c(1, ncol(feature_metadat), 2:(ncol(feature_metadat)-1))]
  rownames(feature_metadat) <- feature_metadat$clean_metabolites
  feature_metadat <- feature_metadat[rownames(dat),]
  
  tmp_alignment <- named_alignment[[filename]]
  tmp_alignment$original <- make_clean_names(tmp_alignment$original)
  feature_metadat <- merge(tmp_alignment, feature_metadat, 
                           by.x = "original", by.y = "clean_metabolites", 
                           all.x = FALSE, all.y = TRUE)
  #feature_metadat$clean_metabolites <- make_clean_names(feature_metadat$metabolite_name)
  
  dat <- convert2int(dat)
  
  moTrPAC_factors <- read.csv(paste0(parent_folder, data_folder, factors_file))[, c(2:5, 10, 16)]
  moTrPAC_factors$pid <- as.character(moTrPAC_factors$pid)
  moTrPAC_factors <- moTrPAC_factors[match(unique(moTrPAC_factors$pid), moTrPAC_factors$pid), ]
  
  moTrPAC_PID_key <- read.csv(file = paste0(parent_folder, data_folder, pid_key_file),
                              header = TRUE)[,seq(3)]
  for(i in seq(ncol(moTrPAC_PID_key))){
    moTrPAC_PID_key[,i] <- as.character(moTrPAC_PID_key[,i])
  }
  rownames(moTrPAC_PID_key) <- moTrPAC_PID_key$vialLabel
  moTrPAC_PID_key <- moTrPAC_PID_key[colnames(dat),]
  sample_metadata <- merge(moTrPAC_factors, moTrPAC_PID_key, by="pid", all=TRUE)
  rownames(sample_metadata) <- sample_metadata$vialLabel
  sample_metadata <- sample_metadata[colnames(dat),]
  
  sample_metadata$timepoint[sample_metadata$timepoint == "0_hr"] <- "Pre"
  sample_metadata$timepoint[sample_metadata$timepoint == "IPE"] <- "IPE"
  sample_metadata$timepoint[sample_metadata$timepoint == "0.5_hr"] <- "30min"
  sample_metadata$timepoint[sample_metadata$timepoint == "1_hr"] <- "60min"
  sample_metadata$timepoint[sample_metadata$timepoint == "4_hr"] <- "240min"
  sample_metadata$timepoint[sample_metadata$timepoint == "7_hr"] <- "420min"
  sample_metadata$timepoint[sample_metadata$timepoint == "24_hr"] <- "1440min"
  sample_metadata$timepoint[sample_metadata$timepoint == "48_hr"] <- "2880min"
  sample_metadata$timepoint <- factor(sample_metadata$timepoint, levels = c("Pre", "IPE", "30min", 
                                                                            "60min", "240min", "420min", 
                                                                            "1440min", "2880min"))
  sample_metadata$timepoint2 <- sample_metadata$timepoint
  sample_metadata$timepoint2[!(sample_metadata$timepoint2 %in% c("Pre", "IPE"))] <- NA
  sample_metadata$timepoint2 <- factor(sample_metadata$timepoint2, levels = c("Pre", "IPE"))
  sample_metadata$study <- "Motrpac_rat"
  sample_metadata$platform <- filename
  sample_metadata <- sample_metadata[, c(1,8,6,2:5,9,10,11)]
  colnames(sample_metadata) <- c("subject", "sample", "age", "sex", 
                                 "weight", "group", "timepoint",
                                 "timepoint2", "study", "platform")
  sample_metadata$age <- NA
  sample_metadata$height <- NA
  sample_metadata$BMI <- NA
  sample_metadata$species <- "rat"
  sample_metadata$dataset <- paste0(unique(sample_metadata$study),
                                    ":", filename)
  sample_metadata$exhaustion <- "submax"
  sample_metadata$has_metadata <- 0
  rownames(sample_metadata) <- sample_metadata$sample
  
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
  
  write.csv(dat, paste0(parent_folder, output_folder, filename, "_MotrPAC_fitered_imputed_log_dat.csv"))
  write.csv(feature_metadat, paste0(parent_folder, output_folder, filename, "_MotrPAC_feature_metadata.csv"))
  
  refmet_metabolites <- feature_metadat[!is.na(feature_metadat$refmet), ]
  refmet_metabolites <- refmet_metabolites[refmet_metabolites$refmet != "", ]
  if(nrow(refmet_metabolites) > 0){
    refmet_dat <- cbind(refmet_metabolites[, c("original", "refmet_name")],
                        dat[refmet_metabolites$original, ])
    rownames(refmet_dat) <- make_clean_names(refmet_dat$refmet_name)
    refmet_dat <- t(refmet_dat[,-seq(2)])
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
  rm(dat, sample_metadata, feature_metadat, tmp_alignment, refmet_dat, refmet_metabolites)
}

rm(list=ls()[-match(c("parent_folder", "table_output_folder","processed_dat"), ls())])
save(processed_dat, file = paste0(parent_folder, table_output_folder, "01-Motrpac_rat.Rdata"))
