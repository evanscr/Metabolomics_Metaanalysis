library(stringr)
library(janitor)
source("./functions.R")
################################################################################
# flags #
################################################################################
#folder paths
parent_folder <- "../"
data_folder <- "01-raw_data/02-MotrPAC_human/"
table_output_folder <- "02-processed_data/00-tables/"
output_folder <- "02-processed_data/02-MotrPAC_human/"
named_alignment_file <- "MotrPAC_refmet_names.csv"
factors_file <- "02-MoTrPAC Precovid Human FACTORS.csv"
file_names <- "motrpac-preCAWG-human-plasma-normalized-log.RData"
exercise_type = list("endurance"="ADUEndur")

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

tryCatch({
  dir.create(paste0(parent_folder, output_folder, names(exercise_type)[1]))
  message(paste0("Creating ",paste0(parent_folder, output_folder, 
                                    names(exercise_type)[1], " folder!")))
},
warning = function(w){
  message("Could not create Folder!")
  message(conditionMessage(w))})

sample_metadata <- read.csv(paste0(parent_folder, data_folder, factors_file))
sample_metadata$study2 <- "Motrpac_human"
sample_metadata$platform <- NA
sample_metadata$BMI <- sample_metadata$weight / ((sample_metadata$height/100)^2)
sample_metadata$species <- "human"
sample_metadata$timepoint[sample_metadata$timepoint == "pre_exercise"] <- "Pre"
sample_metadata$timepoint[sample_metadata$timepoint == "post_10_min"] <- "IPE"
sample_metadata$timepoint[sample_metadata$timepoint == "post_30_min"] <- "30min"
sample_metadata$timepoint[sample_metadata$timepoint == "post_4_hr"] <- "240min"
sample_metadata$timepoint[sample_metadata$timepoint == "post_24_hr"] <- "1440min"
sample_metadata$timepoint[sample_metadata$timepoint == "post_10_min"] <- "IPE"
sample_metadata$timepoint[sample_metadata$timepoint == "during_20_min"] <- "20min_during"
sample_metadata$timepoint[sample_metadata$timepoint == "during_40_min"] <- "40min_during"
sample_metadata$timepoint <- factor(sample_metadata$timepoint, levels = c("Pre", "20min_during", 
                                                                          "40min_during", "IPE", 
                                                                          "30min", "240min", "1440min"))

sample_metadata$timepoint2 <- sample_metadata$timepoint
sample_metadata$timepoint2[!(sample_metadata$timepoint2 %in% c("Pre", "IPE"))] <- NA
sample_metadata$timepoint2 <- factor(sample_metadata$timepoint2, levels = c("Pre", "IPE"))
sample_metadata <- sample_metadata[, c(8, 1, 5, 4, 6, 7, 41, 2, 3, 43, 39, 40, 42)]
colnames(sample_metadata) <- c("subject", "sample", "age", "sex", "weight", 
                               "height", "BMI","group", "timepoint",
                               "timepoint2", "study", "platform", "species")
sample_metadata$dataset <- NA
sample_metadata$subject <- as.character(sample_metadata$subject)
sample_metadata$sample <- as.character(sample_metadata$sample)
rownames(sample_metadata) <- vapply(sample_metadata$sample, function(x) substr(x, 1, nchar(x)-4), character(1))
################################################################################
# Process study #
################################################################################

load(paste0(parent_folder, data_folder, file_names))
plasma.data <- plasma.data[c(13, 12, 8, 9, 11, 10)]
processed_dat <- vector("list", length = 6)
names(processed_dat) <- c("rppos", "rpneg", "hilicpos", 
                              "ionpneg", "lrppos", "lrpneg")
names(plasma.data) <- names(processed_dat)

for(filename in names(plasma.data)){
  tmpdat <- plasma.data[[filename]]
  
  feature_metadat <- tmpdat$row_annot
  feature_metadat <- feature_metadat[feature_metadat$is_named,]
  colnames(feature_metadat)[seq(6)] <- c("metabolite_name", "refmet_name", 
                                 "rt", "mz", "mass", "formula")
  feature_metadat <- cbind(data.frame(original = feature_metadat$metabolite_name,
                                      clean_metabolites = make_clean_names(feature_metadat$metabolite_name)),
                           feature_metadat[, -match("metabolite_name", colnames(feature_metadat))])
  
  tmp_sample_metadata <- sample_metadata
  tmp_sample_metadata$platform <- filename
  tmp_sample_metadata$dataset <- paste0(tmp_sample_metadata$study, ":", 
                                        tmp_sample_metadata$platform)
  tmp_sample_metadata$exhaustion <- "submax"
  tmp_sample_metadata$has_metadata <- 1
  tmp_sample_metadata <- tmp_sample_metadata[tmp_sample_metadata$group == exercise_type[1],]
  
  dat <- tmpdat$normalized_data$raw
  dat <- dat[feature_metadat$original, ]
  rownames(dat) <- feature_metadat$clean_metabolites
  colnames(dat) <- vapply(colnames(dat), function(x) substr(x, 1, nchar(x)-4), character(1))
  keep <- rownames(tmp_sample_metadata)[rownames(tmp_sample_metadata) %in% colnames(dat)]
  # tmp_sample_metadata <- tmp_sample_metadata[tmp_sample_metadata$sample %in% keep,]
  tmp_sample_metadata <- tmp_sample_metadata[keep, ]
  dat <- dat[ ,rownames(tmp_sample_metadata)]
  if(!all(colnames(dat) == rownames(tmp_sample_metadata))) stop("problem matching factors!")
  
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
  
  write.csv(dat, paste0(parent_folder, output_folder, names(exercise_type)[1],"/", 
                        filename,"_",names(exercise_type)[1], 
                        "_MotrPAC_human_fitered_imputed_log_dat.csv"))
  write.csv(feature_metadat, paste0(parent_folder, output_folder, names(exercise_type)[1],"/", 
                                    filename,"_",names(exercise_type)[1], 
                                    "_MotrPAC_human_feature_metadata.csv"))
  
  refmet_metabolites <- feature_metadat[!is.na(feature_metadat$refmet), ]
  refmet_metabolites <- refmet_metabolites[refmet_metabolites$refmet != "", ]
  rownames(refmet_metabolites) <- refmet_metabolites$clean_metabolites
  refmet_metabolites <- refmet_metabolites[rownames(dat), ]
  if(!all(rownames(dat) == rownames(refmet_metabolites))) stop("problem making refmet_dat!")
  if(nrow(refmet_metabolites) > 0){
    refmet_dat <- cbind(refmet_metabolites[, c("original", "refmet_name")],
                        dat[refmet_metabolites$clean_metabolites, ])
    rownames(refmet_dat) <- make_clean_names(refmet_dat$refmet_name)
    refmet_dat <- t(refmet_dat[,-seq(2)])
    refmet_dat <- merge(tmp_sample_metadata, refmet_dat, by = "row.names", all=TRUE)
    rownames(refmet_dat) <- refmet_dat$sample
  }else{
    refmet_dat <- NULL
  }
  
  processed_dat[[filename]] <- list(data=dat,
                                    refmet_dat=refmet_dat,
                                    samples=tmp_sample_metadata,
                                    features=feature_metadat)
  
  #clear
  rm(dat, tmp_sample_metadata, feature_metadat, refmet_dat, refmet_metabolites)
}


rm(list=ls()[-match(c("parent_folder", "table_output_folder", 
                      "exercise_type", "processed_dat"), ls())])
save(processed_dat, file = paste0(parent_folder, table_output_folder, 
                                  "02-Motrpac_human_", names(exercise_type)[1], ".Rdata"))
