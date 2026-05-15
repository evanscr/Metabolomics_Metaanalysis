library(stringr)
library(janitor)
source("./functions.R")
################################################################################
# static inputs #
################################################################################
#folder paths
parent_folder <- "../"
data_folder <- "01-raw_data/03-FUEL/"
table_output_folder <- "02-processed_data/00-tables/"
output_folder <- "02-processed_data//03-FUEL/"
named_alignment_file <- "FUEL_refmet_names.csv"
factors_file <- "03-FUEL-FACTORS.csv"
file_names <- list(rppos="FUEL_Plasma_MetabolomicsData_v2_RPLCpos.csv",
                   hilicneg="FUEL_Plasma_MetabolomicsData_v2_HILICneg.csv")

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
  
  dat <- read.csv(paste0(parent_folder, data_folder, file_names[[filename]]), header = FALSE)
  metab_names <- unlist(dat[1, 5:ncol(dat)])
  colnames(dat) <- dat[1,]
  rownames(dat) <- dat[, 1]
  dat <- dat[-1, ]
  sample_metadata <- dat[, c(2, 1, 4)]
  colnames(sample_metadata) <- c("subject", "sample", "timepoint")
  
  dat <- dat[, -seq(4)]
  colnames(dat) <- make_clean_names(colnames(dat))
  dat <- t(dat)
  dat <- convert2int(dat)
  
  FUEL_factors <- read.csv(paste0(parent_folder, data_folder, factors_file))[, c(seq(6), 10, 11)]
  FUEL_factors <- FUEL_factors[, c(1,4,3,5,6,7,8,2)]
  colnames(FUEL_factors) <- c("subject", "age", "sex", "race", "ethnicity",
                              "weight", "height", "group")
  sample_metadata <- merge(sample_metadata, FUEL_factors, by="subject", all.x = TRUE, all.y = FALSE)
  sample_metadata <- sample_metadata[, c(1,2,4, 5, 6, 7, 8, 9, 10, 3)]
  
  sample_metadata$timepoint[sample_metadata$timepoint == "pre"] <- "Pre"
  sample_metadata$timepoint[sample_metadata$timepoint == "IPE"] <- "IPE"
  sample_metadata$timepoint[sample_metadata$timepoint == "Post-3"] <- "3min"
  sample_metadata$timepoint[sample_metadata$timepoint == "Post-6"] <- "6min"
  sample_metadata$timepoint[sample_metadata$timepoint == "3"] <- "3min_during"
  sample_metadata$timepoint[sample_metadata$timepoint == "6"] <- "6min_during"
  sample_metadata$timepoint[sample_metadata$timepoint == "9"] <- "9min_during"
  sample_metadata$timepoint[sample_metadata$timepoint == "12"] <- "12min_during"
  sample_metadata$timepoint[sample_metadata$timepoint == "15"] <- "15min_during"
  sample_metadata$timepoint[sample_metadata$timepoint == "18"] <- "18min_during"
  sample_metadata$timepoint[sample_metadata$timepoint == "21"] <- "21min_during"
  sample_metadata$timepoint[sample_metadata$timepoint == "24"] <- "24min_during"
  sample_metadata$timepoint[sample_metadata$timepoint == "27"] <- "27min_during"
  sample_metadata$timepoint[sample_metadata$timepoint == "30"] <- "30min_during"
  sample_metadata$timepoint[sample_metadata$timepoint == "33"] <- "33min_during"
  sample_metadata$timepoint[sample_metadata$timepoint == "36"] <- "36min_during"
  sample_metadata$timepoint[sample_metadata$timepoint == "42"] <- "42min_during"
  sample_metadata$timepoint[sample_metadata$timepoint == "45"] <- "45min_during"
  sample_metadata$timepoint[sample_metadata$timepoint == "-3"] <- "3min_Pre"
  
  sample_metadata$timepoint <- factor(sample_metadata$timepoint, levels = c("Pre", "IPE","3min", "6min",
                                                                            "3min_during", "6min_during",
                                                                            "9min_during", "12min_during",
                                                                            "15min_during", "18min_during",
                                                                            "21min_during", "24min_during",
                                                                            "27min_during", "30min_during",
                                                                            "33min_during", "36min_during",
                                                                            "42min_during", "45min_during",
                                                                            "3min_Pre"))
  sample_metadata$timepoint2 <- sample_metadata$timepoint
  sample_metadata$timepoint2[!(sample_metadata$timepoint2 %in% c("Pre", "IPE"))] <- NA
  sample_metadata$timepoint2 <- factor(sample_metadata$timepoint2, levels = c("Pre", "IPE"))
  sample_metadata$study <- "FUEL"
  sample_metadata$platform <- filename
  sample_metadata$BMI <- sample_metadata$weight / (sample_metadata$height/100)^2
  sample_metadata$species <- "human"
  sample_metadata$dataset <- paste0(unique(sample_metadata$study),
                                    ":", filename)
  sample_metadata$exhaustion <- "max"
  sample_metadata$has_metadata <- 1
  rownames(sample_metadata) <- sample_metadata$sample
  
  feature_metadat <- data.frame(original = metab_names,
                                clean_metabolites = make_clean_names(metab_names))
  feature_metadat <- merge(feature_metadat,
                           named_alignment[[filename]][, c("original", "refmet_name")],
                           by = "original", all = TRUE)
  rownames(feature_metadat) <- feature_metadat$clean_metabolites
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
  rm(dat, sample_metadata, feature_metadat, FUEL_factors)
}

rm(list=ls()[-match(c("parent_folder", "table_output_folder","processed_dat"), ls())])
save(processed_dat, file = paste0(parent_folder, table_output_folder, "03-FUEL.Rdata"))




