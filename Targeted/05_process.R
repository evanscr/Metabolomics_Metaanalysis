library(stringr)
library(janitor)
library(readxl)
source("./functions.R")
################################################################################
# static inputs #
################################################################################
#folder paths
parent_folder <- "../"
data_folder <- "01-raw_data/05-Contrepois/"
table_output_folder <- "02-processed_data/00-tables/"
output_folder <- "02-processed_data/05-Contrepois/"
factors_file <- "05-Contrepois_FACTORS.csv"
file_names <- list(rppos="Contrepois et al. 2020 Exercise_Metabolomic.xlsx")
sheet_name <- list("rppos" = "Contrepois et al. 2020 Exercise")
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
  
  dat <- as.data.frame(read_xlsx(path = paste0(parent_folder, data_folder, file_names[[filename]]),
                                 sheet = sheet_name[[filename]]))
  
  #feature metadata
  metab_names <- colnames(dat)[-seq(3)]
  clean_metabs <- make_clean_names(metab_names)
  refmet_names <- dat[1, -seq(3)]
  
  #sample metadata
  sample_names <- dat[-1, 1]
  sample_time <- dat[-1, 2]
  
  #clean up dat
  dat <- dat[-1, -seq(3)]
  colnames(dat) <- clean_metabs
  rownames(dat) <- sample_names
  dat <- t(dat)
  dat <- convert2int(dat)
  unkowns <- rownames(dat)[grepl("^[cC][0-9]+h[0-9]+", rownames(dat))]
  dat <- dat[-match(unkowns, rownames(dat)), ]
  
  sample_metadata <- data.frame(subject = vapply(sample_names,
                                                 function(x) unlist(strsplit(x, "-"))[1],
                                                 FUN.VALUE = character(1)),
                                sample = sample_names,
                                group = NA,
                                timepoint = sample_time)
  colnames(sample_metadata) <- c("subject", "sample","group", "timepoint")
  
 
  
  
  sample_factors <- read.csv(paste0(parent_folder, data_folder, factors_file))
  sample_factors <- sample_factors[, c(1,5, 4, 6)]
  colnames(sample_factors) <- c("subject", "age", "sex", "BMI")
  sample_metadata <- merge(sample_metadata, sample_factors, by="subject", 
                           all.x = TRUE, all.y = FALSE)
  sample_metadata <- sample_metadata[, c(1,2, 5, 6, 7, 3, 4)]
  
  sample_metadata$timepoint[sample_metadata$timepoint == "Pre"] <- "Pre"
  sample_metadata$timepoint[sample_metadata$timepoint == "2min"] <- "IPE"
  sample_metadata$timepoint[sample_metadata$timepoint == "15min"] <- "15min"
  sample_metadata$timepoint[sample_metadata$timepoint == "30min"] <- "30min"
  sample_metadata$timepoint[sample_metadata$timepoint == "60min"] <- "60min"
  sample_metadata$timepoint <- factor(sample_metadata$timepoint, levels = c("Pre", "IPE", "15min",
                                                                            "30min", "60min"))
  sample_metadata$timepoint2 <- sample_metadata$timepoint
  sample_metadata$timepoint2[!(sample_metadata$timepoint2 %in% c("Pre", "IPE"))] <- NA
  sample_metadata$timepoint2 <- factor(sample_metadata$timepoint2, levels = c("Pre", "IPE"))

  sample_metadata$weight <- NA
  sample_metadata$height <- NA
  sample_metadata$study <- "Contrepois"
  sample_metadata$platform <- filename
  sample_metadata$species <- "human"
  sample_metadata$dataset <- paste0(unique(sample_metadata$study),
                                    ":", filename)
  sample_metadata$exhaustion <- "max"
  sample_metadata$has_metadata <- 1
  rownames(sample_metadata) <- sample_metadata$sample
  
  feature_metadat <- data.frame(original = unlist(metab_names),
                                clean_metabolites = unlist(clean_metabs),
                                refmet_name = unlist(refmet_names))

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
  
  write.csv(dat, paste0(parent_folder, output_folder, filename, "_Contrepois_fitered_imputed_log_dat.csv"))
  write.csv(feature_metadat, paste0(parent_folder, output_folder, filename, "_Contrepois_feature_metadata.csv"))
  
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
  rm(dat, sample_metadata, feature_metadat, sample_factors)
}

rm(list=ls()[-match(c("parent_folder", "table_output_folder","processed_dat"), ls())])
save(processed_dat, file = paste0(parent_folder, table_output_folder, "05-Contrepois.Rdata"))



