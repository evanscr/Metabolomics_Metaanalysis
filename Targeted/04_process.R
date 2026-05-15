library(stringr)
library(janitor)
source("./functions.R")
################################################################################
# flags #
################################################################################
#folder paths
parent_folder <- "../"
data_folder <- "01-raw_data/04-ST001789/"
table_output_folder <- "02-processed_data/00-tables/"
output_folder <- "02-processed_data/04-ST001789/"
factors_file <- "04-ST001789-FACTORS.csv"
file_names <- list(rppos="04-ST001790-RPLCdata.txt",
                   hilicneg="04-ST001790-HILICdata.txt")

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

sample_factors <- read.csv(paste0(parent_folder, data_folder, factors_file))
sample_factors <- sample_factors[, c(1, 3, 2, 7, 6, 8)]

colnames(sample_factors) <- c("subject", "age", "sex", "weight", "height", "BMI")
sample_factors$subject <- vapply(sample_factors$subject, 
                                 function(x) str_replace(x, "\\. ", ""),
                                 FUN.VALUE = character(1))
sample_factors$subject <- vapply(sample_factors$subject, 
                                 function(x) str_replace(x, "\\.  ", ""),
                                 FUN.VALUE = character(1))
sample_factors$subject <- vapply(sample_factors$subject, 
                                 function(x) str_replace(x, "\\.", ""),
                                 FUN.VALUE = character(1))
################################################################################
# Process study #
################################################################################
processed_dat <- vector("list", length(names(file_names)))
names(processed_dat) <- names(file_names)
for(filename in names(file_names)){
  
  dat <- read.table(paste0(parent_folder, data_folder, file_names[[filename]]),header = FALSE,
                    sep = "\t",quote = "",fill = FALSE)
  metab_names <- unlist(dat[3:nrow(dat), 1])
  refmet_names <- unlist(dat[3:nrow(dat), 2])
  colnames(dat) <- dat[1,]
  rownames(dat) <- make_clean_names(dat[, 1])
  dat <- dat[,-seq(2)]
  sample_metadata <- data.frame(subject = unlist(dat[1,]),
                                sample = unlist(dat[1,]),
                                group = NA,
                                timepoint = unlist(dat[2,]))
  sample_metadata$subject <- vapply(sample_metadata$subject,
                                    function(x) unlist(strsplit(x, "_"))[1],
                                    FUN.VALUE = character(1))
  sample_metadata$timepoint <- vapply(sample_metadata$timepoint, 
                                      function(x) str_replace(x, "Group:", ""), 
                                      FUN.VALUE = character(1))
  
  dat <- dat[-seq(2), ]
  dat <- convert2int(dat)
  
 
  sample_metadata <- merge(sample_metadata, sample_factors, by="subject", all = TRUE)
  sample_metadata <- sample_metadata[, c(1,2, 5, 6, 7, 8, 9, 3, 4)]
  sample_metadata$timepoint[sample_metadata$timepoint == "Pre"] <- "Pre"
  sample_metadata$timepoint[sample_metadata$timepoint == "Time 0"] <- "IPE"
  sample_metadata$timepoint[sample_metadata$timepoint == "Time 60"] <- "60min"
  sample_metadata$timepoint <- factor(sample_metadata$timepoint, levels = c("Pre", "IPE", "60min"))
  
  sample_metadata$timepoint2 <- sample_metadata$timepoint
  sample_metadata$timepoint2[!(sample_metadata$timepoint2 %in% c("Pre", "IPE"))] <- NA
  sample_metadata$timepoint2 <- factor(sample_metadata$timepoint2, levels = c("Pre", "IPE"))
  
  sample_metadata$study <- "ST001789"
  sample_metadata$platform <- filename
  sample_metadata$species <- "human"
  sample_metadata$dataset <- paste0(unique(sample_metadata$study),
                                    ":", filename)
  sample_metadata$exhaustion <- NA
  sample_metadata$has_metadata <- 1
  rownames(sample_metadata) <- sample_metadata$sample
  
  feature_metadat <- data.frame(original = metab_names,
                                clean_metabolites = make_clean_names(metab_names),
                                refmet_name = refmet_names)
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
  
  write.csv(dat, paste0(parent_folder, output_folder, filename, "_ST001789_fitered_imputed_log_dat.csv"))
  write.csv(feature_metadat, paste0(parent_folder, output_folder, filename, "_ST001789_feature_metadata.csv"))
  
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
save(processed_dat, file = paste0(parent_folder, table_output_folder, "04-ST001789.Rdata"))








