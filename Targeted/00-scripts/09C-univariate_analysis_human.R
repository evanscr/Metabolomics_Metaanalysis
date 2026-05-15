library(stringr)
library(writexl)
library(janitor)
library(lmerTest)
parent_folder <- "../"
data_folder <- "02-processed_data/00-tables/"
output_folder <- "03-results/"
# file_names <- list.files(paste0(parent_folder, data_folder))
# file_names <- file_names[!grepl("ST001789", file_names)]
file_names <- list()
file_names[["02-Motrpac_human_endurance"]] <- "02-Motrpac_human_endurance.Rdata"
file_names[["02-Motrpac_human_resistance"]] <- "02-Motrpac_human_resistance.Rdata"
file_names[["03-FUEL"]] <- "03-FUEL.Rdata"
file_names[["05-Contrepois"]] <- "05-Contrepois.Rdata"
#file_names[["07-ST001749"]] <- "07-ST001749.Rdata"

names(file_names) <- vapply(file_names, FUN = function(x) str_replace(x, ".Rdata", ""),
                            FUN.VALUE = character(1))

model_forms = vector("list", length = length(file_names))
names(model_forms) <- names(file_names)
model_forms[["02-Motrpac_human_endurance"]] <- "~ age + sex + BMI + timepoint + (1 | subject)"
model_forms[["02-Motrpac_human_resistance"]] <- "~ age + sex + BMI + timepoint + (1 | subject)"
model_forms[["03-FUEL"]] <- "~ age + BMI + timepoint + (1 | subject)"
model_forms[["05-Contrepois"]] <- "~ age + sex + BMI + timepoint + (1 | subject)"
# model_forms[["07-ST001749"]] <- "~ timepoint + (1 | subject)"

concat_results <- vector("list", length(file_names))
names(concat_results) <- names(file_names)
for(filename in names(file_names)){
  load(paste0(parent_folder, data_folder, file_names[filename]))
  tryCatch({
    dir.create(paste0(parent_folder, output_folder, filename))
    message(paste0("Creating ",paste0(parent_folder, output_folder, filename,"/ folder!")))
  },
  warning = function(w){
    message("Could not create Folder!")
    message(conditionMessage(w))})
  
  study_results <- data.frame(matrix(data = NA, nrow = 0, ncol = 9, 
                                     dimnames = list(c(), c("clean_metabolites", "variable",
                                                            "Estimate", "Std. Error", 
                                                            "t value", "Pr(>|t|)", "qvalue",
                                                            "platform", "study"))))
  for(platform in names(processed_dat)){
    print(paste0("processing: ", filename, " - ", "platform: ", platform))
    dat <- processed_dat[[platform]][["data"]]
    sampl_metadat <- processed_dat[[platform]][["samples"]]
    sampl_metadat <- sampl_metadat[colnames(dat),]
    feature_metadat <- processed_dat[[platform]][["features"]]
    feature_metadat$clean_metabolites <- make_clean_names(feature_metadat$original)
    rownames(feature_metadat) <- feature_metadat$clean_metabolites
    feat_cols <- c("original", "clean_metabolites", "refmet_name", 
                   "rt", "mz", "mass", "formula")
    feat_cols <-  feat_cols[feat_cols %in% colnames(feature_metadat)]
    feature_metadat <- feature_metadat[, feat_cols]
    
    metabs <- rownames(dat)[rownames(dat) %in% rownames(feature_metadat)]
    feature_metadat <- feature_metadat[metabs,]
    dat <- dat[metabs,]
    if(!all(rownames(feature_metadat) == rownames(dat))) stop("feature metadat does not match data order!")
    if(!all(rownames(sampl_metadat) == colnames(dat))) stop("Sample metadat does not match data order!")
    output <- data.frame(matrix(data = NA, nrow = 0, ncol = length(feat_cols)+8, 
                                dimnames = list(c(), c(feat_cols, "converged", "singular",
                                                       "variable", "Estimate", "Std. Error", 
                                                       "t value", "Pr(>|t|)", "qvalue"))))
    fail_count <- c()
    simple_count <- c()
    for(m in rownames(feature_metadat)){
      tmp <- sampl_metadat
      tmp[[m]] <- unlist(dat[m,])
      tmp_form <- paste0(m, model_forms[[filename]])
      
      warns <- list()
      fit <- withCallingHandlers(lmer(as.formula(tmp_form), data = tmp), 
                                 warning = function(warn){warns <<- append(warns, list(warn))})
      warns <- lapply(warns, function(x) x[["message"]])
      tmodel = "full"
      if(any(grepl('Model failed to converge', warns))){
        print(paste0("metabolite: ", m, " failed to converge. ",
                     "Trying simple model."))
        
        tmp_form <- paste0(m, " ~ timepoint + (1 | subject)")
        warns <- list()
        fit <- withCallingHandlers(lmer(as.formula(tmp_form), data = tmp), 
                                   warning = function(warn){warns <<- append(warns, list(warn))})
        warns <- lapply(warns, function(x) x[["message"]])
        tmodel = "simple"
        if(any(grepl('Model failed to converge:', warns))){
          print(paste0("Simple model failed to converge!"))
          fail_count <- c(fail_count, m)
          fit <- NA
          tmodel = NA
        }else{
          simple_count <- c(simple_count, m)
        }
      }
      
      if(!inherits(fit, "lmerModLmerTest")){
        conv <- FALSE
      }else{
        res <- summary(fit)$coefficients
        conv <- TRUE
        if(isSingular(fit)){
          sing = TRUE
        }else{
          sing = FALSE
        }
      }
      
      tryCatch({
        tmp_feat <- feature_metadat[match(m, feature_metadat$clean_metabolites), ]
        res <- summary(fit)$coefficients[, -match("df", colnames(summary(fit)$coefficients))]
        res <- cbind(tmp_feat[rep(1, nrow(res)),], 
                     data.frame(converged=conv, singular=sing, 
                                variable = rownames(res), res))
        output <- rbind(output, res)
      }, error = function(e){
        message(paste0("Could not estimate: ", m, "\n"))
        message(conditionMessage(e))
      })
      
    }
    output$qvalue <- NA
    colnames(output) <- c(feat_cols, "converged", "singular", 
                          "variable", "estimate", "st_error",
                          "t_value", "pvalue", "qvalue")
    output <- output[output$variable != "(Intercept)",]
    
    for(var in unique(output$variable)){
      tmp <- output$pvalue[output$variable == var]
      output$qvalue[output$variable == var] <- p.adjust(tmp, method = "BH")
    }

    output$platform = platform
    output$study = filename
    write.csv(output, paste0(parent_folder, output_folder, filename, "/",
                             filename, "_", platform, "_", "study_results.csv"), row.names = FALSE)
    
    study_results <- rbind(study_results, output)
  }
  concat_results[[filename]] <- study_results
}
write_xlsx(concat_results, path = paste0(parent_folder, output_folder, "human_DE_by_study.xlsx"))
keep_cols <- c("original", "clean_metabolites", "refmet_name", "variable", 
               "estimate", "st_error", "t_value", "pvalue", "qvalue", 
               "platform", "study")
comb <- list()
for(study in names(concat_results)){
  tmp = concat_results[[study]]
  tmp <- tmp[, keep_cols]
  comb[[study]] <- tmp
}
all_res <- Reduce("rbind", comb)

write.csv(all_res, paste0(parent_folder, output_folder, "human_unistudy_DE_table.csv"), row.names = FALSE)
# for(study in concat_results){
#   study <- study[study$variable != "(Intercept)",]
#   study <- study[order(study)]
# }