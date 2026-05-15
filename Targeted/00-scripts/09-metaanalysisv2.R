library(stringr)
library(lme4)
library(lmerTest)
library(car)
library(writexl)
library(ggplot2)
parent_folder <- "../"
data_folder <- "02-processed_data/00-tables/"
file_names <- list.files(paste0(parent_folder, data_folder))
file_names <- file_names[!grepl("ST001789", file_names)]
names(file_names) <- vapply(file_names, FUN = function(x) str_replace(x, ".Rdata", ""),
                            FUN.VALUE = character(1))


metab_table <- data.frame(matrix(data=NA, nrow=0, ncol=4))
colnames(metab_table) <- c("refmet", "count", "study", "platform")
refmet_out <- vector("list", length(file_names))
names(refmet_out) <- names(file_names)
feature_names <- data.frame(matrix(data=NA, nrow=0, ncol=2))
for(j in names(file_names)){
  load(paste0(parent_folder, data_folder, file_names[j]))
  tmp_out <- vector("list", length(names(processed_dat)))
  names(tmp_out) <- names(processed_dat)
  for(i in names(processed_dat)){
    feature_names <- rbind(feature_names, processed_dat[[i]][["features"]][, c("original", 
                                                                               "clean_metabolites",
                                                                               "refmet_name")])
    tmp <- processed_dat[[i]][["refmet_dat"]]
    tmp <- tmp[tmp$timepoint2 %in% c("Pre", "IPE"), ]
    tmp$sex[tmp$sex == "Male"] <- "M"
    tmp$sex[tmp$sex == "Female"] <- "F"
    tmp_out[[i]] <- tmp
    metab1 <- match("has_metadata", colnames(tmp))+1
    for(y in seq(metab1, ncol(tmp))){
      tmp_metab <- colnames(tmp)[y]
      ind <- match(tmp_metab, metab_table$refmet)
      if(is.na(ind)){
        metab_table <- rbind(metab_table,
                             data.frame(refmet=tmp_metab, 
                                        count = 1,
                                        study = j, 
                                        platform = paste0(j, ":", i)))
      }else{
        metab_table$count[ind] = metab_table$count[ind] + 1
        metab_table$study[ind] <- paste(metab_table$study[ind], j, sep = " , ")
        metab_table$platform[ind] <- paste(metab_table$platform[ind], paste0(j, ":", i), sep = " , ")
      }
    }
  }
  refmet_out[[j]] <- tmp_out
}
rm(tmp_metab, ind)

feature_names <- feature_names[match(unique(feature_names$clean_metabolites), 
                                     feature_names$clean_metabolites),]
for(i in colnames(feature_names)){
  feature_names[[i]] <- unlist(feature_names[[i]])
}

metab_table2 <- metab_table[metab_table$count >=2, ]
metab_table2$study <- vapply(metab_table2$study, 
                             FUN = function(x) paste(unique(unlist(strsplit(x, ","))), collapse = ","),
                             FUN.VALUE = character(1))
metab_table2$study2 <- lapply(seq(nrow(metab_table2)), 
                              function(x) unique(unlist(strsplit(metab_table2$study[x], ","))))
metab_table2$human <- F
metab_table2$rat <- F

for(i in seq(nrow(metab_table2))){
  if(any(grepl("02-Motrpac_human", metab_table2$study[i]),
         grepl("03-FUEL", metab_table2$study[i]),
         grepl("04-ST001789", metab_table2$study[i]),
         grepl("07-ST001749", metab_table2$study[i]),
         grepl("05-Contrepois", metab_table2$study[i]))){
    metab_table2$human[i] <- T
  }
  if(any(grepl("01-Motrpac_rat", metab_table2$study[i]),
         grepl("06-Sato", metab_table2$study[i]),
         grepl("08-LCR_HCR", metab_table2$study[i]))){
    metab_table2$rat[i] <- T
  }
  
}

dvals <- lapply(metab_table2$platform, function(x) unlist(strsplit(x, " , ")))
dvals <- unique(unlist(dvals))
dvals <- dvals[order(dvals)]
dvals <- unlist(lapply(dvals, function(x) str_replace(x, "0[0-9]-", "")))
names(dvals) <- dvals
dvals <- unlist(lapply(dvals, function(x) str_replace(x, ":", "_")))
meta_cols <- c("subject", "sample", "age", "sex", "weight", "BMI", 
               "study", "platform", "species", "timepoint2", "dataset",
               "exhaustion", "has_metadata")
res_cols <- c("A-timepoint2IPE", "B-speciesrat",
              "C-speciesrat:timepoint2IPE",
              "E-age","F-sex1", "G-BMI", "H-exhaustionmax")
#x="01-Motrpac_rat";y="rppos";i="alanine"
foutput <- NULL
routput <- NULL
fail_count <- c()
simple_count <- c()
all_warnings <- vector("list", length(metab_table2$refmet))
names(all_warnings) <- metab_table2$refmet
for(i in metab_table2$refmet){
  print(paste0("metabolite: ", i))
  #gather data
  tmp <- NULL
  for(x in names(refmet_out)){
    for(y in names(refmet_out[[x]])){
      if(i %in% colnames(refmet_out[[x]][[y]])){
        tmp2 <- refmet_out[[x]][[y]][, c(meta_cols, i)]
        tmp <- rbind(tmp, tmp2)
        rm(tmp2)
      }
    }
  }
  
  tmp$BMI[is.na(tmp$BMI)] <- 0
  tmp$age[is.na(tmp$age)] <- 0
  tmp$sex[is.na(tmp$sex)] <- 0
  
  #tmp <- tmp[tmp$study != "Sato", ]
  #tmp$subject <- factor(tmp$subject)
  tmp$species <- factor(tmp$species)
  tmp$sex <- ifelse(tmp$sex  == "M", 1, 0)
  tmp$sex <- factor(tmp$sex, levels = c(0, 1))
  tmp$sex <- factor(tmp$sex)
  tmp$sex[tmp$species == "rat"] <- 0
  tmp$timepoint2 <- factor(tmp$timepoint2, levels = c("Pre", "IPE"))
  #tmp$dataset <- factor(tmp$dataset)
  tmp$exhaustion <- factor(tmp$exhaustion, levels = c("submax", "max"))
  tmp$has_metadata <- factor(tmp$has_metadata, levels = c(0,1),
                             labels = c("No", "Yes"))
  
  #create model equation
  onlyFuel <- length(unique(tmp$study[tmp$species == "human"])) == 1 &
              unique(tmp$study[tmp$species == "human"])[1] == "FUEL"
  onlyST001749 <- length(unique(tmp$study[tmp$species == "human"])) == 1 &
                  unique(tmp$study[tmp$species == "human"])[1] == "ST001749"
  onlyFuelSt001749 <- length(unique(tmp$study[tmp$species == "human"])) == 2 &
                      all(c("FUEL", "ST001749") %in% unique(tmp$study[tmp$species == "human"]))
  if(all(c("human", "rat") %in% unique(tmp$species))){
    check = "human_rat"
  }else if("human" %in% unique(tmp$species)){
    check = "human"
  }else if("rat" %in% unique(tmp$species)){
    check = "rat"
  }else{
    stop("Error!")
  }
  if(check == "human_rat"){
    fit_eq <- paste0(i, " ~ species + timepoint2 + exhaustion + ",
                     "age*I(has_metadata == 'Yes') + sex*I(has_metadata == 'Yes') ",
                     "+ BMI*I(has_metadata == 'Yes') + species:timepoint2 + ",
                     "(1 | dataset/subject)")
    fit_eq2 <- paste0(i, " ~ species + timepoint2 + exhaustion + ",
                      "species:timepoint2 + (1 | dataset/subject)")
    if(onlyFuelSt001749 | onlyFuel){
      fit_eq <- str_replace(fit_eq, "\\+ sex\\*I\\(has_metadata == 'Yes'\\) ", "")
    }else if(onlyST001749){
      fit_eq <- str_replace(fit_eq, "\\+ age\\*I\\(has_metadata == 'Yes'\\) ", "")
      fit_eq <- str_replace(fit_eq, "\\+ sex\\*I\\(has_metadata == 'Yes'\\) ", "")
      fit_eq <- str_replace(fit_eq, "\\+ BMI\\*I\\(has_metadata == 'Yes'\\) ", "")
      fit_eq <- str_replace(fit_eq, "\\(1 \\| dataset/subject\\)", "(1 | dataset)")
      fit_eq2 <- str_replace(fit_eq2, "\\(1 \\| dataset/subject\\)", "(1 | dataset)")
    }
  }else if(check == "human"){
    fit_eq <- paste0(i, " ~ timepoint2 + exhaustion + ",
                     "age*I(has_metadata == 'Yes') + sex*I(has_metadata == 'Yes') ",
                     "+ BMI*I(has_metadata == 'Yes') + ",
                     "(1 | dataset/subject)")
    fit_eq <- paste0(i, " ~ timepoint2 + ",
                     "age*I(has_metadata == 'Yes') + sex*I(has_metadata == 'Yes') ",
                     "+ BMI*I(has_metadata == 'Yes') + ",
                     "(1 | dataset/subject)")
    fit_eq2 <- paste0(i, " ~ timepoint2 + exhaustion + ",
                      "(1 | dataset/subject)")
    if(onlyFuelSt001749){
      fit_eq <- str_replace(fit_eq, "\\+ sex\\*I\\(has_metadata == 'Yes'\\) ", "")
    }
  }else if(check == "rat"){
    fit_eq <- paste0(i, " ~ timepoint2 + exhaustion + (1 | dataset)")
    fit_eq2 <- paste0(i, " ~ timepoint2 + exhaustion + (1 | dataset)")
  }else{
    stop("There was a problem with the formula!")
  }
  #remove exhaustion if all are half-max or max
  if(length(unique(tmp$exhaustion)) == 1){
    fit_eq <- str_replace(fit_eq, "exhaustion \\+ ", "")
    fit_eq2 <- str_replace(fit_eq2, "exhaustion \\+ ", "")
  }
  
  fit_eqf <- as.formula(fit_eq)
  fit_eqf2 <- as.formula(fit_eq2)
  # fit <- lmer(fit_eq, data = tmp)
  warns <- list()
  out_func <- fit_eq
  fit <- withCallingHandlers(lmer(fit_eqf, data = tmp), 
                             warning = function(warn){warns <<- append(warns, list(warn))})
  warns <- lapply(warns, function(x) x[["message"]])
  tmodel = "full"
  if(any(grepl('Model failed to converge', warns))){
    print(paste0("metabolite: ", i, " failed to converge. ",
                 "Trying simple model."))
    warns <- list()
    out_func <- fit_eq2
    fit <- withCallingHandlers(lmer(fit_eqf2, data = tmp), 
                               warning = function(warn){warns <<- append(warns, list(warn))})
    warns <- lapply(warns, function(x) x[["message"]])
    tmodel = "simple"
    if(any(grepl('Model failed to converge:', warns))){
      print(paste0("Simple model failed to converge!"))
      fail_count <- c(fail_count, i)
      fit <- NA
      tmodel = NA
    }else{
      simple_count <- c(simple_count, i)
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
  
  col_values <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  names(col_values) <- c("_coef", "_sd", "_tvalue", "_pvalue")
  restmp <- list(metabolite = i, 
                 total_samples = nrow(tmp),
                 total_studies = length(unique(tmp$dataset)),
                 studies = paste0(unique(tmp$dataset), collapse = ", "),
                 eq = out_func,
                 converged = conv,
                 singular = sing,
                 model = tmodel)
  for(c in res_cols){
    for(c2 in names(col_values)){
      restmp[[paste0(c, c2)]] <- tryCatch({res[str_replace(c, "[A-Z]-", ""), col_values[[c2]]]},
                                          error = function(e){
                                            NA
                                          })
    }
  }
  rm(c, c2)
  if(!is.na(tmodel) & grepl("species:timepoint2", out_func)){
    lh = linearHypothesis(fit, "timepoint2IPE + speciesrat:timepoint2IPE = 0")
    restmp[["D_LH_chi2"]] <- lh[2, 2]
    restmp[["D_LH_pval"]] <- lh[2, 3]
  }
  foutput <- rbind(foutput, unlist(restmp))
  
  raneftmp <- vector("list", length = length(dvals)+1)
  names(raneftmp) <- c("metabolite", dvals)
  raneftmp["metabolite"] <- i
  for(c in seq(length(dvals))){
    raneftmp[dvals[c]] <- tryCatch({ranef(fit)$dataset[names(dvals[c]), ]},
                                   error = function(e){
                                     NA
                                   })
  }
  rm(c)
  routput <- rbind(routput, unlist(raneftmp))
 
}


routput <- data.frame(routput)
foutput <- data.frame(foutput)
foutput2 <- foutput
for(i in colnames(foutput2)){
  foutput2[[i]] <- unlist(foutput2[[i]])
}
res_cols <- str_replace_all(res_cols, ":", ".")
for(i in res_cols){
  foutput2[[paste0(i, "_qvalue")]] <- p.adjust(foutput2[[paste0(i, "_pvalue")]], method = "BH")
}
foutput2 <- merge(feature_names, foutput2, by.x = "clean_metabolites", by.y = "metabolite",
                  all.x = FALSE, all.y = TRUE)
foutput2 <- foutput2[order(foutput2$timepoint2IPE_pvalue), ]
for(i in colnames(foutput2)){
  foutput2[[i]] <- as.character(foutput2[[i]])
}
keep <- vapply(seq(nrow(foutput2)), function(x) foutput2$converged[x] == "TRUE",
               logical(1))
foutput_reduce <- foutput2[keep,]
keepSD <- !grepl("_sd", colnames(foutput_reduce))
keepT <- !grepl("_tvalue", colnames(foutput_reduce))

foutput_reduce <- foutput_reduce[, keepSD & keepT]


write.csv(foutput2, paste0(parent_folder, "PREvsIPE_meta_03052025.csv"), row.names = FALSE)
write.csv(foutput_reduce, paste0(parent_folder, "PREvsIPE_meta_reduced_03052025.csv"), row.names = FALSE)
write.csv(routput, paste0(parent_folder, "PREvsIPE_ranef_03052025.csv"), row.names = FALSE)
write_xlsx(list(valid = foutput_reduce, all = foutput2, random_intercepts = routput),
           path = paste0(parent_folder, "PREvsIPE_metaanalysis_03052025.xlsx"))
foutput_reduce <- foutput_reduce[order(foutput_reduce$timepoint2IPE_pvalue),]

