library(readxl)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(viridis)
filefolder <- "../"
input_folder <- "03-results/"
output_folder <- "04-graphs/"
tryCatch({
  dir.create(paste0(filefolder, output_folder))
  message(paste0("Creating ",paste0(filefolder, output_folder," folder!")))
},
warning = function(w){
  message("Could not create Folder!")
  message(conditionMessage(w))})

dat1 <- read_xlsx(paste0(filefolder, input_folder, "PREvsIPE_metaanalysis_03132025.xlsx"), 
                  sheet = "summary")
dat1 <- data.frame(dat1)
rownames(dat1) <- dat1$clean_metabolites

# top_metabs <- dat1[dat1$rat_qval <= 1e-2 & dat1$Human_qval <= 1e-2,]
top_metabs <- dat1[dat1$Human_qval <= 0.05 & abs(as.numeric(dat1$Human_coef)) >= 0.1,]
top_metabs <- top_metabs[complete.cases(top_metabs$clean_metabolites),]
top_metabs <- top_metabs[, c(seq(10), seq(15,18))]
top_metabs <- top_metabs[order(-abs(top_metabs$Human_coef)),]
top_metabs <- top_metabs[seq(100),]
top_metabs <- top_metabs[order(-top_metabs$Human_coef),]

datR <- data.frame(rat = dat1[top_metabs$clean_metabolites, "rat_coef"],
                   difference = dat1[top_metabs$clean_metabolites, "interaction_coef"],
                   exhaustion_max = NA,
                   row.names = top_metabs$clean_metabolites)
datR$rat <- as.numeric(datR$rat)
pheatmap(datR, color = inferno(500),
         breaks = seq(-1.5, 1.5, length.out = 501),
         cellwidth = 20, cellheight = 10, fontsize = 11,
         labels_row = top_metabs$refmet_name,
         labels_col = c("rodent", "difference", "exhaustion-max"),
         cluster_rows = FALSE, cluster_cols = FALSE,
         legend = TRUE, legend_breaks = c(-1.25, -0.75, 0, 0.75, 1.25),
         show_rownames = TRUE, show_colnames = TRUE,
         filename = paste0(filefolder, output_folder, "rodent_heatmap2.jpeg"))


datH <- read_xlsx(paste0(filefolder, input_folder, "PREvsIPE_metaanalysis_03132025.xlsx"), 
                  sheet = "human")
datH <- data.frame(datH)
rownames(datH) <- datH$clean_metabolites
datH <- datH[top_metabs$clean_metabolites,]
datH <- data.frame(human = datH$timepoint2IPE_coef,
                   species = datH$speciesrat_coef,
                   exhaustion= datH$exhaustionmax_coef,
                   age = datH$age_coef,
                   sex_male = datH$sex1_coef,
                   BMI = datH$BMI_coef,
                   row.names = top_metabs$clean_metabolites)

datH <- datH[, -2]
pheatmap(datH, color = inferno(500),
         breaks = seq(-1.5, 1.5, length.out = 501),
         cellwidth = 20, cellheight = 10, fontsize = 11,
         labels_row = top_metabs$refmet_name,
         labels_col = c("human", "species", "exhaustion-max",
                        "age", "sex-male", "BMI"),
         cluster_rows = FALSE, cluster_cols = FALSE,
         legend = TRUE, legend_breaks = c(-1.25, -0.75, 0, 0.75, 1.25),
         show_rownames = TRUE, show_colnames = TRUE,
         filename = paste0("~/Documents/maom_lab/F31/heatmap.jpeg"))


