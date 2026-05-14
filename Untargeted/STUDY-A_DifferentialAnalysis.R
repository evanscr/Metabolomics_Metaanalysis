rm(list=ls())
setwd('C:/Users/Abrraski/OneDrive - Michigan Medicine/Desktop/DifferentialExpressionAnalysis')
getwd()
#* Now that we can install packages, we have to load them into our working memory each time we open a new
#* R ression

library(DESeq2)
library(tidyverse)
###############
#read in the data set
Pass1AZ <- read.csv('MoTrPACPass1A-ZDrive-NamedAndUnnamed.csv')
#Extract Metadata and metabolites, and get rid of empty columns
Pass1AZMetabolites <- Pass1AZ[,1:5]
Pass1AZMetabolites <- Pass1AZMetabolites[-c(1),-c(2,3)]
Pass1AZMeta <- Pass1AZ[1,]
Pass1AZMeta <- Pass1AZMeta[,-c(1:5)]
#remove metadata from dataframe and place samples in columns and features in rows
Pass1AZ <- Pass1AZ[-c(1),-c(1:5)]
rownames(Pass1AZ)<- Pass1AZMetabolites[,1]
#Check class of frame to ensure numeric
class(Pass1AZ[,5])
#This comes back as character, so I am going to convert it to numeric and filter for missingness
Pass1AZ <- as.data.frame(apply(Pass1AZ,2,as.numeric))
rownames(Pass1AZ)<-Pass1AZMetabolites[,1]
class(Pass1AZ[,5])
any(is.na(Pass1AZ))
#for some reason It was filtering out way too many features with NA's that shouldn't have been, so I replaced the NA with
#Zeroes, then filtered
Pass1AZ <- Pass1AZ %>% replace(is.na(.), 0)
Pass1AZFiltered <- Pass1AZ %>% filter ((rowSums((. > 0)) / ncol(.)) > 0.3, na.rm = TRUE)
#Now I'm going to swtich them back and replace the NA values with 1/5 the min of the row
Pass1AZFiltered[Pass1AZFiltered==0] <- NA
for (i in 1:nrow(Pass1AZFiltered)){
  Pass1AZFiltered[i,][is.na(Pass1AZFiltered[i,])] <- 0.2 * min(Pass1AZFiltered[i,][!(is.na(Pass1AZFiltered[i,]))])
}
#Now I'm going to log transform the data
Pass1AZFiltered <- as.data.frame(apply(Pass1AZFiltered,2, log2))
#We can scale the data and look at some histograms real quick
Pass1AZFiltered <- as.data.frame(t(Pass1AZFiltered))
Pass1AZMeta <- as.data.frame(t(Pass1AZMeta))
Pass1AZFilteredScaled <- as.data.frame(apply(Pass1AZFiltered,1, function(x) (x-mean(x))/sd(x)))
Pass1AZFilteredScaled <- as.data.frame(t(Pass1AZFilteredScaled))
#Histograms
ggplot(Pass1AZFiltered, aes(x=`CAR 16:2`))+geom_histogram()

ggplot(Pass1AZFilteredScaled, aes(x=`CAR 16:2`))+geom_histogram()
#Add metadata and separate by condition
Pass1AZFiltered$Treatment <- Pass1AZMeta$`1`
Pass1AZFilteredScaled$Treatment <- Pass1AZMeta$`1`
#Separate the data first by exercise type, then by sample collection timepoint
Pass1AZBaseline <- filter(Pass1AZFiltered, Treatment == 'Baseline')
Pass1AZIPE <- filter(Pass1AZFiltered, Treatment == 'IPE')
Pass1AZPost30m <- filter(Pass1AZFiltered, Treatment == 'Post-30min')
Pass1AZPost1h <- filter(Pass1AZFiltered, Treatment == 'Post-1h')
Pass1AZPost4h <- filter(Pass1AZFiltered, Treatment == 'Post-4h')
Pass1AZPost7h <- filter(Pass1AZFiltered, Treatment == 'Post-7h')
Pass1AZPost24h <- filter(Pass1AZFiltered, Treatment == 'Post-24h')
Pass1AZPost48h <- filter(Pass1AZFiltered, Treatment == 'Post_48h')


#Scrub MetaData and perform t test
Pass1AZBaseline <- Pass1AZBaseline %>% select(!(Treatment))
Pass1AZIPE <- Pass1AZIPE %>% select(!(Treatment))
Pass1AZPost30m <- Pass1AZPost30m %>% select(!(Treatment))
Pass1AZPost1h <- Pass1AZPost1h %>% select(!(Treatment))
Pass1AZPost4h <- Pass1AZPost4h %>% select(!(Treatment))
Pass1AZPost7h <- Pass1AZPost7h %>% select(!(Treatment))
Pass1AZPost24h <- Pass1AZPost24h %>% select(!(Treatment))
Pass1AZPost48h <- Pass1AZPost48h %>% select(!(Treatment))

varTest <- data.frame(rep(NA,ncol(Pass1AZBaseline)))
for (i in 1:ncol(Pass1AZBaseline)){
  varTest[i,1] <- var(Pass1AZBaseline[,i])
}
rownames(varTest24h) <- colnames(Pass1AZPost24h)

ToBeKept <- varTest %>% filter((.>0))
ToBeKept48h <-
ToBeKept48h <- varTest48h %>% filter(.>0)

reorder48h <- match(rownames(ToBeKept48h),colnames(Pass1AZBaseline))
Baseline48hFiltered <- Pass1AZBaseline[,reorder48h]
Filtered48h <- Pass1AZPost48h[,reorder48h]
reorder24h <- match(rownames(ToBeKept24h),colnames(Pass1AZBaseline))
Baseline24hFiltered <- Pass1AZBaseline[,reorder24h]
Filtered24h <- Pass1AZPost24h[,reorder24h]
#Now we can run the t-test
Pass1AZOutputPost24h <- data.frame(metabolites = colnames(Pass1AZBaselineTest),
                           t_statistic = rep(NA,length(colnames(Pass1AZBaselineTest))),
                           p_value = rep(NA,length(colnames(Pass1AZBaselineTest))))
# we need to remove the groups column since we cant run a t.test on that

for (i in 1:ncol(Pass1AZBaselineTest)){
  
  t_test<- t.test(Pass1AZBaselineTest[,i], Pass1AZPost24hTest[,i], var.equal = FALSE)
  Pass1AZOutputPost24hTest$t_statistic[i] <- t_test$statistic
  Pass1AZOutputPost24hTest$p_value[i] <- t_test$p.value
}

#lets adjust the pvalues

Pass1AZOutputPost48h$adj_pvalues <- p.adjust(Pass1AZOutputPost48h$p_value, method = 'BH')
Pass1AZOutputPost48h <- as.data.frame(Pass1AZOutputPost48h)

BaselineMeans <- as.data.frame(t(colMeans(Pass1AZBaseline)))
IPEMeans <- as.data.frame(t(colMeans(Pass1AZIPE)))
Post30mMeans <- as.data.frame(t(colMeans(Pass1AZPost30m)))
Post1hMeans <- as.data.frame(t(colMeans(Pass1AZPost1h)))
Post4hMeans <- as.data.frame(t(colMeans(Pass1AZPost4h)))
Post7hMeans <- as.data.frame(t(colMeans(Pass1AZPost7h)))
Post24hMeans <- as.data.frame(t(colMeans(Pass1AZPost24h)))
Post48hMeans <- as.data.frame(t(colMeans(Pass1AZPost48h)))

rownames(BaselineMeans) <- 'Mean'
rownames(IPEMeans) <- 'Mean'
rownames(Post30mMeans) <- 'Mean'
rownames(Post1hMeans) <- 'Mean'
rownames(Post4hMeans) <- 'Mean'
rownames(Post7hMeans) <- 'Mean'
rownames(Post24hMeans) <- 'Mean'
rownames(Post48hMeans) <- 'Mean'


Pass1AZBaseline <- union(Pass1AZBaseline,BaselineMeans)
Pass1AZPost4h <- union(Pass1AZPost4h, Post4hMeans)
Pass1AZPost24h <- union(Pass1AZPost24h, Post24hMeans)
Pass1AZPost30m <- union(Pass1AZPost30m, Post30mMeans)
Pass1AZPost48h <- union(Pass1AZPost48h, Post48hMeans)
Pass1AZPost1h <- union(Pass1AZPost1h, Post1hMeans)
Pass1AZPost7h <- union(Pass1AZPost7h, Post7hMeans)
Pass1AZIPE <- union(Pass1AZIPE, IPEMeans)

Pass1AZBaseline <- as.data.frame(t(Pass1AZBaseline))
Pass1AZPost30m <- as.data.frame(t(Pass1AZPost30m))
Pass1AZPost4h <- as.data.frame(t(Pass1AZPost4h))
Pass1AZPost24h <- as.data.frame(t(Pass1AZPost24h))
Pass1AZPost1h <- as.data.frame(t(Pass1AZPost1h))
Pass1AZPost7h <- as.data.frame(t(Pass1AZPost7h))
Pass1AZPost48h <- as.data.frame(t(Pass1AZPost48h))
Pass1AZIPE <- as.data.frame(t(Pass1AZIPE))

Pass1AZOutputIPE$log2_fold_change <- Pass1AZIPE$Mean - Pass1AZBaseline$Mean
Pass1AZOutputPost30m$log2_fold_change <- Pass1AZPost30m$Mean - Pass1AZBaseline$Mean
Pass1AZOutputPost1h$log2_fold_change <- Pass1AZPost1h$Mean - Pass1AZBaseline$Mean
Pass1AZOutputPost4h$log2_fold_change <- Pass1AZPost4h$Mean - Pass1AZBaseline$Mean
Pass1AZOutputPost7h$log2_fold_change <- Pass1AZPost7h$Mean - Pass1AZBaseline$Mean
Pass1AZOutputPost24h$log2_fold_change <- Pass1AZPost24h$Mean - Pass1AZBaseline$Mean
Pass1AZOutputPost48h$log2_fold_change <- Pass1AZPost48h$Mean - Pass1AZBaseline$Mean



write.csv(Pass1AZOutputPost48h, "C:/Users/Abrraski/OneDrive - Michigan Medicine/Desktop/DifferentialExpressionAnalysis/Pass1A_Untargeted_Post48h_DifferentrialAnalysis.csv",row.names =F)

