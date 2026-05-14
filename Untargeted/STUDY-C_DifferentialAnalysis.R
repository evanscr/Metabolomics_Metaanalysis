rm(list=ls())
#* Now that we can install packages, we have to load them into our working memory each time we open a new
#* R ression

library(DESeq2)
library(tidyverse)
setwd("C:/Users/Abrraski/OneDrive - Michigan Medicine/Desktop/DifferentialExpressionAnalysis")
getwd()
Contrepois <- read.csv('Contrepois_normalized_by_SERRF_with_NAs_WithFactors.csv')
ContrepoisMeta <- Contrepois[1:2,]
rownames(Contrepois)<- Contrepois[,1]
Contrepois <- Contrepois[-c(1:2),-c(1)]
ContrepoisMetabolites <- as.data.frame(rownames(Contrepois))
rownames(ContrepoisMeta)<- ContrepoisMeta[,1]
ContrepoisMeta<- ContrepoisMeta[,-c(1)]

#Filter out any metabolites with greater than 30% missingness
any(is.na(Contrepois))
Contrepois <- Contrepois %>% replace(is.na(.),0)
ContrepoisFiltered <- Contrepois %>% filter ((rowSums((.>0))/ ncol(.)) > 0.3, na.rm = TRUE)
#Replace 0 with NA and set to 1/5 min intensity per metabolite

ContrepoisFiltered[ContrepoisFiltered ==0] <- NA
ContrepoisFiltered  <- as.data.frame(apply(ContrepoisFiltered,2,as.numeric))
rownames(ContrepoisFiltered) <- ContrepoisMetabolites[,1]
for (i in 1:nrow(ContrepoisFiltered)){
  ContrepoisFiltered[i,][is.na(ContrepoisFiltered[i,])] <- 0.2 * min(ContrepoisFiltered[i,][!(is.na(ContrepoisFiltered[i,]))])
}

#Now that the data is filtered for missingness and NA values, we can log transform it to normalize
ContrepoisFiltered  <- as.data.frame(apply(ContrepoisFiltered,2,log2))
ContrepoisFiltered <- as.data.frame(t(ContrepoisFiltered))
ContrepoisScaled <- data.frame(apply(ContrepoisFiltered,1, function(x) (x-mean(x))/sd(x)))
ContrepoisScaled <- as.data.frame(t(ContrepoisScaled))

#I want to look at a histogram of the data real quick just to see if it needs to be scaled at all

ggplot(ContrepoisFiltered, aes(x = `UNK_7.783_151.01593`))+ geom_histogram()
ggplot(ContrepoisScaled, aes(x=`UNK_7.783_151.01593`))+ geom_histogram()
#Im going to make a dataframe Just containing the two conditions I want to compare]
ContrepoisMeta <- as.data.frame(t(ContrepoisMeta))

#just gonna check the sample names are in the same order rq
ContrepoisFiltered$Treatment <- ContrepoisMeta$Treatment
ContrepoisFiltered$Id <- ContrepoisMeta$ID
#Ok now Im gonna make a dataframe just containing the two treatments Im comparing from
ContrepoisBaseline <- filter(ContrepoisFiltered, Treatment == 'Baseline')
ContrepoisIPE <- filter(ContrepoisFiltered, Treatment == '2 min')
ContrepoisPost15m <- filter(ContrepoisFiltered, Treatment == '15 min')
ContrepoisPost30m <- filter(ContrepoisFiltered, Treatment == '30 min')
ContrepoisPost1h <- filter(ContrepoisFiltered, Treatment == '1 hour')
ContrepoisPost2h <- filter(ContrepoisFiltered, Treatment == '2 hours')
ContrepoisPost4h <- filter(ContrepoisFiltered, Treatment == '4 hours')
ContrepoisPost6h <- filter(ContrepoisFiltered, Treatment == '6 hours')
ContrepoisPost24h <- filter(ContrepoisFiltered, Treatment == '24 hours')

#Reorder datasets for paired analysis
reorderIPE <- match(ContrepoisIPE$Id,ContrepoisBaseline$Id)
BaselinePairIPE <- ContrepoisBaseline[reorderIPE,]

reorderPost15m <- match(ContrepoisPost15m$Id,ContrepoisBaseline$Id)
BaselinePair15m <- ContrepoisBaseline[reorderPost15m,]

reorderPost30m <- match(ContrepoisPost30m$Id,ContrepoisBaseline$Id)
BaselinePair30m <- ContrepoisBaseline[reorderPost30m,]

reorderPost1h <- match(ContrepoisPost1h$Id,ContrepoisBaseline$Id)
BaselinePair1h <- ContrepoisBaseline[reorderPost1h,]

reorderPost2h <- match(ContrepoisPost2h$Id,ContrepoisBaseline$Id)
BaselinePair2h <- ContrepoisBaseline[reorderPost2h,]

reorderPost4h <- match(ContrepoisPost4h$Id,ContrepoisBaseline$Id)
BaselinePair4h <- ContrepoisBaseline[reorderPost4h,]

reorderPost6h <- match(ContrepoisPost6h$Id,ContrepoisBaseline$Id)
BaselinePair6h <- ContrepoisBaseline[reorderPost6h,]

reorderPost24h <- match(ContrepoisPost24h$Id,ContrepoisBaseline$Id)
BaselinePair24h <- ContrepoisBaseline[reorderPost24h,]

all(BaselinePair2h$Id == ContrepoisPost30m$Id)

#scrub metadata 
ContrepoisBaseline <- ContrepoisBaseline %>% select(!(Treatment:Id))
ContrepoisIPE <- ContrepoisIPE %>% select(!(Treatment:Id))
ContrepoisPost30m <- ContrepoisPost30m %>% select(!(Treatment:Id))
ContrepoisPost1h <- ContrepoisPost1h %>% select(!(Treatment:Id))
ContrepoisPost15m <- ContrepoisPost15m %>% select(!(Treatment:Id))
ContrepoisPost2h <- ContrepoisPost2h %>% select(!(Treatment:Id))
ContrepoisPost4h <- ContrepoisPost4h %>% select(!(Treatment:Id))
ContrepoisPost6h <- ContrepoisPost6h %>% select(!(Treatment:Id))
ContrepoisPost24h <- ContrepoisPost24h %>% select(!(Treatment:Id))


BaselinePairIPE <- BaselinePairIPE %>% select(!(Treatment:Id))
BaselinePair15m <- BaselinePair15m %>% select(!(Treatment:Id))
BaselinePair30m <- BaselinePair30m %>% select(!(Treatment:Id))
BaselinePair1h <- BaselinePair1h %>% select(!(Treatment:Id))
BaselinePair2h <- BaselinePair2h %>% select(!(Treatment:Id))
BaselinePair4h <- BaselinePair4h %>% select(!(Treatment:Id))
BaselinePair6h <- BaselinePair6h %>% select(!(Treatment:Id))
BaselinePair24h <- BaselinePair24h %>% select(!(Treatment:Id))





#lets run a T-test for all of the metabolites in this dataset

#first lets initialize an output dataframe
ContrepoisPairedOutputPost24h <- data.frame(metabolites = colnames(ContrepoisBaseline),
                           t_statistic = rep(NA,length(colnames(ContrepoisBaseline))),
                           p_value = rep(NA,length(colnames(ContrepoisBaseline))))
# we need to remove the groups column since we cant run a t.test on that

for (i in 1:ncol(ContrepoisBaseline)){
  
  t_test<- t.test(BaselinePair30m[,i], ContrepoisPost30m[,i], var.equal = FALSE, paired = TRUE)
  ContrepoisPairedOutputPost30m$t_statistic[i] <- t_test$statistic
  ContrepoisPairedOutputPost30m$p_value[i] <- t_test$p.value
}

#lets adjust the pvalues

ContrepoisPairedOutputIPE$adj_pvalues <- p.adjust(ContrepoisPairedOutputIPE$p_value, method = 'BH')

BaselineMeansIPE <- as.data.frame(t(colMeans(BaselinePairIPE)))
BaselineMeansPost15m <- as.data.frame(t(colMeans(BaselinePair15m)))
BaselineMeansPost30m <- as.data.frame(t(colMeans(BaselinePair30m)))
BaselineMeansPost1h <- as.data.frame(t(colMeans(BaselinePair1h)))
BaselineMeansPost2h <- as.data.frame(t(colMeans(BaselinePair2h)))
BaselineMeansPost4h <- as.data.frame(t(colMeans(BaselinePair4h)))
BaselineMeansPost6h <- as.data.frame(t(colMeans(BaselinePair6h)))
BaselineMeansPost24h <- as.data.frame(t(colMeans(BaselinePair24h)))

IPEMeans <- as.data.frame(t(colMeans(ContrepoisIPE)))
Post15mMeans <- as.data.frame(t(colMeans(ContrepoisPost15m)))
Post30mMeans <- as.data.frame(t(colMeans(ContrepoisPost30m)))
Post1hMeans <- as.data.frame(t(colMeans(ContrepoisPost1h)))
Post2hMeans <- as.data.frame(t(colMeans(ContrepoisPost2h)))
Post4hMeans <- as.data.frame(t(colMeans(ContrepoisPost4h)))
Post6hMeans <- as.data.frame(t(colMeans(ContrepoisPost6h)))
Post24hMeans <- as.data.frame(t(colMeans(ContrepoisPost24h)))

rownames(BaselineMeansIPE) <- "Mean"
rownames(BaselineMeansPost15m) <- "Mean"
rownames(BaselineMeansPost30m) <- "Mean"
rownames(BaselineMeansPost1h) <- "Mean"
rownames(BaselineMeansPost2h) <- "Mean"
rownames(BaselineMeansPost4h) <- "Mean"
rownames(BaselineMeansPost6h) <- "Mean"
rownames(BaselineMeansPost24h) <- "Mean"
rownames(IPEMeans) <-"Mean"
rownames(Post15mMeans) <- "Mean"
rownames(Post30mMeans) <- "Mean"
rownames(Post1hMeans) <- "Mean"
rownames(Post2hMeans) <- "Mean"
rownames(Post4hMeans) <- "Mean"
rownames(Post6hMeans) <- "Mean"
rownames(Post24hMeans) <- "Mean"

BaselinePair15m <- union(BaselinePair15m,BaselineMeansPost15m)
BaselinePair30m <- union(BaselinePair30m,BaselineMeansPost30m)
BaselinePair1h <- union(BaselinePair1h,BaselineMeansPost1h)
BaselinePair2h <- union(BaselinePair2h,BaselineMeansPost2h)
BaselinePair4h <- union(BaselinePair4h,BaselineMeansPost4h)
BaselinePair6h <- union(BaselinePair6h,BaselineMeansPost6h)
BaselinePair24h <- union(BaselinePair24h,BaselineMeansPost24h)
BaselinePairIPE <- union(BaselinePairIPE,BaselineMeansIPE)

ContrepoisIPE <- union(ContrepoisIPE, IPEMeans)
ContrepoisPost15m <- union(ContrepoisPost15m,Post15mMeans)
ContrepoisPost30m <- union(ContrepoisPost30m,Post30mMeans)
ContrepoisPost1h <- union(ContrepoisPost1h,Post1hMeans)
ContrepoisPost2h <- union(ContrepoisPost2h,Post2hMeans)
ContrepoisPost4h <- union(ContrepoisPost4h,Post4hMeans)
ContrepoisPost6h <- union(ContrepoisPost6h,Post6hMeans)
ContrepoisPost24h <- union(ContrepoisPost24h,Post24hMeans)


BaselinePairIPE <- as.data.frame(t(BaselinePairIPE))
BaselinePair15m <- as.data.frame(t(BaselinePair15m))
BaselinePair30m <- as.data.frame(t(BaselinePair30m))
BaselinePair1h <- as.data.frame(t(BaselinePair1h))
BaselinePair2h <- as.data.frame(t(BaselinePair2h))
BaselinePair4h <- as.data.frame(t(BaselinePair4h))
BaselinePair6h <- as.data.frame(t(BaselinePair6h))
BaselinePair24h <- as.data.frame(t(BaselinePair24h))


ContrepoisPost15m <- as.data.frame(t(ContrepoisPost15m))
ContrepoisPost30m <- as.data.frame(t(ContrepoisPost30m))
ContrepoisPost1h <- as.data.frame(t(ContrepoisPost1h))
ContrepoisPost2h <- as.data.frame(t(ContrepoisPost2h))
ContrepoisPost4h <- as.data.frame(t(ContrepoisPost4h))
ContrepoisPost6h <- as.data.frame(t(ContrepoisPost6h))
ContrepoisPost24h <- as.data.frame(t(ContrepoisPost24h))
ContrepoisIPE <- as.data.frame(t(ContrepoisIPE))

ContrepoisPairedOutputIPE$log2_fold_change <- ContrepoisIPE$Mean - BaselinePairIPE$Mean


write.csv(ContrepoisPairedOutputIPE, "C:/Users/Abrraski/OneDrive - Michigan Medicine/Desktop/DifferentialExpressionAnalysis/Contrepois_IPE_UntargetedDifferentrialAnalysis.csv",row.names =T)
