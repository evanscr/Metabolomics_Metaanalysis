rm(list=ls())
#* Now that we can install packages, we have to load them into our working memory each time we open a new
#* R ression
setwd('C:/Users/Abrraski/OneDrive - Michigan Medicine/Desktop/DifferentialExpressionAnalysis')
getwd()
library(DESeq2)
library(tidyverse)
###############
#read in the data set
FUEL2016 <- read.csv('C:/Users/Abrraski/OneDrive - Michigan Medicine/Desktop/DifferentialExpressionAnalysis/FUEL2016_PostSERRFNormalization_Data_For_DifferentialAnalysis.csv')
FUEL2016Meta <- FUEL2016[,c(1:4)]
FUEL2016 <- FUEL2016[,-c(1:4)]
#Extract Metadata and metabolites, and get rid of empty columns
FuelMetabolites <- as.data.frame(colnames(FUEL2016))
rownames(FuelMetabolites) <- FuelMetabolites[,1]
rownames(FUEL2016) <- FUEL2016Meta$File
#remove metadata from dataframe and place samples in columns and features in rows
#Check class of frame to ensure numeric
class(FUEL2016[,163])
any(is.na(FUEL2016))
#for some reason It was filtering out way too many features with NA's that shouldn't have been, so I replaced the NA with
#Zeroes, then filtered
FUEL2016 <- FUEL2016 %>% replace(is.na(.), 0)
FUEL2016<- as.data.frame(t(FUEL2016))
Fuel2016Filtered <- FUEL2016 %>% filter ((rowSums((. > 0)) / ncol(.)) > 0.3, na.rm = TRUE)

print(rowSums(FUEL2016[5,]))
#Now I'm going to swtich them back and replace the NA values with 1/5 the min of the row
Fuel2016Filtered[Fuel2016Filtered==0] <- NA
for (i in 1:nrow(Fuel2016Filtered)){
  Fuel2016Filtered[i,][is.na(Fuel2016Filtered[i,])] <- 0.2 * min(Fuel2016Filtered[i,][!(is.na(Fuel2016Filtered[i,]))])
}
#Now I'm going to log transform the data
Fuel2016Filtered <- as.data.frame(apply(Fuel2016Filtered,2, log2))
#We can scale the data and look at some histograms real quick
Fuel2016Filtered <- as.data.frame(t(Fuel2016Filtered))
FUEL2016Meta <- as.data.frame(t(FUEL2016Meta))
FUEL2016Meta <- as.data.frame(t(FUEL2016Meta))
Fuel2016FilteredScaled <- as.data.frame(apply(Fuel2016Filtered,1, function(x) (x-mean(x))/sd(x)))
Fuel2016FilteredScaled <- as.data.frame(t(Fuel2016FilteredScaled))
#Histograms
z = FuelMetabolites[6,1]
ggplot(Fuel2016Filtered, aes(x =`UNK_4.142_185.13203`))+geom_histogram()

ggplot(Fuel2016FilteredScaled, aes(x=`UNK_4.142_185.13203`))+geom_histogram()

#Add metadata and separate by condition
Fuel2016Filtered$Id <- FUEL2016Meta$Subject
Fuel2016FilteredScaled$Id <- FUEL2016Meta$Subject
Fuel2016Filtered$Timepoint <- FUEL2016Meta$TimepointAdjusted
Fuel2016FilteredScaled$Timepoint <- FUEL2016Meta$TimepointAdjusted

differences <- with(Fuel2016Filtered, UNK_4.142_185.13203[Timepoint == "IPE"] - UNK_4.142_185.13203[Timepoint == "Baseline"])
shapiro.test(differences)

differencesScaled <- with(Fuel2016FilteredScaled, UNK_4.142_185.13203[Timepoint == "IPE"] - UNK_4.142_185.13203[Timepoint == "Baseline"])
shapiro.test(differencesScaled)



#Separate the data first by exercise type, then by sample collection timepoint
Fuel2016Baseline <- Fuel2016Filtered %>% filter(Timepoint == 'Baseline')
Fuel2016IPE <- Fuel2016Filtered %>% filter(Timepoint == 'IPE')
Fuel2016Post3m <- Fuel2016Filtered %>% filter(Timepoint == 'Post 3')
Fuel2016Post6m <- Fuel2016Filtered %>% filter(Timepoint == 'Post 6')
Fuel2016T3 <- Fuel2016Filtered %>% filter(Timepoint == 'T3' | Timepoint == 'T03', Id != "Pre-9")
Fuel2016T6 <- Fuel2016Filtered %>% filter(Timepoint == 'T6' | Timepoint == 'T06', Id != "Pre-9")
Fuel2016T9 <- Fuel2016Filtered %>% filter(Timepoint == 'T9' | Timepoint == 'T09', Id != "Pre-9")
Fuel2016T12 <- Fuel2016Filtered %>% filter(Timepoint == 'T12', Id != "Pre-9")
Fuel2016T15 <- Fuel2016Filtered %>% filter(Timepoint == 'T15', Id != "Pre-9")
Fuel2016T18 <- Fuel2016Filtered %>% filter(Timepoint == 'T18', Id != "Pre-9")
Fuel2016T21 <- Fuel2016Filtered %>% filter(Timepoint == 'T21', Id != "Pre-9")
Fuel2016T24 <- Fuel2016Filtered %>% filter(Timepoint == 'T24', Id != "Pre-9")
Fuel2016T27 <- Fuel2016Filtered %>% filter(Timepoint == 'T27', Id != "Pre-9")
Fuel2016T30 <- Fuel2016Filtered %>% filter(Timepoint == 'T30', Id != "Pre-9")
Fuel2016T33 <- Fuel2016Filtered %>% filter(Timepoint == 'T33', Id != "Pre-9")
Fuel2016T36 <- Fuel2016Filtered %>% filter(Timepoint == 'T36', Id != "Pre-9")
Fuel2016T39 <- Fuel2016Filtered %>% filter(Timepoint == 'T39', Id != "Pre-9")
Fuel2016T42 <- Fuel2016Filtered %>% filter(Timepoint == 'T42', Id != "Pre-9")
Fuel2016T45 <- Fuel2016Filtered %>% filter(Timepoint == 'T45', Id != "Pre-9")


ReorderIPE <- match(Fuel2016IPE$Id, Fuel2016Baseline$Id)
BaselinePairIPE <- Fuel2016Baseline[ReorderIPE,]

ReorderT3 <- match(Fuel2016T3$Id, Fuel2016Baseline$Id)
BaselinePairT3 <- Fuel2016Baseline[ReorderT3,]

ReorderT6 <- match(Fuel2016T6$Id, Fuel2016Baseline$Id)
BaselinePairT6 <- Fuel2016Baseline[ReorderT6,]

ReorderT9 <- match(Fuel2016T9$Id, Fuel2016Baseline$Id)
BaselinePairT9 <- Fuel2016Baseline[ReorderT9,]

ReorderT12 <- match(Fuel2016T12$Id, Fuel2016Baseline$Id)
BaselinePairT12 <- Fuel2016Baseline[ReorderT12,]

ReorderT15 <- match(Fuel2016T15$Id, Fuel2016Baseline$Id)
BaselinePairT15 <- Fuel2016Baseline[ReorderT15,]

ReorderT18 <- match(Fuel2016T18$Id, Fuel2016Baseline$Id)
BaselinePairT18 <- Fuel2016Baseline[ReorderT18,]

ReorderT21 <- match(Fuel2016T21$Id, Fuel2016Baseline$Id)
BaselinePairT21 <- Fuel2016Baseline[ReorderT21,]

ReorderT24 <- match(Fuel2016T24$Id, Fuel2016Baseline$Id)
BaselinePairT24 <- Fuel2016Baseline[ReorderT24,]

ReorderT27 <- match(Fuel2016T27$Id, Fuel2016Baseline$Id)
BaselinePairT27 <- Fuel2016Baseline[ReorderT27,]

ReorderT30 <- match(Fuel2016T30$Id, Fuel2016Baseline$Id)
BaselinePairT30 <- Fuel2016Baseline[ReorderT30,]

ReorderT33 <- match(Fuel2016T33$Id, Fuel2016Baseline$Id)
BaselinePairT33 <- Fuel2016Baseline[ReorderT33,]

ReorderPost3 <- match(Fuel2016Post3m$Id, Fuel2016Baseline$Id)
BaselinePairPost3 <- Fuel2016Baseline[ReorderPost3,]

ReorderPost6 <- match(Fuel2016Post6m$Id, Fuel2016Baseline$Id)
BaselinePairPost6 <- Fuel2016Baseline[ReorderPost6,]

all(BaselinePairT27$Id == Fuel2016T27$Id)
#Scrub MetaData and perform t test
BaselinePairIPE <- BaselinePairIPE %>% select(!(Id:Timepoint))
BaselinePairT3 <- BaselinePairT3 %>% select(!(Id:Timepoint))
BaselinePairT6 <- BaselinePairT6 %>% select(!(Id:Timepoint))
BaselinePairT9 <- BaselinePairT9 %>% select(!(Id:Timepoint))
BaselinePairT12 <- BaselinePairT12 %>% select(!(Id:Timepoint))
BaselinePairT15 <- BaselinePairT15 %>% select(!(Id:Timepoint))
BaselinePairT18 <- BaselinePairT18 %>% select(!(Id:Timepoint))
BaselinePairT21 <- BaselinePairT21 %>% select(!(Id:Timepoint))
BaselinePairT24 <- BaselinePairT24 %>% select(!(Id:Timepoint))
BaselinePairT27 <- BaselinePairT27 %>% select(!(Id:Timepoint))
BaselinePairT30 <- BaselinePairT30 %>% select(!(Id:Timepoint))
BaselinePairT33 <- BaselinePairT33 %>% select(!(Id:Timepoint))
BaselinePairPost3 <- BaselinePairPost3 %>% select(!(Id:Timepoint))
BaselinePairPost6 <- BaselinePairPost6 %>% select(!(Id:Timepoint))

Fuel2016IPE <- Fuel2016IPE %>% select(!(Id:Timepoint))
Fuel2016T3 <- Fuel2016T3 %>% select(!(Id:Timepoint))
Fuel2016T6 <- Fuel2016T6 %>% select(!(Id:Timepoint))
Fuel2016T9 <- Fuel2016T9 %>% select(!(Id:Timepoint))
Fuel2016T12 <- Fuel2016T12 %>% select(!(Id:Timepoint))
Fuel2016T15 <- Fuel2016T15 %>% select(!(Id:Timepoint))
Fuel2016T18 <- Fuel2016T18 %>% select(!(Id:Timepoint))
Fuel2016T21 <- Fuel2016T21 %>% select(!(Id:Timepoint))
Fuel2016T24 <- Fuel2016T24 %>% select(!(Id:Timepoint))
Fuel2016T27 <- Fuel2016T27 %>% select(!(Id:Timepoint))
Fuel2016T30 <- Fuel2016T30 %>% select(!(Id:Timepoint))
Fuel2016T33 <- Fuel2016T33 %>% select(!(Id:Timepoint))
Fuel2016Post3m <- Fuel2016Post3m %>% select(!(Id:Timepoint))
Fuel2016Post6m <- Fuel2016Post6m %>% select(!(Id:Timepoint))




#Now we can run the t-test
Fuel2016OutputT33 <- data.frame(metabolites = colnames(BaselinePairIPE),
                            t_statistic = rep(NA,length(colnames(BaselinePairIPE))),
                            p_value = rep(NA,length(colnames(BaselinePairIPE))))
# we need to remove the groups column since we cant run a t.test on that

for (i in 1:ncol(BaselinePairIPE)){
  
  t_test<- t.test(BaselinePairT33[,i], Fuel2016T33[,i], var.equal = FALSE, paired = TRUE)
  Fuel2016OutputT33$t_statistic[i] <- t_test$statistic
  Fuel2016OutputT33$p_value[i] <- t_test$p.value
}

#lets adjust the pvalues

Fuel2016OutputPost6m$adj_pvalues <- p.adjust(Fuel2016OutputPost6m$p_value, method = 'BH')


IPEBaselineMeans <- as.data.frame(t(colMeans(BaselinePairIPE))) 
T3BaselineMeans <- as.data.frame(t(colMeans(BaselinePairT3)))
T6BaselineMeans <- as.data.frame(t(colMeans(BaselinePairT6)))
T9BaselineMeans <- as.data.frame(t(colMeans(BaselinePairT9)))
T12BaselineMeans <- as.data.frame(t(colMeans(BaselinePairT12)))
T15BaselineMeans <- as.data.frame(t(colMeans(BaselinePairT15)))
T18BaselineMeans <- as.data.frame(t(colMeans(BaselinePairT18)))
T21BaselineMeans <- as.data.frame(t(colMeans(BaselinePairT21)))
T24BaselineMeans <- as.data.frame(t(colMeans(BaselinePairT24)))
T27BaselineMeans <- as.data.frame(t(colMeans(BaselinePairT27)))
T30BaselineMeans <- as.data.frame(t(colMeans(BaselinePairT30)))
T33BaselineMeans <- as.data.frame(t(colMeans(BaselinePairT33)))
Post3BaselineMeans <- as.data.frame(t(colMeans(BaselinePairPost3)))
Post6BaselineMeans <- as.data.frame(t(colMeans(BaselinePairPost6)))

IPEMeans <- as.data.frame(t(colMeans(Fuel2016IPE)))
T3Means <- as.data.frame(t(colMeans(Fuel2016T3)))
T6Means <- as.data.frame(t(colMeans(Fuel2016T6)))
T9Means <- as.data.frame(t(colMeans(Fuel2016T9)))
T12Means <- as.data.frame(t(colMeans(Fuel2016T12)))
T15Means <- as.data.frame(t(colMeans(Fuel2016T15)))
T18Means <- as.data.frame(t(colMeans(Fuel2016T18)))
T21Means <- as.data.frame(t(colMeans(Fuel2016T21)))
T24Means <- as.data.frame(t(colMeans(Fuel2016T24)))
T27Means <- as.data.frame(t(colMeans(Fuel2016T27)))
T30Means <- as.data.frame(t(colMeans(Fuel2016T30)))
T33Means <- as.data.frame(t(colMeans(Fuel2016T33)))
Post3Means <- as.data.frame(t(colMeans(Fuel2016Post3m)))
Post6Means <- as.data.frame(t(colMeans(Fuel2016Post6m)))

rownames(IPEBaselineMeans) <-  'Mean'
rownames(T3BaselineMeans) <-  'Mean'
rownames(T6BaselineMeans) <-  'Mean'
rownames(T9BaselineMeans) <-  'Mean'
rownames(T12BaselineMeans) <-  'Mean'
rownames(T15BaselineMeans) <-  'Mean'
rownames(T18BaselineMeans) <-  'Mean'
rownames(T21BaselineMeans) <-  'Mean'
rownames(T24BaselineMeans) <-  'Mean'
rownames(T27BaselineMeans) <-  'Mean'
rownames(T30BaselineMeans) <-  'Mean'
rownames(T33BaselineMeans) <-  'Mean'
rownames(Post3BaselineMeans) <-  'Mean'
rownames(Post6BaselineMeans) <-  'Mean'
rownames(IPEMeans) <-  'Mean'
rownames(T3Means) <-  'Mean'
rownames(T6Means) <-  'Mean'
rownames(T9Means) <-  'Mean'
rownames(T12Means) <-  'Mean'
rownames(T15Means) <-  'Mean'
rownames(T18Means) <-  'Mean'
rownames(T21Means) <-  'Mean'
rownames(T24Means) <-  'Mean'
rownames(T27Means) <-  'Mean'
rownames(T30Means) <-  'Mean'
rownames(T33Means) <-  'Mean'
rownames(Post3Means) <-  'Mean'
rownames(Post6Means) <-  'Mean'

Fuel2016IPE <- union(Fuel2016IPE, IPEMeans)
Fuel2016T3 <- union(Fuel2016T3, T3Means)
Fuel2016T6 <- union(Fuel2016T6, T6Means)
Fuel2016T9 <- union(Fuel2016T9, T9Means)
Fuel2016T12 <- union(Fuel2016T12, T12Means)
Fuel2016T15 <- union(Fuel2016T15, T15Means)
Fuel2016T18 <- union(Fuel2016T18, T18Means)
Fuel2016T21 <- union(Fuel2016T21, T21Means)
Fuel2016T24 <- union(Fuel2016T24, T24Means)
Fuel2016T27 <- union(Fuel2016T27, T27Means)
Fuel2016T30 <- union(Fuel2016T30, T30Means)
Fuel2016T33 <- union(Fuel2016T33, T33Means)
Fuel2016Post3m <- union(Fuel2016Post3m, Post3Means)
Fuel2016Post6m <- union(Fuel2016Post6m, Post6Means)

BaselinePairIPE <- union(BaselinePairIPE,IPEBaselineMeans)
BaselinePairT3 <- union(BaselinePairT3,T3BaselineMeans)
BaselinePairT6 <- union(BaselinePairT6,T6BaselineMeans)
BaselinePairT9 <- union(BaselinePairT9,T9BaselineMeans)
BaselinePairT12 <- union(BaselinePairT12,T12BaselineMeans)
BaselinePairT15 <- union(BaselinePairT15,T15BaselineMeans)
BaselinePairT18 <- union(BaselinePairT18,T18BaselineMeans)
BaselinePairT21 <- union(BaselinePairT21,T21BaselineMeans)
BaselinePairT24 <- union(BaselinePairT24,T24BaselineMeans)
BaselinePairT27 <- union(BaselinePairT27,T27BaselineMeans)
BaselinePairT30 <- union(BaselinePairT30,T30BaselineMeans)
BaselinePairT33 <- union(BaselinePairT33,T33BaselineMeans)
BaselinePairPost3 <- union(BaselinePairPost3,Post3BaselineMeans)
BaselinePairPost6 <- union(BaselinePairPost6,Post6BaselineMeans)


Fuel2016IPE <- as.data.frame(t(Fuel2016IPE))
Fuel2016T3 <- as.data.frame(t(Fuel2016T3))
Fuel2016T6 <- as.data.frame(t(Fuel2016T6))
Fuel2016T9 <- as.data.frame(t(Fuel2016T9))
Fuel2016T12 <- as.data.frame(t(Fuel2016T12))
Fuel2016T15 <- as.data.frame(t(Fuel2016T15))
Fuel2016T18 <- as.data.frame(t(Fuel2016T18))
Fuel2016T21 <- as.data.frame(t(Fuel2016T21))
Fuel2016T24 <- as.data.frame(t(Fuel2016T24))
Fuel2016T27 <- as.data.frame(t(Fuel2016T27))
Fuel2016T30 <- as.data.frame(t(Fuel2016T30))
Fuel2016T33 <- as.data.frame(t(Fuel2016T33))
Fuel2016Post3m <- as.data.frame(t(Fuel2016Post3m))
Fuel2016Post6m <- as.data.frame(t(Fuel2016Post6m))

BaselinePairIPE <- as.data.frame(t(BaselinePairIPE))
BaselinePairT3 <- as.data.frame(t(BaselinePairT3))
BaselinePairT6 <- as.data.frame(t(BaselinePairT6))
BaselinePairT9 <- as.data.frame(t(BaselinePairT9))
BaselinePairT12 <- as.data.frame(t(BaselinePairT12))
BaselinePairT15 <- as.data.frame(t(BaselinePairT15))
BaselinePairT18 <- as.data.frame(t(BaselinePairT18))
BaselinePairT21 <- as.data.frame(t(BaselinePairT21))
BaselinePairT24 <- as.data.frame(t(BaselinePairT24))
BaselinePairT27 <- as.data.frame(t(BaselinePairT27))
BaselinePairT30 <- as.data.frame(t(BaselinePairT30))
BaselinePairT33 <- as.data.frame(t(BaselinePairT33))
BaselinePairPost3 <- as.data.frame(t(BaselinePairPost3))
BaselinePairPost6 <- as.data.frame(t(BaselinePairPost6))


Fuel2016OutputPost6m$log2_fold_change <- Fuel2016Post6m$Mean - BaselinePairPost6$Mean


write.csv(Fuel2016OutputPost6m, "C:/Users/Abrraski/OneDrive - Michigan Medicine/Desktop/DifferentialExpressionAnalysis/Fuel2016_DifferentialAnalysis_Post6.csv",row.names =F)


