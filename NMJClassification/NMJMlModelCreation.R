library(caret)
library(tidyverse)
library(skimr)
library(ggplot2)
library(RANN) 
library(plotROC)
library(MLeval)
set.seed(1337)

#########################################################################################################################
#########################################Initialise Data#################################################################
#########################################################################################################################

#Load in initial data
MLDataset <- read_csv("AlanMLDatabaseTrueFinal.csv")
#Rename All columns and remove some unwanted values
colnames(MLDataset) <- c("Sample", "Muscle", "Fenotype", "Genotype", "Timepoint", "Maturity",
                         "NMJCounting", "Bounds0Red", "Bounds1Red", "Bounds2Red", 
                         "CCDistRed", "CCSize0Red", "CCSize1Red", "CCSize2Red", 
                         "CCRed", "FragmentationRed", "Compactness0Red", "Compactness1Red", 
                         "MeanIntensityRed","ShapeFactor0Red", "ShapeFactor1Red", "ShapeFactor2Red", 
                         "SkeletonLengthRed", "SurfaceVolumeRatio0Red", "SurfaceVolumeRatio1Red", "Surface0Red", 
                         "Surface1Red", "Surface2Red", "Surface3Red", "SurfaceDil0Red", 
                         "SurfaceDil1Red", "Volume0Red", "Volume1Red",  "Bounds0Gr", 
                         "Bounds1Gr", "Bounds2Gr", "CCDistGr", "CCSize0Gr", 
                         "CCSize1Gr",  "CCSize2Gr", "CCGr", "Compactness0Gr",  
                         "Compactness1Gr", "FragmentationGr", "ShapeFactor0Gr", "ShapeFactor1Gr", 
                         "ShapeFactor2Gr", "SkeletonLengthGr", "SurfaceVolumeRatio0Gr", "SurfaceVolumeRatio1Gr", 
                         "Surface0Gr", "Surface1Gr", "Surface2Gr", "Surface3Gr", 
                         "SurfaceDil0Gr", "SurfaceDil1Gr", "Volume0Gr", "Volume1Gr", 
                         "IoUGr", "AveDist", "ComDist", "Coverage", 
                         "HausDist", "NIntersection", "NUnion")

MLDataset <- MLDataset[rowSums(is.na(MLDataset)) != ncol(MLDataset),] %>%
  select(-Sample, -Timepoint) %>% select( NMJCounting, everything()) %>%
  filter(Muscle != "MUSCLE", Fenotype != "FENOTYPE")

#Change Variables from Text so caret can parse them
MLDatasetChanged <- MLDataset
MLDatasetChanged <- drop_na(MLDatasetChanged, NMJCounting)

#Drop average dist column and merge the two types of unhealthy NMJs
MLDatasetChangeNA <- MLDatasetChanged %>% select(-AveDist, -Muscle, -Fenotype, -Genotype, -Maturity) %>%
  drop_na()

MLDatasetChangeNAMerge <- MLDatasetChangeNA %>% 
  mutate(NMJCounting = replace(NMJCounting, NMJCounting == "Denervated" | NMJCounting == "Partial" , "Degenerating")) %>%
  mutate(NMJCounting = replace(NMJCounting, NMJCounting == "Fully", "Healthy")) 

MLDatasetTrimmed <- MLDatasetChangeNAMerge

#########################################################################################################################
#########################################Create ML Models################################################################
#########################################################################################################################

#Create a ML dataset without k-fold validation
#Partition the data into training and test data
trainRowNumbers <- createDataPartition(MLDatasetTrimmed$NMJCounting, p=0.8, list=FALSE)

trainData <- MLDatasetTrimmed[trainRowNumbers,]

testData <- MLDatasetTrimmed[-trainRowNumbers,]

##########################################################################################################################

#Format data for filtering
x = select(trainData, -NMJCounting)
y = trainData$NMJCounting

trainData$NMJCounting <- y

##########################################################################################################################

#Fit the ML model with k-fold validation
RFFitNoKFold <- train(as.factor(NMJCounting) ~ ., 
                                   data = trainData, 
                                   method = "ranger")

RFPredNoKFold <- predict(RFFitNoKFold, testData) 
# compare predicted outcome and true outcome
confusionMatrix(RFPredNoKFold, as.factor(testData$NMJCounting))


saveRDS(RFFitNoKFoldNoNormalParam, file="RFFit.rds")

#########################################################################################################################

#Create a ML dataset with k-fold validation
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)


RFFitKFold <- train(as.factor(NMJCounting) ~ ., 
                                 data = trainData, 
                                 method = "ranger",
                                 trControl = fitControl)

RFPredKFold <- predict(RFFitKFold, testData) 
confusionMatrix(RFPredKFold, as.factor(testDataDropped$NMJCounting))

saveRDS(RFFitKFoldNoNormalParam, file="RFFitKFold.rds")

#########################################################################################################################

# Create a ML dataset that saves predictions so AUC can be used 
trainData <- select(trainData, -ShapeFactor1Gr)
ctrl <- trainControl(method="cv", 
                     summaryFunction=twoClassSummary, 
                     classProbs=T,
                     savePredictions = T)
rfFitWctrl <- train(NMJCounting ~ ., data=trainData, 
                    method="rf", preProc=c("center", "scale"), 
                    trControl=ctrl)

res <- evalm(rfFitWctrl,gnames='rf')


saveRDS(rfFitWctrl, file="RFFitKFoldAUC.rds")

#########################################################################################################################
##########################################Equalise the Datasets##########################################################
#########################################################################################################################

MLDatasetTrimmedEqualised <- MLDatasetTrimmed %>% 
  group_by(NMJCounting) %>%
  do(sample_n(.,1000))

trainRowNumbersEqualised <- createDataPartition(MLDatasetTrimmedEqualised$NMJCounting, p=0.8, list=FALSE)

trainDataEqualised <- MLDatasetTrimmedEqualised[trainRowNumbersEqualised,]

testDataEqualised <- MLDatasetTrimmedEqualised[-trainRowNumbersEqualised,]

#Format data for filtering
xEqualised = select(trainDataEqualised, -NMJCounting)
yEqualised = trainDataEqualised$NMJCounting

trainDataEqualised$NMJCounting <- yEqualised


RFFitNoKFoldEqualised <- train(as.factor(NMJCounting) ~ ., 
                               data = trainDataEqualised, 
                               method = "ranger")


RFPredNoKFoldEqualised <- predict(RFFitNoKFoldEqualised, testDataEqualised) 
# compare predicted outcome and true outcome
confusionMatrix(RFPredNoKFoldEqualised, as.factor(testDataEqualised$NMJCounting))

saveRDS(RFFitNoKFoldEqualised, file="NMJRFFitEqualised.rds")

#########################################################################################################################

fitControl <- trainControl(## 10-fold CV)
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)

RFFitKFoldEqualised <- train(as.factor(NMJCounting) ~ ., 
                             data = trainDataEqualised, 
                             method = "ranger",
                             trControl = fitControl)

RFPredKFoldEqualised <- predict(RFFitKFoldEqualised, testDataEqualised) 
confusionMatrix(RFPredKFoldEqualised, as.factor(testDataEqualised$NMJCounting))


saveRDS(RFFitKFoldEqualised, file="NMJRFKfoldFitEqualised.rds")