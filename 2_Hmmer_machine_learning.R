
#Load the required R libraries
getwd()
library(ggplot2)
library(caret)
library(factoextra)
library(Cairo)
library(randomForest)
library(NeuralNetTools)
library(gtools)


#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################


#Open a new window to store plots
x11()


#Bacterial genera that are reduced in neuropsychiatric disorders
reduced_taxa <- c("Akkermansia", "Faecalibacterium", "Lactobacillus",
                    "Prevotella", "Dialister", "Coprococcus")


#Bacterial genera that are increased in neuropsychiatric disorders
increased_taxa <- c("Bacteroides", "Alkaliflexus", "Sutterella",
                  "Haemophilus", "Eggerthella", "Holdemania",
                  "Turicibacter", "Paraprevotella", "Anaerofilm")


#Read R objects that will be used as inputs for machine learning models
domain_data_numeric <- readRDS("domain_data_numeric.rds")
#Transpose the data frame
domain_data_numeric <- as.data.frame(t(domain_data_numeric))
#Add taxa information for each sample as a label
domain_data_numeric$Taxa <- as.factor(gsub("_GCA.*", "", rownames(domain_data_numeric)))
#Reorder the data frame so taxa label corresponds to the first column
domain_data_numeric <- domain_data_numeric[,c(ncol(domain_data_numeric), 
                                                     1:(ncol(domain_data_numeric)-1))]


#Machine learning models will discriminated between taxa increased and reduce in neuropsychiatric disorders
#Rename those taxa that are reduced in neuropsychiatric disorders as "Reduced"
domain_data_numeric$Taxa <- gsub(paste0(increased_taxa, collapse="|"), 
                                      "Increased", domain_data_numeric$Taxa)
#Rename those taxa that are increased in neuropsychiatric disorders as "Increased"
domain_data_numeric$Taxa <- gsub(paste0(reduced_taxa, collapse="|"), 
                                      "Reduced", domain_data_numeric$Taxa)
#Set the label as factor
domain_data_numeric$Taxa <- factor(domain_data_numeric$Taxa)


#Prepare the different inputs for supervised machine learning models
#Generate an empty list to store train and test data sets
train_test_datasets <- list()
#Generate empty data frames to store train and test variables and classes
train_variables <- data.frame()
train_labels <- data.frame()
test_variables <- data.frame()
test_labels <- data.frame()


#Iterate all levels from the factor variable (i.e. microbial taxa)
for (i in 1:length(levels(domain_data_numeric$Taxa))) {
  
  #Subset samples for i level
  i_level <- domain_data_numeric[which(domain_data_numeric$Taxa == levels(domain_data_numeric$Taxa)[i]),]
  #Randomly select 70% samples from i_level
  set.seed(600)
  i_random_selection <- sample(1:nrow(i_level),round(0.3*nrow(i_level)))
  #Generate train and test sets from 70 and 30% of samples, respectively
  i_train <- i_level[-i_random_selection,]
  i_test <- i_level[i_random_selection,]
  
  #Select train variables
  i_train_variables <- i_train[,2:ncol(i_train)]
  #Select train label (i.e. microbial taxa)
  i_train_label <- factor(i_train[,1], levels=unique(i_train[,1]))
  
  #Select train variables
  i_test_variables <- i_test[,2:ncol(i_test)]
  #Select train class
  i_test_label <- factor(i_test[,1], levels=unique(i_test[,1]))
  
  #Append train and test data sets to a list
  i_list <- list(i_train_variables, i_train_label, 
                 i_test_variables, i_test_label)
  names(i_list) <- c("train_variables", "train_label", 
                     "test_variables", "test_label")
  
  #Append this new list to the list that was previously created before the loop
  train_test_datasets[[i]] <- i_list
  names(train_test_datasets)[length(train_test_datasets)] <- levels(domain_data_numeric$Taxa)[i]
  
  #Append train and test variables corresponding to i taxa to the data frames that were previously created before the loop
  train_variables <- rbind.data.frame(train_variables, i_train_variables)
  test_variables <- rbind.data.frame(test_variables, i_test_variables)
  
  #Append train and test labels corresponding to i taxa to the data frames that were previously created before the loop
  train_labels <- unlist(list(train_labels, i_train_label))
  test_labels <- unlist(list(test_labels, i_test_label))
  
}


#Prune near zero variance predictors
#Not necessary: neurotransmitter-related functional domains were already selected by manual curation in the first script
#train_variables <- train_variables[,-nearZeroVar(train_variables)]
#test_variables <- test_variables[,names(train_variables)]


#Predictors can also be pruned according to their importance
# filter_variable_importance <- filterVarImp(train_variables, train_labels)
# most_important_variables <- which(rowMeans(filter_variable_importance) > 0.65)
#Prune train and test set considering only most important variables
# train_variables <- train_variables[,names(most_important_variables)]
# test_variables <- test_variables[,names(train_variables)]




#Perform a principal components analysis (PCA) on the train set
#Prepare data for PCA model
pca_variables <- rbind.data.frame(train_variables, test_variables)
names(pca_variables) <- gsub("_", " ", names(pca_variables))
pca_labels <- unlist(list(train_labels, test_labels))
#This initial PCA model reveals some differences between clades increased and reduced in neuropsychiatric disorders
#PCA plot showing variables names
pca_plot_with_var_names <- fviz_pca_biplot(prcomp(pca_variables), label="var", 
                            habillage=pca_labels,
                            addEllipses=T, ellipse.level=0.95, 
                            pointsize=4.4, labelsize=5, arrowsize=1) + 
                            theme(text=element_text(size=20))
#PCA plot without variables names
pca_plot_without_var_names <- fviz_pca_biplot(prcomp(pca_variables), label="var", 
                                              habillage=pca_labels,
                                              addEllipses=T, ellipse.level=0.95, 
                                              pointsize=4.4, labelsize=0, arrowsize=1) + 
                                              theme(text=element_text(size=20))
#Format axis titles
pca_plot_with_var_names$labels$x <- gsub("Dim1", "PC 1", pca_plot_with_var_names$labels$x)
pca_plot_with_var_names$labels$y <- gsub("Dim2", "PC 2", pca_plot_with_var_names$labels$y)
pca_plot_without_var_names$labels$x <- gsub("Dim1", "PC 1", pca_plot_without_var_names$labels$x)
pca_plot_without_var_names$labels$y <- gsub("Dim2", "PC 2", pca_plot_without_var_names$labels$y)
#Export high quality figures
ggsave(file="PCA_plot_with_var_names.png", pca_plot_with_var_names, device='png', dpi=600) 
ggsave(file="PCA_plot_without_var_names.png", pca_plot_without_var_names, device='png', dpi=600) 




#Perform supervised learning
#Train a partial least squares (pls) model
#Results: samples are not properly discriminated, profiles of neurotransmitter-related domains are very similar among groups
set.seed(600) 
pls <- train(train_variables, train_labels, method="pls", 
             tuneLength=15)


#Calculate the train accuracy rate
pls_predictions_train <- predict(pls, train_variables)
pls_confusion_matrix_train <- confusionMatrix(as.factor(train_labels), pls_predictions_train)
#Calculate the test accuracy rate
pls_predictions_test <- predict(pls, test_variables)
pls_confusion_matrix_test <- confusionMatrix(as.factor(test_labels), pls_predictions_test)
#Export pls model as an R object
saveRDS(pls, "pls.rds")
#Read pls model stored in an R object
pls <- readRDS("pls.rds")




#Train a boosted logistic regression model (LogitBoost) model
#Results: samples are not properly discriminated, profiles of neurotransmitter-related domains are very similar among groups
set.seed(600) 
logitboost <- train(train_variables, train_labels, method="LogitBoost", 
                    tuneLength=5)


#Calculate the train accuracy rate
logitboost_predictions_train <- predict(logitboost, train_variables)
logitboost_confusion_matrix_train <- confusionMatrix(as.factor(train_labels), logitboost_predictions_train)
#Calculate the test accuracy rate
logitboost_predictions_test <- predict(logitboost, test_variables)
logitboost_confusion_matrix_test <- confusionMatrix(as.factor(test_labels), logitboost_predictions_test)
#Export pls model as an R object
saveRDS(logitboost, "logitboost.rds")
#Read pls model stored in an R object
logitboost <- readRDS("logitboost.rds")




#Train a generalized linear model elastic-net (glmnet) model
#Results: samples are not properly discriminated, profiles of neurotransmitter-related domains are very similar among groups
set.seed(600)
glmnet <- train(train_variables, train_labels, method="glmnet", 
                tuneLength=2)


#Calculate the train accuracy rate
glmnet_predictions_train <- predict(glmnet, train_variables)
glmnet_confusion_matrix_train <- confusionMatrix(as.factor(train_labels), glmnet_predictions_train)
#Calculate the test accuracy rate
glmnet_predictions_test <- predict(glmnet, test_variables)
glmnet_confusion_matrix_test <- confusionMatrix(as.factor(test_labels), glmnet_predictions_test)
#Export pls model as an R object
saveRDS(glmnet, "glmnet.rds")
#Read pls model stored in an R object
glmnet <- readRDS("glmnet.rds")




#Train a random forest model (rf) model
#Results: samples are not properly discriminated, profiles of neurotransmitter-related domains are very similar among groups
set.seed(600)
rf <- train(train_variables, train_labels, method="rf", 
            importance=T,
            tuneLength=2)


#Calculate the train accuracy rate
rf_predictions_train <- predict(rf, train_variables)
rf_confusion_matrix_train <- confusionMatrix(as.factor(train_labels), rf_predictions_train)
#Calculate the test accuracy rate
rf_predictions_test <- predict(rf, test_variables)
rf_confusion_matrix_test <- confusionMatrix(as.factor(test_labels), rf_predictions_test)
#Export pls model as an R object
saveRDS(rf, "rf.rds")
#Read pls model stored in an R object
rf <- readRDS("rf.rds")




#Train an artificial neural network model (ann) model
#Results: samples are not properly discriminated, profiles of neurotransmitter-related domains are very similar among groups
set.seed(600)
ann <- train(train_variables, train_labels, method="mlp", 
            tuneLength=2)


#Calculate the train accuracy rate
ann_predictions_train <- predict(ann, train_variables)
ann_confusion_matrix_train <- confusionMatrix(as.factor(train_labels), ann_predictions_train)
#Calculate the test accuracy rate
ann_predictions_test <- predict(ann, test_variables)
ann_confusion_matrix_test <- confusionMatrix(as.factor(test_labels), ann_predictions_test)
#Export pls model as an R object
saveRDS(ann, "ann.rds")
#Read pls model stored in an R object
ann <- readRDS("ann.rds")




#Compare the performance metrics of different machine learning models studied
#Selected parameters: Accuracy and Kappa (instead of ROC, sensitivity and specificity which are only for binary classifications)
#Make statistical statements about models performance differences by first collect the resampling results 
resamples <- resamples(list(PLS=pls, LOGITBOOST=logitboost, GLMNET=glmnet,
                            RF=rf, ANN=ann))
#Summarise resampling distributions
summary(resamples)
#Visualize model comparison
splom(resamples)
densityplot(resamples)


#Boxplot of model performance metrics
theme1 <- trellis.par.get()
theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
theme1$plot.symbol$pch = 30
theme1$plot.line$col = rgb(1, 0, 0, .7)
theme1$plot.line$lwd <- 2
trellis.par.set(theme1)
model_comparison_boxplot <- bwplot(resamples, layout=c(3, 1), cex=3, cex.axis=10)



#Confidence intervals for accuracy and kappa metrics
trellis.par.set(caretTheme())
model_comparison_accuracy_plot <- dotplot(resamples, metric="Accuracy")
model_comparison_kappa_plot <- dotplot(resamples, metric="Kappa")


#Calculate statistical differences between models, simple t-test
model_comparison_differences <- diff(resamples)
summary(model_comparison_differences)


#Export high quality figures
#Model comparison - boxplot
Cairo(dev.size("cm")[[1]], dev.size("cm")[[2]],
      file="Model_comparison_boxplot.png",
      type="png", bg="transparent", dpi=600, units="cm")
model_comparison_boxplot
dev.off()


#Model comparison - accuracy plot
Cairo(dev.size("cm")[[1]], dev.size("cm")[[2]],
      file="Model_comparison_accuracy_plot.png",
      type="png", bg="transparent", dpi=600, units="cm")
model_comparison_accuracy_plot
dev.off()


#Model comparison - kappa plot
Cairo(dev.size("cm")[[1]], dev.size("cm")[[2]],
      file="Model_comparison_kappa_plot.png",
      type="png", bg="transparent", dpi=600, units="cm")
model_comparison_kappa_plot
dev.off()




#Generate performance metrics tables
model_metrics <- rbind(pls_confusion_matrix_train$byClass,
                       logitboost_confusion_matrix_train$byClass,
                       glmnet_confusion_matrix_train$byClass,
                       rf_confusion_matrix_train$byClass,
                       ann_confusion_matrix_train$byClass,
                       pls_confusion_matrix_test$byClass,
                       logitboost_confusion_matrix_test$byClass,
                       glmnet_confusion_matrix_test$byClass,
                       rf_confusion_matrix_test$byClass,
                       ann_confusion_matrix_test$byClass)
#Convert the table to a data frame
model_metrics <- data.frame(model_metrics)
#Set row names
rownames(model_metrics) <- c("PLS - Train phase",
                             "LOGITBOOST - Train phase",
                             "GLMNET - Train phase",
                             "RF - Train phase",
                             "ANN - Train phase",
                             "PLS - Test phase",
                             "LOGITBOOST - Test phase",
                             "GLMNET - Test phase",
                             "RF - Test phase",
                             "ANN - Test phase")
#set column names
names(model_metrics) <- paste0(names(model_metrics), " (%)")
names(model_metrics) <- gsub("[.]", " ", names(model_metrics))
#Round performance metrics values
model_metrics <- format(round(model_metrics*100, digits=2), nsmall=2)
#Write xlsx file
xlsx::write.xlsx(model_metrics, "Model_metrics.xlsx")




#Perform a variable importance analysis
#Model-specific variable importance metrics are available for all models except LOGITBOOST
#General variable importance method implemented in caret package will be used to analyse LOGITBOOST model
#PLS variable importance
pls_importance <- abs(data.frame(varImp(pls)$importance))
#Decreasing order of variable importance coefficients
pls_importance <- pls_importance[order(pls_importance$Overall, decreasing=T),, drop=F]
#Transform variable importance coefficients to relative percentages (%)
pls_importance <- round(pls_importance*100/max(pls_importance), digits=2)
#Format row names
rownames(pls_importance) <- gsub("_", " ", rownames(pls_importance))


#LOGITBOOST variable importance
logitboost_importance <- abs(data.frame(varImp(logitboost)$importance))
#In binary classifications, the absolute values of variable importance coefficients are the same for both categories
logitboost_importance <- logitboost_importance[,1, drop=F]
#Decreasing order of variable importance coefficients
logitboost_importance <- logitboost_importance[order(logitboost_importance[,1], decreasing=T),, drop=F]
#Transform variable importance coefficients to relative percentages (%)
logitboost_importance <- round(logitboost_importance*100/max(logitboost_importance), digits=2)
#Format row names
rownames(logitboost_importance) <- gsub("_", " ", rownames(logitboost_importance))


#GLMNET variable importance
glmnet_importance <- abs(data.frame(varImp(glmnet)$importance))
#Decreasing order of variable importance coefficients
glmnet_importance <- glmnet_importance[order(glmnet_importance[,1], decreasing=T),, drop=F]
#Transform variable importance coefficients to relative percentages (%)
glmnet_importance <- round(glmnet_importance*100/max(glmnet_importance), digits=2)
#Format row names
rownames(glmnet_importance) <- gsub("_", " ", rownames(glmnet_importance))


#RF variable importance
rf_importance <- as.data.frame(importance(rf$finalModel))
#Select the global mean decrease accuracy indicator as variable importance estimator
rf_importance <- abs(rf_importance[,3, drop=F])
#Decreasing order of variable importance coefficients
rf_importance <- rf_importance[order(rf_importance[,1], decreasing=T),, drop=F]
#Transform variable importance coefficients to relative percentages (%)
rf_importance <- round(rf_importance*100/max(rf_importance), digits=2)
#Format row names
rownames(rf_importance) <- gsub("_", " ", rownames(rf_importance))


#ANN variable importance
ann_importance_label_1 <- olden(ann$finalModel, bar_plot=F, out_var=paste0("Output_", ann$levels[[1]]))
ann_importance_label_2 <- olden(ann$finalModel, bar_plot=F, out_var=paste0("Output_", ann$levels[[2]]))
ann_importance_labels <- cbind.data.frame(ann_importance_label_1, ann_importance_label_2)
#In binary classifications, the absolute values of variable importance coefficients are the same for both categories
ann_importance <- abs(ann_importance_labels[,1, drop=F])
#Decreasing order of variable importance coefficients
ann_importance <- ann_importance[order(ann_importance[,1], decreasing=T),, drop=F]
#Format row names
rownames(ann_importance) <- gsub("Input_", "", rownames(ann_importance))
rownames(ann_importance) <- gsub("_", " ", rownames(ann_importance))


#Transform variable importance coefficients to relative percentages (%)
ann_importance <- round(ann_importance*100/max(ann_importance), digits=2)


#Create a new table containing variable importance coefficients for all models
model_variable_importance <- cbind.data.frame(rownames(pls_importance), pls_importance,
                                              rownames(logitboost_importance), logitboost_importance,
                                              rownames(glmnet_importance), glmnet_importance,
                                              rownames(rf_importance), rf_importance,
                                              rownames(ann_importance), ann_importance)
#Set row names
rownames(model_variable_importance) <- seq(1, nrow(model_variable_importance))
#Set column names
names(model_variable_importance) <- c("Domain", "PLS: Importance (%)",
                                      "Domain", "LOGITBOOST: Importance (%)",
                                      "Domain", "GLMNET: Importance (%)",
                                      "Domain", "RF: Importance (%)",
                                      "Domain", "ANN: Importance (%)")
#Write xlsx file
xlsx::write.xlsx(model_variable_importance, "Model_variable_importance.xlsx")




#Generate a table summarising the distribution of relevant domains
summary_domains <- data.frame(matrix(nrow=0, ncol=4))
#Iterate domain domains
for(i in 1:length(names(train_variables))) {
  
  #Select i domain domain
  i_domain <- names(train_variables)[i]
  
  #Iterate the levels of factor variable used as label for classification models
  for (j in 1:length(levels(domain_data_numeric$Taxa))) {
    
    #Subset data corresponding to j level
    j_subset <- subset(domain_data_numeric, Taxa == levels(domain_data_numeric$Taxa)[j])
    
    #Subset data corresponding to j level and i domain domain
    i_j_subset <- j_subset[,i_domain, drop=F]
    
    #Get sample names corresponding to j level showing i domain domain
    i_j_genomes <- rownames(i_j_subset)[which(i_j_subset[,1] > 0)]
    
    #Iterate sample names corresponding to j level showing i domain domain
    if (length(i_j_genomes) > 0) {
      
      #Generate a table showing the number of taxa containing i domain domain
      i_j_frequency <- data.frame(table(gsub("_GCA.*", "", i_j_genomes)))
      i_j_frequency <- i_j_frequency[order(i_j_frequency$Freq, decreasing=T),]
      
      #Calculate the number of genomes corresponding to j level showing i domain domain
      i_j_n_genomes <- length(i_j_genomes)
      
      #Combine relevant information on domain distribution
      i_j_taxa <- paste0(i_domain, ": ",
                         paste(i_j_frequency$Var1, 
                               i_j_frequency$Freq, 
                               sep=" n=", collapse=" and "))
      
      #Append the results from each iteration to the summary table previously created
      summary_domains <- rbind(summary_domains, 
                                    cbind(levels(domain_data_numeric$Taxa)[j], 
                                          i_domain, i_j_n_genomes, i_j_taxa))
  
    }
  }
}


#Format the table summarising the distribution of relevant domains
names(summary_domains) <- c("Group", "Domain", "Total genomes", "Description")
summary_domains$Domain <- as.factor(gsub("_", " ", summary_domains$Domain))
summary_domains$Description  <- gsub("_", " ", summary_domains$Description)
summary_domains$Group <- as.factor(summary_domains$Group)
summary_domains$`Total genomes` <- as.numeric(summary_domains$`Total genomes`)




#Separate the data from domain more abundant in different groups in different data frames
#Create an empty list of empty data frames
differences_domains <- vector(mode="list", 
                                   length=length(levels(summary_domains$Group)))
names(differences_domains) <- levels(summary_domains$Group)
# differences_domains[1:length(differences_domains)] <- data.frame()

#Iterate domain domains
for (i in 1:length(levels(summary_domains$Domain))) {
  
  #Subset data corresponding to i domain domain
  i_subset <- subset(summary_domains, Domain == levels(summary_domains$Domain)[i])
  
  #Determine the group showing a higher number of genomes containing i domain 
  i_max_group <- i_subset$Group[which(i_subset$`Total genomes` == max(i_subset$`Total genomes`))]
  
  #Append domain data to the corresponding element list
  i_max_element <- which(names(differences_domains) == i_max_group)
  differences_domains[i_max_element][[1]] <- rbind(differences_domains[i_max_element][[1]],
                                                        i_subset)
  
}


#Format the results corresponding to each group
differences_increased_group <- differences_domains$Increased
differences_reduced_group <- differences_domains$Reduced


#Order data frames according to domain domains
#General summary
summary_domains <- summary_domains[order(factor(summary_domains$Domain, 
                                             levels=unique(mixedsort(summary_domains$Domain)))),]
#domains more abundant in increased group
differences_increased_group <- differences_increased_group[order(factor(differences_increased_group$Domain, 
                                                          levels=unique(mixedsort(differences_increased_group$Domain)))),]
#domains more abundant in reduced group
differences_reduced_group <- differences_reduced_group[order(factor(differences_reduced_group$Domain, 
                                                       levels=unique(mixedsort(differences_reduced_group$Domain)))),]
#Write xlsx files
xlsx::write.xlsx(summary_domains, "Summary_domains.xlsx")
xlsx::write.xlsx(differences_increased_group, "Differences_increased_group.xlsx")
xlsx::write.xlsx(differences_reduced_group, "Differences_reduced_group.xlsx")


#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################


#Finish
