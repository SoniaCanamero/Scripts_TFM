
#Load the required R libraries
getwd()
library(knitr)
suppressPackageStartupMessages(library(googleVis))
library(googleVis)
library(plotly)
library(tidyverse)
library(ggplot2)
library(ggdendro)
library(reshape2)
library(grid)
library(stringi)
library(gtools)
library(Cairo)
rm(list=ls())


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


#Create a file list comprising cazy functional domains
file_list <- list.files(path = "cazy_functional_domains/")


#Read individual csv files from the file list and store them in a raw file list
raw_files <- list()
for (i in 1:length(file_list)) {
  raw_files[[i]] <- read.csv(paste0("cazy_functional_domains/", file_list[i]), header=T, sep="\t")
}


#Select high-quality functional domains based on a coverage threshold (i.e. 0.95)
coverage_threshold <- 0.95
#Iterate individual files to select high-quality functional domains
filtered_domains <- list()
for (i in 1:length(raw_files)) {
  filtered_domains[[i]] <- raw_files[[i]][which(raw_files[[i]]$Coverage >= coverage_threshold), ]
}


#The length of lists containing raw csv files and filtered domain data must be the same
length(file_list) == length(filtered_domains)
length(raw_files) == length(filtered_domains)


#Set the names of the list containing filtered domain data
names(filtered_domains) <- sub("(.)", "\\U\\1", gsub("_hmmer.*", "", file_list), perl=T)


#Gather codes indicating glycosidase families from all samples present in the list containing filtered domain data
glycosidase_data <- data.frame()
for (i in 1:length(filtered_domains)) {
  glycosidase_data <- rbind(glycosidase_data, as.data.frame(filtered_domains[[i]]$HMM.Profile))
}


#Remove duplicated glycosidase codes from the glycosidase data frame
glycosidase_data <- unique(glycosidase_data)
names(glycosidase_data) <- c("Glycosidase")


#Generate a new data frame summarising the number of glycosidase families per genome sequence
for (i in 1:length(filtered_domains)) {
  
  #Determine glycosidase families present in i genome sequence
  i_glycosidases <- intersect(glycosidase_data$Glycosidase, filtered_domains[[i]]$HMM.Profile)
  i_glycosidases <- data.frame(i_glycosidases)
  
  #Continue the analysis if any glycosidases are found in i genome sequence
  if (length(i_glycosidases) > 0) {
    names(i_glycosidases) <- names(filtered_domains)[i]
    
    #If the number of glycosidases found in i genome is lower than total number of glycosidases found in the complete dataset, 
    #Then, add NAs corresponding to the missing glycosidases
    if (nrow(i_glycosidases) < nrow(glycosidase_data)) {
      i_add_nas <- nrow(glycosidase_data) - nrow(i_glycosidases)
      i_glycosidases[nrow(i_glycosidases)+i_add_nas,] <- NA
    }
    
    #Align glycosidases found in i genome to match the total number of glycosidases found in the complete dataset
    i_match_index <- match(glycosidase_data$Glycosidase, i_glycosidases[,1])
    i_glycosidases_match  <- i_glycosidases[i_match_index,]
    
    #glycosidases found in i genomes to the complete dataset
    glycosidase_data <- cbind(glycosidase_data, i_glycosidases_match)
    #Set the column name corresponding to i genome
    names(glycosidase_data)[ncol(glycosidase_data)] <- names(i_glycosidases)
    
  }
}


#Format the the data frame summarising the number of glycosidase families per genome sequence
#Set row names (i.e. glycosidase domain) 
rownames(glycosidase_data) <- gsub(".hmm", "", glycosidase_data$Glycosidase)
#Remove column containing glycosidase names 
glycosidase_data <- glycosidase_data[,2:ncol(glycosidase_data)]
#Replace characters in cells by 1
glycosidase_data_numeric <- glycosidase_data
glycosidase_data_numeric <- glycosidase_data_numeric %>% mutate_all(funs(str_replace_all(., ".*hmm", "1")))
#Replace NAs by 0
glycosidase_data_numeric[is.na(glycosidase_data_numeric)] <- 0
#Convert the data frame to numeric
glycosidase_data_numeric <- mutate_all(glycosidase_data_numeric, function(x) as.numeric(as.character(x)))


#Generate a list of accession codes corresponding to each genome included in the dataset
accesion_codes <- as.data.frame(gsub(" GCA ", ": GCA_", gsub("_", " ", names(glycosidase_data_numeric))))
names(accesion_codes) <- c("Taxa: accession code")
#Write xlsx file
xlsx::write.xlsx(accesion_codes, "Accession_codes_table.xlsx")


#Get taxa names and generate a new list
taxa <- unique(gsub("_GCA.*", "", names(glycosidase_data_numeric)))
glycosidase_presence_by_taxa <- list()
#Create a new data frame to store glycosidase percentages by taxa
glycosidase_percentages_by_taxa <- data.frame(matrix(nrow=0, ncol=0))


#Generate a new data frame summarising the number of glycosidase families per taxa
for (i in 1:length(taxa)) {
  
  #Select genomes corresponding to i taxa
  i_taxa <- grep(taxa[i], colnames(glycosidase_data_numeric))
  #Select glycosidases corresponding to genomes from i taxa
  i_glycosidases <- data.frame(glycosidase_data_numeric[ ,i_taxa, drop=F])
  
  #Determine presence/absence of a domain glycosidase considering all genomes from i taxa
  i_presence <- data.frame()
  for (j in 1:nrow(i_glycosidases)) {
    if (any(i_glycosidases[j,] == 1) == TRUE) {
      #Calculate the percentage of genomes from i taxa containing glycosidase domain j
      j_percentage <- (sum(i_glycosidases[j,])*100)/length(i_glycosidases[j,])
      i_presence <- rbind(i_presence, j_percentage)
      #Append j glycosidase percentages by i taxa to a general data frame
      glycosidase_percentages_by_taxa <- rbind.data.frame(glycosidase_percentages_by_taxa,
                                                          paste0(taxa[i], " - ", gsub("_", " ", rownames(i_glycosidases)[j]), 
                                                                 " ", format(round(j_percentage, 2), nsmall =2), 
                                                                 "% (n=", sum(i_glycosidases[j,]), "/", 
                                                                 length(i_glycosidases[j,]), " genomes)"))
    } else {
      i_presence <- rbind(i_presence, 0)
    }
  }
  
  #Append he percentage of genomes from i taxa containing different glycosidases to the list
  glycosidase_presence_by_taxa[[i]] <- i_presence
  names(glycosidase_presence_by_taxa)[i] <- paste0(taxa[i], " (n=", ncol(i_glycosidases), " genomes)")
  
}


#Convert the list showing glycosidase presence by taxa to data frame
glycosidase_data_by_taxa <- data.frame(matrix(nrow=nrow(glycosidase_data_numeric), ncol=0))
for (i in 1:length(glycosidase_presence_by_taxa)) {
  glycosidase_data_by_taxa <- cbind(glycosidase_data_by_taxa, glycosidase_presence_by_taxa[[i]])
}


#Format the data frame showing glycosidase presence by taxa
#Set column names
names(glycosidase_data_by_taxa) <- names(glycosidase_presence_by_taxa)
#Set row names
rownames(glycosidase_data_by_taxa) <- gsub("_", " ", rownames(glycosidase_data_numeric))


#Note: glycosyl transferase and CBM35inCE17 might be annotated with different names resulting in duplicated annotations
# glycosidase_data_by_taxa <- glycosidase_data_by_taxa[-grep("Glyco trans ", rownames(glycosidase_data_by_taxa)),]
glycosidase_data_by_taxa <- glycosidase_data_by_taxa[-grep("CBM35inCE17", rownames(glycosidase_data_by_taxa)),]
#Format "Glyco tranf" domain - standarised enzyme codes
rownames(glycosidase_data_by_taxa) <- gsub("Glyco tranf ", "", rownames(glycosidase_data_by_taxa))
rownames(glycosidase_data_by_taxa) <- gsub("Glyco trans ", "", rownames(glycosidase_data_by_taxa))
rownames(glycosidase_data_by_taxa) <- gsub("Glycos transf ", "", rownames(glycosidase_data_by_taxa))


#Order the data frame showing glycosidase presence by taxa
#Order rows by activities
glycosidase_data_by_taxa <- glycosidase_data_by_taxa[order(factor(rownames(glycosidase_data_by_taxa), 
                                                                   levels=unique(mixedsort(rownames(glycosidase_data_by_taxa))))),]
#Order columns by microbial species
glycosidase_data_by_taxa <- glycosidase_data_by_taxa[,order(factor(colnames(glycosidase_data_by_taxa), 
                                                       levels=unique(mixedsort(colnames(glycosidase_data_by_taxa)))))]


#Format the data frame showing glycosidase percentages by taxa
names(glycosidase_percentages_by_taxa) <- c("Glycosidase percentages (%) by taxa") 
#Write xlsx file
xlsx::write.xlsx(glycosidase_percentages_by_taxa, "Glycosidase_percentages_by_taxa.xlsx")


#Generate a heatmap plot with ggplot2
#A high number of glycosidase functional domains were annotated
#Three different heatmaps showing these domains will be plotted for visualization purposes
#Select a cutoff to subset data (i.e. one third of total glycosidase domains)
cutoff <- round(nrow(glycosidase_data_by_taxa)/3,1)


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


#Subset glycosidase domain data to plot heatmap 1
subset <- glycosidase_data_by_taxa[1:cutoff,]


#Generate dendrogram to include in the heatmap
dendrogram <- as.dendrogram(hclust(d=dist(x=t(subset))))
dendrogram_plot <- ggdendrogram(data=dendrogram, rotate=T) +  
                     theme(axis.text.y=element_text(size=10.5, face="italic"))


#Format data for heatmap plotting
heatmap_data <- melt(t(subset), id=colnames(t(subset)))
#Re-order heatmap rows to match dendrogram branches
dendrogram_order <- order.dendrogram(dendrogram)
#Order the levels according to their position in the cluster
heatmap_data$Var1 <- factor(x=heatmap_data$Var1,
                              levels=heatmap_data$Var1[dendrogram_order], 
                              ordered=T)


#Create the heatmap plot
heatmap_plot <- ggplot(data=heatmap_data, aes(x=Var2, y=Var1, fill=value)) +
                  geom_tile(aes(fill=value)) +
                  scale_fill_gradient2(name="Percentage (%) of genomes containing each domain", 
                                       low="white", high="black") + labs(x="") +
                  scale_y_discrete(position="right") + 
                  geom_tile(colour="gray50") +
                  theme(text=element_text(size=8), legend.text=element_text(size=12), legend.title=element_text(size=12),
                        axis.text.y=element_blank(),
                        axis.title.y=element_blank(),
                        legend.position="top", 
                        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, face="bold"))


#Print all together
x11()
grid.newpage()
print(heatmap_plot, vp=viewport(x=0.37, y=0.5, width=0.75, height=1.0))
print(dendrogram_plot, vp=viewport(x=0.87, y=0.46, width=0.25, height=0.88))


#Export high quality figures
recorded_heatmap <- recordPlot()
Cairo(dev.size("cm")[[1]], dev.size("cm")[[2]],
             file="Glycosidase_heatmap_1.png",
             type="png", bg="transparent", dpi=600, units="cm")
recorded_heatmap
dev.off()


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


#Subset glycosidase domain data to plot heatmap 2
subset <- glycosidase_data_by_taxa[(cutoff+1):(cutoff*2),]


#Generate dendrogram to include in the heatmap
dendrogram <- as.dendrogram(hclust(d=dist(x=t(subset))))
dendrogram_plot <- ggdendrogram(data=dendrogram, rotate=T) +  
  theme(axis.text.y=element_text(size=10.5, face="italic"))


#Format data for heatmap plotting
heatmap_data <- melt(t(subset), id=colnames(t(subset)))
#Re-order heatmap rows to match dendrogram branches
dendrogram_order <- order.dendrogram(dendrogram)
#Order the levels according to their position in the cluster
heatmap_data$Var1 <- factor(x=heatmap_data$Var1,
                            levels=heatmap_data$Var1[dendrogram_order], 
                            ordered=T)


#Create the heatmap plot
heatmap_plot <- ggplot(data=heatmap_data, aes(x=Var2, y=Var1, fill=value)) +
                  geom_tile(aes(fill=value)) +
                  scale_fill_gradient2(name="Percentage (%) of genomes containing each domain", 
                                       low="white", high="black") + labs(x="") +
                  scale_y_discrete(position="right") + 
                  geom_tile(colour="gray50") +
                  theme(text=element_text(size=8), legend.text=element_text(size=12), legend.title=element_text(size=12),
                        axis.text.y=element_blank(),
                        axis.title.y=element_blank(),
                        legend.position="top", 
                        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, face="bold"))


#Print all together
x11()
grid.newpage()
print(heatmap_plot, vp=viewport(x=0.37, y=0.5, width=0.75, height=1.0))
print(dendrogram_plot, vp=viewport(x=0.87, y=0.46, width=0.25, height=0.88))


#Export high quality figures
recorded_heatmap <- recordPlot()
Cairo(dev.size("cm")[[1]], dev.size("cm")[[2]],
             file="Glycosidase_heatmap_2.png",
             type="png", bg="transparent", dpi=600, units="cm")
recorded_heatmap
dev.off()


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


#Subset glycosidase domain data to plot heatmap 3
subset <- glycosidase_data_by_taxa[(cutoff*2+1):(cutoff*3),]


#Generate dendrogram to include in the heatmap
dendrogram <- as.dendrogram(hclust(d=dist(x=t(subset))))
dendrogram_plot <- ggdendrogram(data=dendrogram, rotate=T) +  
  theme(axis.text.y=element_text(size=10.5, face="italic"))


#Format data for heatmap plotting
heatmap_data <- melt(t(subset), id=colnames(t(subset)))
#Re-order heatmap rows to match dendrogram branches
dendrogram_order <- order.dendrogram(dendrogram)
#Order the levels according to their position in the cluster
heatmap_data$Var1 <- factor(x=heatmap_data$Var1,
                            levels=heatmap_data$Var1[dendrogram_order], 
                            ordered=T)


#Create the heatmap plot
heatmap_plot <- ggplot(data=heatmap_data, aes(x=Var2, y=Var1, fill=value)) +
  geom_tile(aes(fill=value)) +
  scale_fill_gradient2(name="Percentage (%) of genomes containing each domain", 
                       low="white", high="black") + labs(x="") +
  scale_y_discrete(position="right") + 
  geom_tile(colour="gray50") +
  theme(text=element_text(size=8), legend.text=element_text(size=12), legend.title=element_text(size=12),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        legend.position="top", 
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, face="bold"))


#Print all together
x11()
grid.newpage()
print(heatmap_plot, vp=viewport(x=0.37, y=0.5, width=0.75, height=1.0))
print(dendrogram_plot, vp=viewport(x=0.87, y=0.455, width=0.25, height=0.89))


#Export high quality figures
recorded_heatmap <- recordPlot()
Cairo(dev.size("cm")[[1]], dev.size("cm")[[2]],
             file="Glycosidase_heatmap_3.png",
             type="png", bg="transparent", dpi=600, units="cm")
recorded_heatmap
dev.off()


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


#Export glycosidase data frames as R objects (these objects will be used as inputs for machine learning models)
saveRDS(glycosidase_data_by_taxa, "glycosidase_data_by_taxa.rds")
saveRDS(glycosidase_data_numeric, "glycosidase_data_numeric.rds")


#Read R objects containing the glycosidase data frames
glycosidase_data_by_taxa <- readRDS("glycosidase_data_by_taxa.rds")
glycosidase_data_numeric <- readRDS("glycosidase_data_numeric.rds")


#Finish
