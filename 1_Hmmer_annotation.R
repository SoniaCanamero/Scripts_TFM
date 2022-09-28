
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


#Create a file list comprising hmmer functional domains
file_list <- list.files(path = "tbl/")


#Read individual csv files from the file list and store them in a raw file list
raw_files <- list()
for (i in 1:length(file_list)) {
  raw_files[[i]] <- read.table(paste0("tbl/", file_list[i]), header=F, skip=3)
}


#Export raw data as an R object
saveRDS(raw_files, "raw_files.rds")
#Read R object containing the raw data
raw_files <- readRDS("raw_files.rds")


#Set the names of the list containing raw data
names(raw_files) <- sub("(.)", "\\U\\1", gsub("_aaORFs.*", "", file_list), perl=T)


#Gather codes indicating glycosidase families from all samples present in the list containing filtered domain data
domain_data <- data.frame()
for (i in 1:length(raw_files)) {
  domain_data <- rbind(domain_data, as.data.frame(raw_files[[i]]$V4))
}


#Remove duplicated domain codes from the domain data frame
domain_data <- unique(domain_data)
names(domain_data) <- c("Domain")


#Export domain data as an R object
saveRDS(domain_data, "domain_data.rds")
#Read R object containing domain data
domain_data <- readRDS("domain_data.rds")


#Search for neurotransmitter-related functional domains using manually-curated pfam codes
pfam_neurotransmitter_codes <- c("PF00209", "PF02931", "PF02932",
                                 "PF03222", "PF03491", "PF04622",
                                 "PF04664", "PF04680", "PF06387",
                                 "PF06540", "PF08035", "PF10208",
                                 "PF18455")
  

#Prune domain data to select neurotransmitter-related functional domains
domain_data <- domain_data[grep(paste0(pfam_neurotransmitter_codes, 
                                  collapse="|"),
                                  domain_data$Domain),, drop=F]


#Generate a new data frame summarising the number of glycosidase families per genome sequence
for (i in 1:length(raw_files)) {
  
  #Determine glycosidase families present in i genome sequence
  i_domains <- intersect(domain_data$Domain, raw_files[[i]]$V4)
  i_domains <- data.frame(i_domains)
  
  #Continue the analysis if any glycosidases are found in i genome sequence
  if (length(i_domains) > 0) {
    names(i_domains) <- names(raw_files)[i]
    
    #If the number of domains found in i genome is lower than total number of domains found in the complete dataset, 
    #Then, add NAs corresponding to the missing domains
    if (nrow(i_domains) < nrow(domain_data)) {
      i_add_nas <- nrow(domain_data) - nrow(i_domains)
      i_domains[nrow(i_domains)+i_add_nas,] <- NA
    }
    
    #Align domains found in i genome to match the total number of domains found in the complete dataset
    i_match_index <- match(domain_data$Domain, i_domains[,1])
    i_domains_match  <- i_domains[i_match_index,]
    
    #glycosidases found in i genomes to the complete dataset
    domain_data <- cbind(domain_data, i_domains_match)
    #Set the column name corresponding to i genome
    names(domain_data)[ncol(domain_data)] <- names(i_domains)
    
  }
}


#Format the the data frame summarising the number of domains per genome sequence
#Set row names (i.e. functional domain) 
rownames(domain_data) <- gsub("[.].*", "", domain_data$Domain)
#Remove column containing domain names 
domain_data <- domain_data[,2:ncol(domain_data)]
#Replace characters in cells by 1
domain_data_numeric <- domain_data
domain_data_numeric <- domain_data_numeric %>% mutate_all(funs(str_replace_all(., "PF.*", "1")))
#Replace NAs by 0
domain_data_numeric[is.na(domain_data_numeric)] <- 0
#Convert the data frame to numeric
domain_data_numeric <- mutate_all(domain_data_numeric, function(x) as.numeric(as.character(x)))


#Generate a list of accession codes corresponding to each genome included in the dataset
accesion_codes <- as.data.frame(gsub(" GCA ", ": GCA_", gsub("_", " ", names(domain_data_numeric))))
names(accesion_codes) <- c("Taxa: accession code")
#Write xlsx file
xlsx::write.xlsx(accesion_codes, "Accession_codes_table.xlsx")


#Get taxa names and generate a new list
taxa <- unique(gsub("_GCA.*", "", names(domain_data_numeric)))
domain_presence_by_taxa <- list()
#Create a new data frame to store domain percentages by taxa
domain_percentages_by_taxa <- data.frame(matrix(nrow=0, ncol=0))


#Generate a new data frame summarising the number of functional domains per taxa
for (i in 1:length(taxa)) {
  
  #Select genomes corresponding to i taxa
  i_taxa <- grep(taxa[i], colnames(domain_data_numeric))
  #Select domains corresponding to genomes from i taxa
  i_domains <- data.frame(domain_data_numeric[ ,i_taxa, drop=F])
  
  #Determine presence/absence of a functional domain considering all genomes from i taxa
  i_presence <- data.frame()
  for (j in 1:nrow(i_domains)) {
    if (any(i_domains[j,] == 1) == TRUE) {
      #Calculate the percentage of genomes from i taxa containing functional domain j
      j_percentage <- (sum(i_domains[j,])*100)/length(i_domains[j,])
      i_presence <- rbind(i_presence, j_percentage)
      #Append j domain percentages by i taxa to a general data frame
      domain_percentages_by_taxa <- rbind.data.frame(domain_percentages_by_taxa,
                                                          paste0(taxa[i], " - ", gsub("_", " ", rownames(i_domains)[j]), 
                                                                 " ", format(round(j_percentage, 2), nsmall =2), 
                                                                 "% (n=", sum(i_domains[j,]), "/", 
                                                                 length(i_domains[j,]), " genomes)"))
    } else {
      i_presence <- rbind(i_presence, 0)
    }
  }
  
  #Append he percentage of genomes from i taxa containing different glycosidases to the list
  domain_presence_by_taxa[[i]] <- i_presence
  names(domain_presence_by_taxa)[i] <- paste0(taxa[i], " (n=", ncol(i_domains), " genomes)")
  
}


#Convert the list showing glycosidase presence by taxa to data frame
domain_data_by_taxa <- data.frame(matrix(nrow=nrow(domain_data_numeric), ncol=0))
for (i in 1:length(domain_presence_by_taxa)) {
  domain_data_by_taxa <- cbind(domain_data_by_taxa, domain_presence_by_taxa[[i]])
}


#Format the data frame showing domain presence by taxa
#Set column names
names(domain_data_by_taxa) <- names(domain_presence_by_taxa)
#Set row names
rownames(domain_data_by_taxa) <- gsub("_", " ", rownames(domain_data_numeric))


#Order the data frame showing glycosidase presence by taxa
#Order rows by activities
domain_data_by_taxa <- domain_data_by_taxa[order(factor(rownames(domain_data_by_taxa), 
                                                                   levels=unique(mixedsort(rownames(domain_data_by_taxa))))),]
#Order columns by microbial species
domain_data_by_taxa <- domain_data_by_taxa[,order(factor(colnames(domain_data_by_taxa), 
                                                       levels=unique(mixedsort(colnames(domain_data_by_taxa)))))]


#Format the data frame showing glycosidase percentages by taxa
names(domain_percentages_by_taxa) <- c("Domain percentages (%) by taxa") 
#Write xlsx file
xlsx::write.xlsx(domain_percentages_by_taxa, "Domain_percentages_by_taxa.xlsx")


#Generate a heatmap plot with ggplot2


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


#Select functional domain data to plot a heatmap
subset <- domain_data_by_taxa


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
                  theme(text=element_text(size=12), legend.text=element_text(size=12), legend.title=element_text(size=12),
                        axis.text.y=element_blank(),
                        axis.title.y=element_blank(),
                        legend.position="top", 
                        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, face="bold"))


#Print all together
x11()
grid.newpage()
print(heatmap_plot, vp=viewport(x=0.38, y=0.5, width=0.75, height=1.0))
print(dendrogram_plot, vp=viewport(x=0.87, y=0.475, width=0.27, height=0.84))


#Export high quality figures
recorded_heatmap <- recordPlot()
Cairo(dev.size("cm")[[1]], dev.size("cm")[[2]],
             file="Domain_heatmap.png",
             type="png", bg="transparent", dpi=600, units="cm")
recorded_heatmap
dev.off()




##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


#Export glycosidase data frames as R objects (these objects will be used as inputs for machine learning models)
saveRDS(domain_data_by_taxa, "domain_data_by_taxa.rds")
saveRDS(domain_data_numeric, "domain_data_numeric.rds")


#Read R objects containing the glycosidase data frames
domain_data_by_taxa <- readRDS("domain_data_by_taxa.rds")
domain_data_numeric <- readRDS("domain_data_numeric.rds")


#Finish
