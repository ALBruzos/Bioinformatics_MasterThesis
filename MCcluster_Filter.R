#!/usr/bin/env

###############################################################################################################
#                                   Script MasterCopyFilter                                                   #
#                                        version 2.0                                                          #
###############################################################################################################

#Read the data from a text file
MC <- read.table("~/Scripts/201605_2_MasterCopy_Clusters/output/somatic2_2_Clustered.txt", sep="\t", header=TRUE)

#Creating an empty dataframe that will be the output of this script
cluster.df <- data.frame(Cluster=character(), Chr.of.MC=character(), Higher.Coord1.MC=numeric(), Lower.Coord2.MC=numeric(), Coord1.Tr=numeric(), Coord2.Tr=numeric(),
                         No.of.Samples=numeric(), No.of.Samples.No.Repeated=numeric(), Samples=character(), stringsAsFactors=FALSE) 

#Creating a vector with the names of the different clusters without the NAs or 0s.
CLnames <- unique(MC$Clusters)
#Deleting the NA
CLnames <- CLnames[!is.na(CLnames)]
#Deleting the 0
end <- length(CLnames); CLnames <- c(as.character(CLnames[1]), as.character(CLnames[3:end]))

#Loop to iterate among the different clusters
for (cl in CLnames){
    #Subsetting as an small dataframe called cluster.subset the rows belonging to the cl cluster
    cluster.subset <- subset(MC, MC$Clusters==cl)  
    
    #Select the Tissue Type and Tumour Code of all the samples of that cluster; write them in a vector as one element.
    samples <- paste(cluster.subset$Tissue.Type, "#", cluster.subset$Tumour.Code, sep = "")
    #Compute the number of samples of that cluster
    total.samples.number <- length(samples)
    #Select only no repeated samples of that cluster
    samples.no.repeated <- unique(samples)
    #and store the number of samples no repeated of that clusterq
    no.repeated.samples.number <- length(samples.no.repeated)
    #Collapse all the elements of the vector with the samples no repeated in one element
    samples.no.repeated <- paste(samples.no.repeated[1:length(samples.no.repeated)], collapse=',')
   
    #Creating a dataframe (info) with the values that we want to have in the dataframe for each cluster. 
    info <- data.frame(as.character(cluster.subset$Clusters[1]), as.character(cluster.subset$Chromosome.of.MC[1]), as.numeric(max(cluster.subset$Coordinate.1.of.MC)), 
              as.numeric(min(cluster.subset$Coordinate.2.of.MC)), as.numeric(cluster.subset$Coordinate.1.of.transduced.region[1]), 
              as.integer(as.character(cluster.subset$Coordinate.2.of.transduced.region[1])), as.numeric(total.samples.number), as.numeric(no.repeated.samples.number),
              as.character(samples.no.repeated))
    #Binding the info dataframe to cluster.df output dataframe
    cluster.df <- rbind(cluster.df, info)
    
}

#Naming the columns of the output cluster.df dataframe
colnames(cluster.df) <- c("Cluster", "Choromosome of MC", "HIGHER Coordinate 1 of MC", "LOWER Coordinate 2 of MC", "1st row of Coordinate 1 of transduced region",
                          "1st row of Coordinate 2 of transduced region","Total number of samples", "Number of samples no repeated", "Samples (TissueType#TumourCode)") 

#Exporting cluster.df object into a tab delimited text file
write.table(cluster.df, "~/Scripts/201605_2_MasterCopy_Clusters/output/somatic2_4_ClustersInfo.txt", row.names = FALSE, sep="\t")

###########
# THE END #
###########
