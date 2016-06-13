#!/usr/bin/env

###############################################################################################################
#                                Script MasterCopyFinder                                                      #
#                                        version 5                                                            #
###############################################################################################################

#VARIABLE TO BE CHANGED:
threshold=500

#Read file as a data.frame (separated by tabs)
MC <- read.table("~/Scripts/201605_2_MasterCopy_Clusters/output/somatic2_1_Filtered.txt", sep="\t", header=FALSE)

######################################################################
###   Function to compare the lines                                ###
######################################################################
MCFinder <- function(chr){
  #Selecting data for just one chromosome
  MC.chr <- data.frame(MC[MC$V1==chr,])
  #Calculating the dimensions of the dataframe
  totalLen <- dim(MC.chr)[1]
  #Create a new data.frame to compare by two consecutives lines
  MC.comp <- data.frame(id1=MC.chr$V10[1:totalLen-1],val1=MC.chr$V3[1:totalLen-1],id2=MC.chr$V10[2:totalLen],val2=MC.chr$V3[2:totalLen])
  #Compare the values two by two and naming the ones that correspond to a cluster
  MC.comp$diff <- abs(MC.comp$val2-MC.comp$val1)
  #Creating a new column where the clusters would be written
  MC.comp[MC.comp$diff>threshold, "cluster"] <- 0
  MC.comp[MC.comp$diff<=threshold, "cluster"] <- paste("Cluster.chr", chr, sep="")

  # Get all the rows that have a cluster assigned
  clusterIds<-which(MC.comp$cluster!=0)
  # initialize the cluster id counter
  clusterIdCounter<-1
  # Iterate over all the rows that belong to a cluster. Initialize first element, which belongs to the first cluster (counter default value).
  # Then, the loop check current row vs. current row +1, this relative to the MC table
  MC.comp[MC.comp$cluster!=0, "cluster"][1]<-paste(MC.comp[MC.comp$cluster!=0, "cluster"][1],".",clusterIdCounter, sep="")
  if (length(clusterIds) > 1){
    for(item in 1:(length(clusterIds)-1)){
      # check the next one
      if((clusterIds[item]+1)!=clusterIds[item+1]){
        # create new cluster
        clusterIdCounter<-clusterIdCounter+1
      }
      # Assign value
      MC.comp[MC.comp$cluster!=0, "cluster"][item+1]<-paste(MC.comp[MC.comp$cluster!=0, "cluster"][item+1],".",clusterIdCounter, sep="")
    }
  }

  #Creating a vector with the clusters results
  clusters <- as.vector(MC.comp$cluster)
  #Adding last value of the comparisons to the vector
  clusters <- c(clusters, clusters[length(clusters)])
  i=1
  while(i < (length(clusters)-1)){
    if((clusters[i+1]==0) && (clusters[i]!=clusters[i+1])){
      clusters[i+1]<-clusters[i]
      i=i+1
    }
  i=i+1
  }
  return(clusters)
}

######################################################################
###   Function to do the MCFinder function for every chromosome    ###
######################################################################
MCFinder.main <- function(){
  # Vector with all the chromosomes available in the input folder (all the elements are characters)
  chrNames<-c('MT','X',1:22)
  # Vector to save the results of the doing MCFinder for each row
  MC.output<<-c()
  #Loop to iterate among the different chromosomes
  for (chr in chrNames){
    #print(chr)
    if (is.na(as.numeric(chr))){ #For the cases of MT and chr X 
      MC.output<<-c(MC.output,MCFinder(chr))
    }else{ #For chromosomes from 1 until 22
      MC.output<<-c(MC.output,MCFinder(as.numeric(chr))) #It will use the name of the chrs as numeric
    }
  }
  #Concatenate the vector with cluster information to the input data creating a new column at the end
  MC<<-cbind(MC,MC.output)
}

#Do the task
MCFinder.main()
#Naming the columns
colnames(MC) <- c("Chromosome of MC", "Coordinate 1 of MC", "Coordinate 2 of MC", "MC orientation", "Coordinate 1 of transduced region",
  "Coordinate 2 of transduced region","Chromosome of sample", "Coordinate 1 of sample", "Coordinate 2 of sample","Transduced Material Distance", 
  "Transduced Material Size", "Tissue Type", "Tumour Code", "Clusters")

#Exporting MC object into a tab delimited text file
write.table(MC, "~/Scripts/201605_2_MasterCopy_Clusters/output/somatic2_2_Clustered.txt", row.names = FALSE, sep="\t", eol="\n")

#Subsetting the samples which do not belong to any cluster
MC.noCluster <- subset(MC, MC$Clusters==0)
#Exporting MC.noCluster object into a tab delimited text file
write.table(MC.noCluster, "~/Scripts/201605_2_MasterCopy_Clusters/output/somatic2_3_NoClusterSamples.txt", row.names = FALSE, sep="\t", eol="\n")

###########
# THE END #
###########