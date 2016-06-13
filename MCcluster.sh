#!/bin/bash

#######################################
##
##  MasterCopyCluster Script
##	author: Alicia L. Bruzos
##	date: May 10, 2016
##	version: 1.0
##
#######################################


# Step 1: Filtering somatic2.txt with bash commands
###################################################
printf "\n\e[1;30m\e[34mSTEP 1\n  \e[0;37mStarting to filter somatic2.txt with bash commands\n  ...\e[0m\n"
mkdir output
awk 'OFS="\t" {if(($9 == "td1") || ($9 == "td2")){print$10,$11,$12,$13,$14,$15,$4,$5,$6,$16,$17,$1,$3}}' somatic2.txt | sort -k1,1n -k2,2n | grep putative | awk '$1!="Y"' > output/somatic2_1_Filtered.txt
printf "  \e[0;37mEnd of filtering somatic2.txt with bash commands \e[0m\nOutput file created: \e[1;30msomatic2_1_Filtered.txt\n\n"

# Step 2: Creating clusters in a new column of somatic2.txt using "MCFinder.R" script
#####################################################################################
printf "\n\e[34mSTEP 2\n  \e[0;37mCreating clusters in a new column of somatic2.txt\n  ...\e[0m\n"
Rscript src/MCcluster_Finder.R
printf "\n  \e[0;37m...\n  Clusters created as a new column of somatic2.txt\e[0m\nTwo output files created: \e[1;30msomatic2_2_Clustered.txt\e[0m  &  \e[1;30msomatic2_3_NoClusterSamples.txt\n\n"

# Step 3: Create output file with clusters information using "MCFilter.R" script
################################################################################
printf "\n\e[34mSTEP 3\n  \e[0;37mCreating output file with clusters information\n  ...\n\e[0m"
Rscript src/MCcluster_Filter.R 
printf "\n  \e[0;37m...\n  End of creating a file with clusters information\e[0m\nOutput file created: \e[1;30msomatic2_4_ClustersInfo.txt\n\n\e[0m\e[5mTake a look to the \e[4moutput file\e[0m\n\n\e[1;30m\e[34mTHE END\n\n"

###########
# THE END #
###########