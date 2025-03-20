# This script is designed to download the mRNA, miRNA, and CpG data
# available for the samples that were identified in the 1.Get_data script.
# Make sure to update the script for the data of your interest and, again, 
# be careful with possible updates in the TCGA platform that could change
# how this whole script works

######################################################################
#------------------------------Load libraries-------------------------
######################################################################
library(TCGAbiolinks)          # Version: 2.35.3
library(GEOquery)              # Version: 2.74.0



######################################################################
#------------------------------Load previous objects------------------
######################################################################
# See script 1.Get_data to understand how and which objects were made
samples_data <- read.csv("3_Data/samples_data.tsv", sep = "\t", header = T)



######################################################################
#------------------------Fetch mRNA expression data-------------------
######################################################################
# Check for available mRNA HCC data
exp_data <- GDCquery(project = "TCGA-LIHC",                          # Liver hepatocellular carcinoma
                     data.category = "Transcriptome Profiling",      # Self explanatory 
                     data.type = "Gene Expression Quantification",   # Self explanatory
                     workflow.type="STAR - Counts",                  # Read counts
                     barcode = samples_data$Experiment_id)           # Barcodes used to filter the download files

# Do not run if you've already downloaded your data
GDCdownload(exp_data, directory = "3_Data/")                         # Download data (will create GDCdata folder unless directory is set)



######################################################################
#--------------------Fetch miRNA expression data----------------------
######################################################################
# Check for available mRNA HCC data
mirna_data <- GDCquery(project = "TCGA-LIHC",                        # Liver hepatocellular carcinoma
                       data.category = "Transcriptome Profiling",    # Self explanatory 
                       data.type ="miRNA Expression Quantification", # Self explanatory 
                       barcode = samples_data$Experiment_id)         # Barcodes used to filter the download files

# Do not run if you've already downloaded your data
GDCdownload(mirna_data, directory = "3_Data/")                       



######################################################################
#--------------------Fetch CpG expression data------------------------
######################################################################
# Check for available mRNA HCC data
meth_data <-  GDCquery(project = "TCGA-LIHC",                         # Liver hepatocellular carcinoma
                       data.category = "DNA Methylation",             # Self explanatory
                       platform="Illumina Human Methylation 450",     # Self explanatory
                       data.type = "Methylation Beta Value",          # Methylation data type (this one is already pre-processed)
                       barcode = samples_data$Experiment_id)          # Barcodes used to filter the download files

# Do not run if you've already downloaded your data
GDCdownload(meth_data, directory = "3_Data/", files.per.chunk = 3)    # Download data (will create GDCdata folder unless directory is set)


