#This is the main script to call
#Configure it to fit your file names

################PACKAGES AND OPTIONS################
options(stringsAsFactors = F)
library("reshape2")
library("stringr")
library("ape")
library("readr")
library("dplyr")
library("lubridate")
library("igraph")
################################

################GLOBAL VARIABLES################
## 'barcodefile' is a tab-separated file with 
# the 1st column corresponding to barcode IDs
# and the remaining columns corresponding to individual SNP calls.
# Barcode IDs should not be one of REF", "ALT", "ALT2" or "ALT3
# Column names 2 and above should be be formatted as
# 'chromosome.position'.
barcodefile <- "" #path/to/barcode file

## 'metadatafile' is a tab-separated file with 
# the 1st column corresponding to barcode IDs
# and the remaining columns corresponding to 
# individual ID, date (YYYY-MM-DD) and location.
metadatafile <- "" #path/to/metadata file

# Set to NA if you don't want to specify an order
# This order will be used to transform chromosome names into
# integers. If you don't mention a chromosome here that is in your
# dataset, it will not be included in the pairwise comparisons
chrOrder <- c("Pf3D7_01_v3", "Pf3D7_02_v3", "Pf3D7_03_v3", 
              "Pf3D7_04_v3", "Pf3D7_05_v3", "Pf3D7_06_v3", 
              "Pf3D7_07_v3", "Pf3D7_08_v3", "Pf3D7_09_v3", 
              "Pf3D7_10_v3", "Pf3D7_11_v3", "Pf3D7_12_v3",
              "Pf3D7_13_v3", "Pf3D7_14_v3")

# Directory to output all data
outDir <- "test"

hmmIBD <- "/opt/hmmIBD/hmmIBD" # path/to/hmmIBD binary
runInBackground <- TRUE #Run hmmIBD in a background process

minsites <- 30 #Minimal number of comparable sites to consider a pair
clusterIBDmin <- 0.5 #Value of IBD above which barcodes are considered similar
################################

################STEPS################
## Formatting barcodes
source("src/barcodeRead.R")
source("src/getCompSites.R")
## 

## Pairwise similarity
source("src/pairwiseSim.R")
##

## hmmIBD
source("src/formatBarcodesForHmmIBD.R")
# hmmIBD could take a long time to run, 
#you could send hmmIBD command to a server to speed up the process.
# To calculate all IBD values of all pairs, you can change 'min_inform = 10' to
# 'min_inform = 0' in hmmIBD source code
system2("./src/run-hmmibd-local.sh", args = c(hmmIBD, outDir, runInBackground), 
        wait=TRUE)
##

##Format network
source("src/makeIBDNetwork.R")
##
################################




