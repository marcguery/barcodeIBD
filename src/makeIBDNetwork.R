###########################HMM IBD FILE###########################
databarcodes <- read_delim(file = paste0(outDir, "/barcodes-IBD.hmm_fract.txt.zip"), 
                           delim = "\t")
######################################################

###########################ATTRIBUTING IDs###########################
prevcolnames <- colnames(databarcodes)

databarcodes$Var1 <- pmin(databarcodes$sample1, databarcodes$sample2)
databarcodes$Var2 <- pmax(databarcodes$sample1, databarcodes$sample2)
databarcodes$sample1 <- databarcodes$Var1
databarcodes$sample2 <- databarcodes$Var2
databarcodes <- databarcodes[,c(1:(ncol(databarcodes)-2))]

all(table(c(databarcodes$sample1, databarcodes$sample2))==(nrow(metadata)-1))

databarcodes <- merge(databarcodes,
                      barcodes.comp.num, all = T,
                      by.x = c("sample1", "sample2"),
                      by.y = c("ID1", "ID2"))
databarcodes <- databarcodes[,order(factor(colnames(databarcodes), levels = c(prevcolnames, "N_comp_sites")))]

databarcodes$fract_sites_IBD[is.na(databarcodes$fract_sites_IBD)] <- -1
databarcodes$fract_sites_IBD[databarcodes$N_comp_sites<minsites] <- -1

databarcodes <- databarcodes[!databarcodes$sample1%in%c("REF", "ALT", "ALT2", "ALT3") & 
               !databarcodes$sample2%in%c("REF", "ALT", "ALT2", "ALT3"),]
length(which(databarcodes$fract_sites_IBD==-1))/nrow(databarcodes)
summary(databarcodes$fract_sites_IBD[databarcodes$fract_sites_IBD>-1])
######################################################

###########################SPACE/TIME INFORMATION###########################

##DATES
databarcodes <- merge(databarcodes, metadata,
                      by.x = "sample1", by.y = "ID",
                      all.x = T)
databarcodes <- merge(databarcodes, metadata,
                      by.x = "sample2", by.y = "ID",
                      all.x = T, suffixes = c("1", "2"))

databarcodes$date1 <- as.Date(databarcodes$date1)
databarcodes$date2 <- as.Date(databarcodes$date2)
databarcodes$daysElapsed <- round(abs(difftime(databarcodes$date1,databarcodes$date2, units="days")), digits = 2)

##LOCATIONS
databarcodes$sameIndividual <- databarcodes$individual1==databarcodes$individual2
databarcodes$commonIndividual <- NA
databarcodes$commonIndividual[databarcodes$sameIndividual] <- databarcodes$individual1[databarcodes$sameIndividual]

databarcodes$sameLocation <- databarcodes$location1==databarcodes$location2
databarcodes$commonLocation <- NA
databarcodes$commonLocation[databarcodes$sameLocation] <- databarcodes$shortlocation1[databarcodes$sameLocation]

combos <- strsplit(paste0(databarcodes$shortlocation1[!databarcodes$sameLocation], databarcodes$shortlocation2[!databarcodes$sameLocation]), split="")
combos <- sapply(combos, sort)
combos <- paste0(combos[1,], combos[2,])
databarcodes$commonLocation[!databarcodes$sameLocation] <- combos

##INDIVIDUALS
combos <- strsplit(paste(databarcodes$individual1[!databarcodes$sameIndividual], databarcodes$individual2[!databarcodes$sameIndividual], sep =";"), split=";")
combos <- sapply(combos, sort)
combos <- paste0(combos[1,], combos[2,])
databarcodes$commonIndividual[!databarcodes$sameIndividual] <- combos

#Number of samples
length(unique(c(databarcodes$sample1, databarcodes$sample2)))
#Number of individuals
length(unique(c(databarcodes$individual1, databarcodes$individual2)))
######################################################

#############COMPARABLE SAMPLES############
sample1comparison <- databarcodes %>%
  group_by(sample1)%>%
  summarise(length(which(fract_sites_IBD!=-1)))
colnames(sample1comparison) <- c("sample", "comparisons1")
sample2comparison <- databarcodes %>%
  group_by(sample2)%>%
  summarise(length(which(fract_sites_IBD!=-1)))
colnames(sample2comparison) <- c("sample", "comparisons2")

samplecomparison <- merge(sample1comparison, sample2comparison, all = T)
samplecomparison$comparisons1[is.na(samplecomparison$comparisons1)] <- 0
samplecomparison$comparisons2[is.na(samplecomparison$comparisons2)] <- 0

metadata$comparable <- TRUE
metadata$comparable[metadata$ID%in%samplecomparison$sample[samplecomparison$comparisons1+samplecomparison$comparisons2==0]] <- FALSE
##########################

#############CLUSTERING############
# Get rid of loops and ensure right naming of vertices
databarcodes.topology <- databarcodes[databarcodes$fract_sites_IBD>=clusterIBDmin,c("sample1", "sample2")]
databarcodes.topology <- simplify(graph.data.frame(databarcodes.topology[order(databarcodes.topology[[1]]),],directed = FALSE))

# Find all components
comps <- components(databarcodes.topology)
cluster.df <- data.frame(comps$membership)
cluster.df <- cbind(row.names(cluster.df), cluster.df)
colnames(cluster.df) <- c("ID", "cluster")
metadata <- merge(metadata, cluster.df, all.x = T)
metadata$cluster[is.na(metadata$cluster)] <- 0

metadata$cluster <- sprintf(paste0("%0",max(nchar(metadata$cluster)),"d"), metadata$cluster)

databarcodes <- merge(databarcodes, metadata[,c("ID", "cluster")],
               by.x = "sample1", by.y = "ID", all.x = TRUE)
databarcodes <- merge(databarcodes, metadata[,c("ID", "cluster")],
               by.x = "sample2", by.y = "ID", all.x = TRUE,
               suffixes = c("1", "2"))

databarcodes$sameCluster <- databarcodes$cluster1 == databarcodes$cluster2
databarcodes$commonCluster <- NA
databarcodes$commonCluster[databarcodes$sameCluster] <- databarcodes$cluster1[databarcodes$sameCluster]

combos <- strsplit(paste(databarcodes$cluster1[!databarcodes$sameCluster], databarcodes$cluster2[!databarcodes$sameCluster], sep =";"), split=";")
combos <- sapply(combos, sort)
if (length(combos) > 0){
  combos <- paste0(combos[1,], combos[2,])
  databarcodes$commonCluster[!databarcodes$sameCluster] <- combos
}

##################################

###########################SAVE DATA###########################
write.table(metadata,
            file = paste0(outDir, "/nodes.tsv"),
            sep ="\t",
            quote = F, row.names = F)
write.table(databarcodes[databarcodes$fract_sites_IBD>-1,],
            file = paste0(outDir, "/edges.tsv"),
            sep ="\t",
            quote = F, row.names = F)
##################################