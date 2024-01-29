#startPosindex and endPosIndex delimit the locations of the columns storing
#the bases called by SNP genotyping or WGS
startPosindex <- 2
endPosIndex <- ncol(barcodes)-1
######################CONVERT FORMAT FOR HMMIBD###############
#All positions that are unknown or mixed are ignored
cv4ibd <- function(sample, ref, alt, alt2, alt3){
  sample <- as.vector(sample)
  sample[sample=="N"] <- -1
  sample[sample=="-"] <- -1
  sample[sample=="X"] <- -1
  sample[sample==ref] <- 0
  sample[sample==alt] <- 1
  sample[sample==alt2] <- 2
  sample[sample==alt3] <- 3
  return(sample)
}

#FOR REF AND ALT
hmmibdrefalt <- t(barcodes[barcodes$ID%in%c("REF", "ALT", "ALT2", "ALT3"),])

colnames(hmmibdrefalt) <- hmmibdrefalt[1,]
hmmibdrefalt <- data.frame(hmmibdrefalt[startPosindex:endPosIndex,])
hmmibdrefalt$chrom <- unlist(lapply(rownames(hmmibdrefalt),
                                    FUN = function(x) { unlist(strsplit(x, split = "[.]"))[1] }))
hmmibdrefalt$pos <- unlist(lapply(rownames(hmmibdrefalt),
                                  FUN = function(x) { unlist(strsplit(x, split = "[.]"))[2] }))

if(any(is.na(chrOrder))){
  chrOrder <- sort(unique(hmmibdrefalt$chrom))
}

hmmibdrefalt$chrom <- as.numeric(factor(hmmibdrefalt$chrom, levels = chrOrder))
hmmibdrefalt <- hmmibdrefalt[!is.na(hmmibdrefalt$chrom),]
hmmibdrefalt$pos <- as.integer(hmmibdrefalt$pos)

#FOR BARCODES
hmmibdbarcodes <- t(barcodes)
colnames(hmmibdbarcodes) <- hmmibdbarcodes[1,]
hmmibdbarcodes <- data.frame(hmmibdbarcodes[startPosindex:endPosIndex,])
hmmibdbarcodes$chrom <- unlist(lapply(rownames(hmmibdbarcodes),
                                      FUN = function(x) { unlist(strsplit(x, split = "[.]"))[1] }))
hmmibdbarcodes$pos <- unlist(lapply(rownames(hmmibdbarcodes),
                                    FUN = function(x) { unlist(strsplit(x, split = "[.]"))[2] }))

hmmibdbarcodes$chrom <- as.numeric(factor(hmmibdbarcodes$chrom, levels = chrOrder))
hmmibdbarcodes <- hmmibdbarcodes[!is.na(hmmibdbarcodes$chrom),]
hmmibdbarcodes$pos <- as.integer(hmmibdbarcodes$pos)

hmmibdbarcodes <- cbind(hmmibdbarcodes[,c((ncol(hmmibdbarcodes)-1):ncol(hmmibdbarcodes))],
                        hmmibdbarcodes[,-c((ncol(hmmibdbarcodes)-1):ncol(hmmibdbarcodes))])

hmmibdbarcodes[,-c(1:2)] <- apply(X = hmmibdbarcodes[,-c(1:2)], MARGIN = 2,
                                  FUN = cv4ibd, ref=hmmibdrefalt$REF, alt=hmmibdrefalt$ALT, 
                                  alt2=hmmibdrefalt$ALT2, alt3=hmmibdrefalt$ALT3)

hmmibdbarcodes <- hmmibdbarcodes[with(hmmibdbarcodes, order(chrom, pos)), ]
hmmibdbarcodes.refalt <- hmmibdbarcodes
hmmibdbarcodes <- hmmibdbarcodes[,!colnames(hmmibdbarcodes)%in%c("REF", "ALT", "ALT2", "ALT3")]
########################

############WRITE DATA############
write.table(hmmibdbarcodes, file = paste0(outDir, "/hmmIBD-barcodes.txt"), sep="\t", quote = F, col.names = T, row.names = F)
write.table(hmmibdbarcodes.refalt, file = paste0(outDir, "/hmmIBD-barcodes-refalt.txt"), sep="\t", quote = F, col.names = T, row.names = F)
########################
