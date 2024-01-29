##############NUMBER OF COMPARIBLE SITES##############
barcodes.num <- gsub("X", "0", barcodes$Barcode)
barcodes.num <- gsub("N", "0", barcodes.num)
barcodes.num <- gsub("(A|T|G|C)", "1", barcodes.num)
names(barcodes.num) <- barcodes$ID
barcodes.num <- melt(barcodes.num)
barcodes.comp.num <- expand.grid(row.names(barcodes.num), row.names(barcodes.num))
barcodes.comp.num$ID1 <- pmin(as.character(barcodes.comp.num$Var1), as.character(barcodes.comp.num$Var2))
barcodes.comp.num$ID2 <- pmax(as.character(barcodes.comp.num$Var1), as.character(barcodes.comp.num$Var2))
barcodes.comp.num <- barcodes.comp.num[,c(3,4)]
barcodes.comp.num <- barcodes.comp.num[!duplicated(barcodes.comp.num),]
barcodes.comp.num <- barcodes.comp.num[barcodes.comp.num$ID1!=barcodes.comp.num$ID2,]
barcodes.comp.num$N_comp_sites <- apply(barcodes.comp.num,1,
                                        FUN = function(x){
                                          bc1 <- as.numeric(do.call("rbind", str_split(barcodes.num[row.names(barcodes.num)==x[1],1], "")))
                                          bc2 <- as.numeric(do.call("rbind", str_split(barcodes.num[row.names(barcodes.num)==x[2],1], "")))
                                          return(length(which(bc1+bc2==2)))
                                        })
length(which(barcodes.comp.num$N_comp_sites>40))/nrow(barcodes.comp.num)
############################

##############SAVE FILES##############
write.csv(barcodes.comp.num, paste0(outDir, "/barcodes-comp-sites.csv"), quote = F, row.names = F)
############################
