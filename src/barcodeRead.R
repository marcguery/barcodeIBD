###############BARCODES###############
#Loading raw data files
barcodes <- read.csv(barcodefile)
names(barcodes)[1] <- "ID"
##############################

##################ADD REF AND ALT##################
refbarcode <- apply(barcodes[,c(2:ncol(barcodes))],2,
                    FUN = function(x){
                      alleles <- x[x!="N" & x != "X"]
                      if(length(unique(alleles))>=1){
                        base <- sort(unique(alleles))[1]
                      }else{
                        base <- "-"
                      }
                      base
                    })
altbarcode <- apply(barcodes[,c(2:ncol(barcodes))],2,
                    FUN = function(x){
                      alleles <- x[x!="N" & x != "X"]
                      if(length(unique(alleles))>=2){
                        base <- sort(unique(alleles))[2]
                      }else{
                        base <- "-"
                      }
                      base
                    })
alt2barcode <- apply(barcodes[,c(2:ncol(barcodes))],2,
                    FUN = function(x){
                      alleles <- x[x!="N" & x != "X"]
                      if(length(unique(alleles))>=3){
                        base <- sort(unique(alleles))[3]
                      }else{
                        base <- "-"
                      }
                      base
                    })
alt3barcode <- apply(barcodes[,c(2:ncol(barcodes))],2,
                    FUN = function(x){
                      alleles <- x[x!="N" & x != "X"]
                      if(length(unique(alleles))>=4){
                        base <- sort(unique(alleles))[4]
                      }else{
                        base <- "-"
                      }
                      base
                    })

refalt <- data.frame(t(data.frame(refbarcode,
                                  altbarcode,
                                  alt2barcode,
                                  alt3barcode)))
refalt$ID <- c("REF", "ALT", "ALT2", "ALT3")
barcodes <- rbind(barcodes, refalt)
####################################

############BUILD BARCODE############
barcodes$Barcode <- apply(barcodes[,-1], 1, function(x){paste(x, collapse = "")})

#Checking if a character in a barcode is unexpected
unexpectedchar <- apply(X = barcodes[,-which(colnames(barcodes)%in%c("ID", "Barcode"))], 
                        FUN = function (x) { x!="-" & x!="X" & x!="N" & x!="A" & x!="T" & x!="G" & x!="C" }, 
                        MARGIN = 2)
boolunexpectedchar <- apply(X = unexpectedchar, 
                            FUN = function(x) {any(x)}, 
                            MARGIN = 1)
which(boolunexpectedchar) #it should be empty
##############################

###############ID MAPPING###############
#Loading raw data files
metadata <- read.csv(metadatafile)
names(metadata)[1:4] <- c("ID", "individual", "date", "location")
#Number of blood samples
nrow(metadata)
#Number of participants
length(unique(metadata$individual))

#Showing duplicated IDs
which(duplicated(metadata$ID, fromLast = F))
which(duplicated(metadata$Sample.Internal.ID, fromLast = F))

#Simplify
maxcharnumoflocations<- nchar(as.character(max(as.numeric(factor(metadata$location)))))
metadata$shortlocation <- sprintf(paste0("%0",maxcharnumoflocations,"d"), 
                                  as.numeric(factor(metadata$location)))

#Add number of Ns and Xs
n.num <- apply(X = barcodes[,-which(colnames(barcodes)%in%c("ID", "Barcode"))], 
                        FUN = function (x) { length(which(x=="N"))}, 
                        MARGIN = 1)
x.num <- apply(X = barcodes[,-which(colnames(barcodes)%in%c("ID", "Barcode"))], 
               FUN = function (x) { length(which(x=="X"))}, 
               MARGIN = 1)
barcodenum <- data.frame("ID" = barcodes$ID,
                         "Ns" = n.num, "Xs" = x.num)

metadata <- merge(metadata, barcodenum)

#IDs that are not barcodes (should be only REF, ALT, ALT2 and ALT3)
barcodes$ID[!barcodes$ID%in%metadata$ID]
###############################################

##################SAVING FILES##############
write.csv(barcodes, paste0(outDir, "/barcodes.csv"),
          quote = F, row.names = F) #Barcodes without any filter
write.csv(metadata, file = paste0(outDir, "/barcodes-metadata.csv"), quote = F, row.names = F)
###############################################
