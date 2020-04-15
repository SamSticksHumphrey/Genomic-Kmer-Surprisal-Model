args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

library(data.table)
species <-  args[1]
seqType <- args[2]

getRedundancy <- function(k, species){
  
  print(paste0("K = ", k, " ---------------------"))
  
  data <- as.data.frame(fread(paste0("/data/rnabiology/shumphrey/Project/3-Distributions/", seqType, "/", seqType, "dist_", species, "_K", k, ".txt")))
  colnames(data) <- c("seq", "occurance")
  totalOcc <- sum(data$occurance)
  data$prob <- data$occurance/totalOcc
  data$surp <- -log2(data$prob)
  
  entropy <- sum(data$prob * data$surp)
  redundancy <- 1 - (entropy/log2(20^(k/3)))
  
  individualKmers <- length(data$seq)
  
  speciesOut <- c(species, k, individualKmers, totalOcc, entropy,  redundancy)
  
  return(speciesOut)
}

output <- lapply(seq(3,15, 3), getRedundancy, species = species) 


speciesInfo <- data.frame(matrix(unlist(output), nrow=length(output), byrow=T),stringsAsFactors=FALSE)
names(speciesInfo) <- c("species", "k", "individualKmers", "totalKmers", "entropy", "redundancy")


write.table(speciesInfo, paste0("/home/shumphrey/GenomicInformationAnalysis/Data/Redundancy/", species, "_", seqType, "_RedundancyInfoTable_AA.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


