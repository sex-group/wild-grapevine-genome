###################################################
# Analyse SNP and INDEL data from genome sequence #
# Functions Set                                   #
###################################################

# Perform a density plot analasis for all chromossomes of both sexes
snpDensity <- function(chrList, snpListFemale, snpListMale) {
  
  options(scipen=7)
  
  for (chr in chrList) {
    cat(paste("\n ...", chr))
    
    snpDensityFemale <- density(snpListFemale$POSITION[snpListFemale$FILTER == "PASS" & snpListFemale$CHROMOSOME == chr])
    snpDensityMale <- density(snpListMale$POSITION[snpListMale$FILTER == "PASS" & snpListMale$CHROMOSOME == chr])
    
    png(filename = paste("density/snp/", chr, ".png", sep = ""))
    plot(snpDensityFemale, xlab = "Position", ylim = range(snpDensityFemale$y, snpDensityMale$y), col="red", main = paste("SNP density plot to", chr))
    lines(snpDensityMale, col = "blue")
    legend("topright", legend = c("Female", "Male"), fill = c("red", "blue"))
    dev.off()
  }
}

chrHeatmap <- function (chr, chrSize, femaleSnp, maleSnp, steps) {
  femaleSNPchr <- femaleSnp[femaleSnp$CHROMOSOME == chr & femaleSnp$FILTER == "PASS" ,]
  maleSNPchr <- maleSnp[maleSnp$CHROMOSOME == chr & maleSnp$FILTER == "PASS" ,]
  
  ocurrencesFemale <- c()
  min <- 0
  for (i in 1:steps) {
    max <- chrSize / steps * i
    ocurrencesFemale <- c(ocurrencesFemale, nrow(femaleSNPchr[femaleSNPchr$POSITION >= min & femaleSNPchr$POSITION < max,]))
    min <- max
  }
  
  ocurrencesMale <- c()
  min <- 0
  for (i in 1:steps) {
    max <- chrSize / steps * i
    ocurrencesMale <- c(ocurrencesMale, nrow(maleSNPchr[maleSNPchr$POSITION >= min & maleSNPchr$POSITION < max,]))
    min <- max
  }
  
  join <- data.frame(Female = ocurrencesFemale, Male = ocurrencesMale)
  return(join)
}


# Perform a modification type analysis (number of ocurrences from A to C, etc), for both sexes
# Modification table should be empty
snpModificationTypes <- function(modificationsTable, snpListFemale, snpListMale) {
  
  for (modification in row.names(modificationsTable)) {
    cat(paste("\n ...", modification))
    
    modificationType <- strsplit(modification, "")
    
    modificationsTable[modification,] <- c(nrow(snpListFemale[snpListFemale$FILTER == "PASS" & snpListFemale$REFERENCE.BASE == modificationType[[1]][1] & snpListFemale$OBSERVED.BASE == modificationType[[1]][2], ]), nrow(snpListMale$FILTER == "PASS" & snpListMale[snpListMale$REFERENCE.BASE == modificationType[[1]][1] & snpListMale$OBSERVED.BASE == modificationType[[1]][2], ]))
  }
  
  return(modificationsTable)
}



# Perform a modification type analysis in termos of aminoacid change, for both sexes
# It will use only events whose effect are NON_SYNONYMOUS_CODING, STOP_GAINED or STOP_LOST
aaModificationTypes <- function(modificationsTableAA, snpListFemale, snpListMale) {
  
  snpSubListFemale <- snpListFemale[snpListFemale$FILTER == "PASS" & (snpListFemale$EFFECT == "NON_SYNONYMOUS_CODING" | snpListFemale$EFFECT == "STOP_GAINED" | snpListFemale$EFFECT == "STOP_LOST"),]
  snpSubListMale <- snpListMale[snpListMale$FILTER == "PASS" & (snpListMale$EFFECT == "NON_SYNONYMOUS_CODING" | snpListMale$EFFECT == "STOP_GAINED" | snpListMale$EFFECT == "STOP_LOST"),]
  
  cat("\n ... Executing for Females ...")
  modificationsTableAA <- aaModificationTypesForOneSex(modificationsTableAA, sexName = "Female", snpSubList = snpSubListFemale)

  cat("\n ... Executing for Males ...")
  modificationsTableAA <- aaModificationTypesForOneSex(modificationsTableAA, sexName = "Male", snpSubList = snpSubListMale)
  
  return(modificationsTableAA)
}

aaModificationTypesForOneSex <- function(modificationsTableAA, sexName, snpSubList) {
  
  for (row in 1:nrow(snpSubList)) {
    origin <- substr(snpSubList$AMINO.ACID.CHANGE[row],1,1)
    destination <- substr(snpSubList$AMINO.ACID.CHANGE[row],nchar(as.character(snpSubList$AMINO.ACID.CHANGE[row])),nchar(as.character(snpSubList$AMINO.ACID.CHANGE[row])))
    
    change <- paste(origin,destination, sep="")
    
    if (change != "NANA") {
      if (is.na(modificationsTableAA[change, sexName])) {modificationsTableAA[change, sexName] <- 1}
      else {modificationsTableAA[change, sexName] <- modificationsTableAA[change, sexName] + 1}
    }
  }
  
  return(modificationsTableAA)
}