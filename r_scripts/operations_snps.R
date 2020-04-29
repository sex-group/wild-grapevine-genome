###################################################
# Analyse SNP and INDEL data from genome sequence #
# Operations Set                                  #
###################################################

# Clear all environment
rm(list=ls())

# Change working directory
setwd("~/Genome_SNP_INDEL")



########
# SNPs #
########

# Load SNPs file
femaleSnp <- read.table("Female.snps.tsv", header = T, sep = "\t")
maleSnp <- read.table("Male.snps.tsv", header = T, sep = "\t")
save.image()



# Total number of SNPs
nrow(femaleSnp[femaleSnp$FILTER == "PASS",])
# [1] 6,197,145
nrow(maleSnp[maleSnp$FILTER == "PASS",])
# [1] 6,575,026



# Ocurrences per chromossome
femaleSnpChromossome <- data.frame(summary(femaleSnp$CHROMOSOME[femaleSnp$FILTER == "PASS"]))
femaleSnpChromossome

maleSnpChromossome <- data.frame(summary(maleSnp$CHROMOSOME[maleSnp$FILTER == "PASS"]))
maleSnpChromossome

snpChromossome <- cbind(femaleSnpChromossome, maleSnpChromossome)
names(snpChromossome) <- c("Female", "Male")
snpChromossome

row.names(snpChromossome)
snpChromossome$chrLength <- c(23037639, 18140952, 789605, 19818926, 282498, 22702307,
                              1566225, 24396255, 3268264, 30274277, 20304914, 22053297,
                              740079, 17126926, 829735, 29360087, 5170003, 24021853,
                              568933, 18779844, 19341862, 1220746, 23867706, 76237,
                              25021643, 421237, 21508407, 21026613, 1447032, 22385789,
                              23006712, 487831, 43154196)
snpChromossome

snpChromossome$FemaleFreq <- snpChromossome$Female/snpChromossome$chrLength
snpChromossome$MaleFreq <- snpChromossome$Male/snpChromossome$chrLength
snpChromossome

rm(femaleSnpChromossome, maleSnpChromossome)

# Total number
barplot(
  t(data.frame(snpChromossome$Female,snpChromossome$Male, row.names = row.names(snpChromossome))),
  beside = T,
  main = "Number os SPNs per Chromossome",
  ylab = "Ocurrences",
  xlab = "Chromossome",
  col = c("Red", "Blue"),
  las = 3
)
# Major differences are observed on chromossomes 1 and 10 (small differences)
# 13, 18, 5 and 6 (major differences)
# 19 (female has much more than male)

# Frequencies
barplot(
  t(data.frame(snpChromossome$FemaleFreq,snpChromossome$MaleFreq, row.names = row.names(snpChromossome))),
  beside = T,
  main = "Number os SPNs per Chromossome",
  ylab = "Ocurrences",
  xlab = "Chromossome",
  col = c("Red", "Blue"),
  las = 3
)
# Main differences in 11_random, 12_random, 13, 13_random, 16_random, 17_random, 3_random, 5_random, 6 and 9_random

femaleFreqQ <- quantile(snpChromossome$FemaleFreq)

femaleFreqOutSup <- femaleFreqQ["75%"] + 1.5 * (femaleFreqQ["75%"] - femaleFreqQ["25%"])
names(femaleFreqOutSup) <- NULL
femaleFreqOutSup

femaleFreqOutInf <- femaleFreqQ["25%"] - 1.5 * (femaleFreqQ["75%"] - femaleFreqQ["25%"])
names(femaleFreqOutInf) <- NULL
femaleFreqOutInf

barplot(
  t(data.frame(snpChromossome$FemaleFreq, row.names = row.names(snpChromossome))),
  beside = T,
  main = "Number os SPNs per Chromossome (Female Only)",
  ylab = "Ocurrences",
  xlab = "Chromossome",
  las = 3
)
abline(h=femaleFreqOutSup, col="red")
abline(h=femaleFreqOutInf, col="blue")
# Alerts on 16_random (too low) and 9_random (too hight)

maleFreqQ <- quantile(snpChromossome$MaleFreq)

maleFreqOutSup <- maleFreqQ["75%"] + 1.5 * (maleFreqQ["75%"] - maleFreqQ["25%"])
names(maleFreqOutSup) <- NULL
maleFreqOutSup

maleFreqOutInf <- maleFreqQ["25%"] - 1.5 * (maleFreqQ["75%"] - maleFreqQ["25%"])
names(maleFreqOutInf) <- NULL
maleFreqOutInf

barplot(
  t(data.frame(snpChromossome$MaleFreq, row.names = row.names(snpChromossome))),
  beside = T,
  main = "Number os SPNs per Chromossome (Male Only)",
  ylab = "Ocurrences",
  xlab = "Chromossome",
  las = 3
)
abline(h=maleFreqOutSup, col="red")
abline(h=maleFreqOutInf, col="blue")
# Alerts on 16_random (too low) and 9_random (too hight)


snpFreqQ <- quantile(c(snpChromossome$FemaleFreq,snpChromossome$MaleFreq))

snpFreqOutSup <- snpFreqQ["75%"] + 1.5 * (snpFreqQ["75%"] - snpFreqQ["25%"])
names(snpFreqOutSup) <- NULL
snpFreqOutSup

snpFreqOutInf <- snpFreqQ["25%"] - 1.5 * (snpFreqQ["75%"] - snpFreqQ["25%"])
names(snpFreqOutInf) <- NULL
snpFreqOutInf

order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr1_random", "chr3_random", "chr4_random", "chr5_random", "chr7_random", "chr9_random", "chr10_random", "chr11_random", "chr12_random", "chr13_random", "chr16_random", "chr17_random", "chr18_random", "chrUn")

snpChromossome <- snpChromossome[match(order, row.names(snpChromossome)),]

jpeg(
  filename = "SNPs_per_Chr.jpeg",
  width = 700,
  height = 700,
  quality = 100
  )
par(mar = c(5, 9, 1, 1) + 0.1)
barplot(
  t(data.frame(snpChromossome$FemaleFreq, snpChromossome$MaleFreq, row.names = row.names(snpChromossome))),
  beside = T,
  horiz = T,
  ylab = "",
  xlab = "",
  xlim = c(0, 0.025),
  col = c("#aea36c", "#984807"),
  las = 1,
  cex.names = 0.9,
  cex.axis = 0.9
)
title(ylab="Chromossome", line=7, cex.lab=0.9)
title(xlab="Frequency", line=2, cex.lab=0.9)
title(xlab="(number of SNPs per chromossome bp)", line=3, cex.lab=0.8)
abline(v = maleFreqOutSup, col = "red")
abline(v = maleFreqOutInf, col = "black")
legend(
  x = "topright",
  pch = 15,
  c("Male", "Female"),
  col = c("#984807", "#aea36c"),
  cex = 0.9,
  bty = "n"
)
dev.off()

# Alerts on 16_random (too low) and 9_random (too hight)


##############################
# Male vs Female fold change #
##############################
snpChromossome$diference <- snpChromossome$Female / snpChromossome$Male

snpChromossome$diferencePlot <- rep(0, 33)
for (i in 1:nrow(snpChromossome)) {
  if (snpChromossome$diference[i] < 1) {snpChromossome$diferencePlot[i] <- -1/snpChromossome$diference[i]}
  else {snpChromossome$diferencePlot[i] <- snpChromossome$diference[i]}
}

snpChromossome

snpDifQ <- quantile(snpChromossome$diference)

snpDifOutSup <- snpDifQ["75%"] + 1.5 * (snpDifQ["75%"] - snpDifQ["25%"])
names(snpDifOutSup) <- NULL
snpDifOutSup

snpDifOutInf <- snpDifQ["25%"] - 1.5 * (snpDifQ["75%"] - snpDifQ["25%"])
names(snpDifOutInf) <- NULL
snpDifOutInf <- -1/snpDifOutInf



jpeg(
  filename = "foldchange_per_Chr.jpeg",
  width = 700,
  height = 700,
  quality = 100
)
par(mar = c(8, 4, 1, 1) + 0.1)
barplot(
  t(data.frame(snpChromossome$diferencePlot, row.names = row.names(snpChromossome))),
  ylab = "",
  xlab = "",
  ylim = c(-2, 2),
  las = 2,
  cex.names = 0.9,
  cex.axis = 0.9
)
title(xlab="Chromossome", line=6, cex.lab=0.9)
title(ylab="Female/male fold change", line=2.5, cex.lab=0.9)
abline(h = snpDifOutSup, col = "red")
abline(h = snpDifOutInf, col = "blue")
text(40, 1.75, "♀", cex = 2)
text(40, -1.75, "♂", cex = 2)
dev.off()

save.image()



# Chromossome heatmap to all chromossomes, on both sexes
snpDensity(chrList = row.names(snpChromossome), snpListFemale = femaleSnp, snpListMale = maleSnp)
save.image()
# Plots are saved under its proper directory.
# There are some chromossomes that have some differences between male/female and may be further investigated
# The best way to draw those plots is using excel!


# SNP density by heatmap!
chromossomes <- row.names(snpChromossome)

join <- data.frame(a = 1:100)
for (i in 1:length(chromossomes)) {
  cat(paste("\n ...", chromossomes[i]))
  join <- cbind(join, chrHeatmap(chromossomes[i], snpChromossome$chrLength[i], femaleSnp, maleSnp, 100), Na = rep (0, 100))
}

join$a <- NULL

jpeg(
  filename = "chr_heatmaps_scale_none.jpeg",
  width = 700,
  height = 700,
  quality = 100
)
heatmap(
  x = as.matrix(t(join)),
  scale = "none",
  Colv = NA,
  Rowv = NA,
  col = rev(heat.colors(256)),
  labRow = NA,
  labCol = NA
)
dev.off()

jpeg(
  filename = "chr_heatmaps_scale_row.jpeg",
  width = 700,
  height = 700,
  quality = 100
)
heatmap(
  x = as.matrix(t(join)),
  scale = "row",
  Colv = NA,
  Rowv = NA,
  col = rev(heat.colors(256)),
  labRow = NA,
  labCol = NA
)
dev.off()

jpeg(
  filename = "chr_heatmaps_scale.jpeg",
  width = 700,
  height = 700,
  quality = 100
)
plot(
  x = seq(
    from = 1,
    to = max(join),
    by = max(join)/256
  ),
  y = rep(1, 256),
  col = rev(heat.colors(256)),
  pch = 15,
  cex = 3
)
dev.off()

join


chr <- 17
start <- 79.5
end <- 80.5
paste(round(snpChromossome$chrLength[chr]/100 * start), round(snpChromossome$chrLength[chr]/100 * end), sep="..")



# Type of SNPs
modificationsTable <- data.frame(
  female = rep(0, 12),
  male = rep(0, 12)
)

row.names(modificationsTable) <- c("AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG")

modificationsTable <- snpModificationTypes(modificationsTable, snpListFemale = femaleSnp, snpListMale = maleSnp)
row.names(modificationsTable) <- c("A»C", "A»G", "A»T", "C»A", "C»G", "C»T", "G»A", "G»C", "G»T", "T»A", "T»C", "T»G")
t(modificationsTable)
save.image()

jpeg(
  filename = "SNPs_per_modification_type.jpeg",
  width = 700,
  height = 700,
  quality = 100
)
par(mar = c(5, 5, 1, 1) + 0.1)
barplot(
  t(modificationsTable),
  beside = T,
  ylab = "Number of modifications",
  xlab = "Modification",
  col = c("#aea36c", "#984807"),
  ylim = c(0, 1600000),
  las = 0,
  cex.names = 0.9,
  cex.axis = 0.9
)
legend(
  x = "topright",
  pch = 15,
  c("Female", "Male"),
  col = c("#aea36c", "#984807"),
  cex = 0.9,
  bty = "n"
)
dev.off()

# No differences are observed between male/female but male has allways more
# As usual, transitions (A <-> G + C <-> T) are more frequent than transversions



# According to impact
snpImpact <- cbind(data.frame(summary(femaleSnp$IMPACT[femaleSnp$FILTER == "PASS"])), data.frame(summary(maleSnp$IMPACT[maleSnp$FILTER == "PASS"])))
names(snpImpact) <- c("Female", "Male")
snpImpact

barplot(
  t(snpImpact),
  beside = T,
  main = "Number os SPNs according the impact",
  ylab = "Number of ocurrences",
  xlab = "Impact",
  col = c("Red", "Blue")
)

# It seems that there are no differences between male and female
# but it is super clear that MODIFIER is the most relevant class.



# According Effect
snpEffect <- cbind(data.frame(summary(femaleSnp$EFFECT[femaleSnp$FILTER == "PASS"])), data.frame(summary(maleSnp$EFFECT[maleSnp$FILTER == "PASS"])))
names(snpEffect) <- c("Female", "Male")
row.names(snpEffect) <- c("Downstream", "Intergenic", "Intron", "Non synonumous", "Splice site acceptor", "Splice site donor", "Start gained", "Start lost", "Stop gained", "Stop lost", "Synonymous (coding)", "Synonymous (stop)", "Upstream", "UTR (3')", "UTR (5')", "")
snpEffect[1:15,]
# No significative differences are found between sexes


####
#### ATENÇÃO AOS NAS QUE FORAM REMOVIDOS
###

jpeg(
  filename = "SNPs_per_effect.jpeg",
  width = 700,
  height = 700,
  quality = 100
)
par(mar = c(4, 10, 1, 2) + 0.1)
barplot(
  t(snpEffect[1:15,]),
  beside = T,
  horiz = T,
  xlim = c(0,2500000),
  col = c("#aea36c", "#984807"),
  las = 1,
  cex.names = 0.9,
  cex.axis = 0.9
)
title(xlab="Number of ocurrences", line=2.5, cex.lab=0.9)
title(ylab="Effect", line=8.5, cex.lab=0.9)
legend(
  x = "topright",
  pch = 15,
  c("Male", "Female"),
  col = c("#984807", "#aea36c"),
  cex = 0.9,
  bty = "n"
)
dev.off()






# These results may be showed as a staced table or graphic, as effects belongs to a specific impact category



# Check the modification table for Amino acids
# Only to NON_SYNONYMOUS_CODING, STOP_GAINED or STOP_LOST
modificationsTableAA <- data.frame(
  Female = integer(),
  Male = integer()
)

modificationsTableAA <- aaModificationTypes(modificationsTableAA, snpListFemale = femaleSnp, snpListMale = maleSnp)
save.image()

modificationsTableAA

jpeg(
  filename = "SNPs_per_amino_acid.jpeg",
  width = 700,
  height = 700,
  quality = 100
)
par(mar = c(4, 4, 1, 2) + 0.1)
barplot(
  t(modificationsTableAA),
  beside = T,
  horiz = T,
  ylab = "Number of modifications",
  xlab = "Modification",
  col = c("#aea36c", "#984807"),
  las = 1,
  cex.names = 0.9,
  cex.axis = 0.9
)
legend(
  x = "topright",
  pch = 15,
  c("Male", "Female"),
  col = c("#984807", "#aea36c"),
  cex = 0.9,
  bty = "n"
)
dev.off()

####
#### Representar isto vai ser muito difícil. Talvez possa entrar como tabela em addional
####



# Check content of FUNCTIONAL.CLASS column
summary(femaleSnp$FUNCTIONAL.CLASS[femaleSnp$FILTER == "PASS"])

# Check content of GENE.NAME column
summary(femaleSnp$GENE.NAME[femaleSnp$FILTER == "PASS"])



# Check genotype proportions
genotype <- data.frame(Female = summary(femaleSnp$Female.GENOTYPE[femaleSnp$FILTER == "PASS"]), Male = summary(maleSnp$Male.GENOTYPE[maleSnp$FILTER == "PASS"]))
genotype

#     Female    Male
#0/1 3897545 4253285
#1/1 2271185 2291547
#1/2   28415   30194

jpeg(
  filename = "homozigotic_levels.jpeg",
  width = 700,
  height = 700,
  quality = 100
)
par(mar = c(4, 6, 1, 2) + 0.1)
barplot(
  t(genotype),
  beside = T,
  col = c("#aea36c", "#984807"),
  las = 1,
  cex.names = 0.9,
  cex.axis = 0.9
)
title(ylab="Number of loci", line=4.5, cex.lab=0.9)
title(xlab="Genotype", line=2.5, cex.lab=0.9)
legend(
  x = "topright",
  pch = 15,
  c("Female", "Male"),
  col = c("#aea36c", "#984807"),
  cex = 0.9,
  bty = "n"
)
dev.off()

save.image()




#####
##### VENN
#####

library(VennDiagram)
?draw.pairwise.venn

femaleSNP <- paste(femaleSnp$CHROMOSOME[femaleSnp$FILTER == "PASS"], femaleSnp$POSITION[femaleSnp$FILTER == "PASS"], sep = "_")
maleSNP <- paste(maleSnp$CHROMOSOME[maleSnp$FILTER == "PASS"], maleSnp$POSITION[maleSnp$FILTER == "PASS"], sep = "_")

overlap <- calculate.overlap(
  x = list(femaleSNP, maleSNP)
)

length(overlap$a1)
length(femaleSNP)

length(overlap$a2)
length(maleSNP)

length(overlap$a3)

jpeg(
  filename = "snps_venn_ocurrences.jpeg",
  width = 500,
  height = 400,
  quality = 100
)
draw.pairwise.venn(
  area1 = length(overlap$a1),
  area2 = length(overlap$a2),
  cross.area = length(overlap$a3),
  scaled = FALSE,
  category = c("Female", "Male"),
  cat.pos = c(0, 0),
  fill = c("#aea36c", "#984807")
)
dev.off()


################
## 2019.01.28 ##
################

# Frequency
names(femaleSnp)
chromosomes <- unique(femaleSnp$CHROMOSOME)
chromosomes <- as.character(chromosomes)

calculate_frequency <- function(positions) {
  distances <- c()
  
  for (chr in 1:length(positions)) {
    positions_chr <- positions[[chr]]
    
    # Loop to calculate distances
    for (snp in 1:(length(positions_chr)-1)) {
      next_snp <- snp + 1
      distance <- positions_chr[next_snp] - positions_chr[snp]
      distances <- c(distances, distance)
    }
  }
  
  # Compute
  avg_distance <- mean(distances, na.rm = TRUE)
  max_distance <- max(distances, na.rm = TRUE)
  min_distance <- min(distances, na.rm = TRUE)
  
  # Return elements
  output <- list("avg" = avg_distance,
                 "max" = max_distance,
                 "min" = min_distance,
                 "dist" = distances
                 )
  
  return(output)
}


statistics_snp_frequency <- function(sample_data, chromosomes) {
  positions <- list()
  for (chr in 1:length(chromosomes)) {
    # Compute
    positions[[chr]] <- sort(sample_data$POSITION[sample_data$CHROMOSOME == chromosomes[chr]])
  }
  
  # Calculate frequency
  statistics <- calculate_frequency(positions)
  
  # Output results
  return(statistics)
}



female_frequency <- statistics_snp_frequency(femaleSnp[femaleSnp$FILTER == "PASS",], chromosomes)
print(paste("For females, the (average) frequency is", female_frequency$avg, sep = " "))

male_frequency <- statistics_snp_frequency(maleSnp[maleSnp$FILTER == "PASS",], chromosomes)
print(paste("For males, the (average) frequency is", male_frequency$avg, sep = " "))



