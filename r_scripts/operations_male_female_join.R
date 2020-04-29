###################################################
# Analyse SNP and INDEL data from genome sequence #
# Operations Set -- Male and Female join          #
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

# Keep only the pass SNPs
femaleSnpPass <- femaleSnp[femaleSnp$FILTER == "PASS",]
maleSnpPass <- maleSnp[maleSnp$FILTER == "PASS",]
save.image()

# Number of SNPs
nrow(femaleSnpPass) # 6,197,145
nrow(maleSnpPass) # 6,575,026

# See head
head(femaleSnpPass)

# Remove non required columns
femaleSnpPass['DBSNP.ID'] <- NULL
femaleSnpPass['QD'] <- NULL
femaleSnpPass['Female.GQ'] <- NULL

maleSnpPass['DBSNP.ID'] <- NULL
maleSnpPass['QD'] <- NULL
maleSnpPass['Male.GQ'] <- NULL


# Rename columns
names(femaleSnpPass) <- c('chromosome', 'position', 'reference', 'observed', 'score', 'filter', 'A.C.G.T', 'coverage', 'amino_acid', 'codon', 'effect', 'exon', 'functional_class', 'gene', 'impact', 'transcript', 'genotype', 'allele_depth', 'allele_balance', 'coverage')
head(femaleSnpPass)

names(maleSnpPass) <- c('chromosome', 'position', 'reference', 'observed', 'score', 'filter', 'A.C.G.T', 'coverage', 'amino_acid', 'codon', 'effect', 'exon', 'functional_class', 'gene', 'impact', 'transcript', 'genotype', 'allele_depth', 'allele_balance', 'coverage')
head(maleSnpPass)



# Add sex information
femaleSnpPass['sex'] <- "Female"
maleSnpPass['sex'] <- "Male"
head(femaleSnpPass)
head(maleSnpPass)
save.image()


# Join Female + Male
wildSnpPass <- rbind(femaleSnpPass, maleSnpPass)
head(wildSnpPass)
tail(wildSnpPass)
save.image()


# Code SNPs
wildSnpPass['snp_code'] <- paste(wildSnpPass$chromosome, wildSnpPass$position, sep="_")
head(wildSnpPass)
tail(wildSnpPass)
save.image()


# 






femaleSnpChromossome <- data.frame(summary(femaleSnpPass$chromosome))
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

