#########
# Indel #
#########



rm(list=ls())
setwd("~/Genome_SNP_INDEL")


##########


femaleIndels <- read.table("Female.indels.tsv", header = T, sep = "\t")
maleIndels <- read.table("Male.indels.tsv", header = T, sep = "\t")
save.image()



##########



nrow(femaleIndels[femaleIndels$FILTER == "PASS",]) # 1,151,997
nrow(maleIndels[maleIndels$FILTER == "PASS",]) # 1,177,918



##########



femaleIndelsChromossome <- data.frame(summary(femaleIndels$CHROMOSOME[femaleIndels$FILTER == "PASS"]))
femaleIndelsChromossome

maleIndelsChromossome <- data.frame(summary(maleIndels$CHROMOSOME[maleIndels$FILTER == "PASS"]))
maleIndelsChromossome

inDelsChromossome <- cbind(femaleIndelsChromossome, maleIndelsChromossome)
names(inDelsChromossome) <- c("Female", "Male")
inDelsChromossome
inDelsChromossome$chrLength <- c(23037639, 18140952, 789605, 19818926, 282498, 22702307,
                              1566225, 24396255, 3268264, 30274277, 20304914, 22053297,
                              740079, 17126926, 829735, 29360087, 5170003, 24021853,
                              568933, 18779844, 19341862, 1220746, 23867706, 76237,
                              25021643, 421237, 21508407, 21026613, 1447032, 22385789,
                              23006712, 487831, 43154196)
inDelsChromossome

inDelsChromossome$FemaleFreq <- inDelsChromossome$Female/inDelsChromossome$chrLength
inDelsChromossome$MaleFreq <- inDelsChromossome$Male/inDelsChromossome$chrLength
inDelsChromossome

rm(femaleIndelsChromossome, maleIndelsChromossome)

InDelsFreqQ <- quantile(c(inDelsChromossome$FemaleFreq,inDelsChromossome$MaleFreq))

IndelsFreqOutSup <- InDelsFreqQ["75%"] + 1.5 * (InDelsFreqQ["75%"] - InDelsFreqQ["25%"])
names(IndelsFreqOutSup) <- NULL
IndelsFreqOutSup # 0.003773585

IndelsFreqOutInf <- InDelsFreqQ["25%"] - 1.5 * (InDelsFreqQ["75%"] - InDelsFreqQ["25%"])
names(IndelsFreqOutInf) <- NULL
IndelsFreqOutInf # 0.00088069

order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr1_random", "chr3_random", "chr4_random", "chr5_random", "chr7_random", "chr9_random", "chr10_random", "chr11_random", "chr12_random", "chr13_random", "chr16_random", "chr17_random", "chr18_random", "chrUn")

inDelsChromossome <- inDelsChromossome[match(order, row.names(inDelsChromossome)),]
inDelsChromossome

jpeg(
  filename = "InDels_per_Chr.jpeg",
  width = 700,
  height = 700,
  quality = 100
)
par(mar = c(5, 9, 1, 1) + 0.1)
barplot(
  t(data.frame(inDelsChromossome$FemaleFreq, inDelsChromossome$MaleFreq, row.names = row.names(inDelsChromossome))),
  beside = T,
  horiz = T,
  ylab = "",
  xlab = "",
  xlim = c(0, 0.0040),
  col = c("#aea36c", "#984807"),
  las = 1,
  cex.names = 0.9,
  cex.axis = 0.9
)
title(ylab="Chromossome", line=7, cex.lab=0.9)
title(xlab="Frequency", line=2, cex.lab=0.9)
title(xlab="(number of InDels per chromossome bp)", line=3, cex.lab=0.8)
abline(v = IndelsFreqOutSup, col = "red")
abline(v = IndelsFreqOutInf, col = "black")
legend(
  x = "topright",
  pch = 15,
  c("Male", "Female"),
  col = c("#984807", "#aea36c"),
  cex = 0.9,
  bty = "n"
)
dev.off()
# Alerts on 16_random (too low) There are any too hight



##########



inDelsChromossome$diference <- inDelsChromossome$Female / inDelsChromossome$Male

inDelsChromossome$diferencePlot <- rep(0, 33)
for (i in 1:nrow(inDelsChromossome)) {
  if (inDelsChromossome$diference[i] < 1) {inDelsChromossome$diferencePlot[i] <- -1/inDelsChromossome$diference[i]}
  else {inDelsChromossome$diferencePlot[i] <- inDelsChromossome$diference[i]}
}

inDelsChromossome

inDelsDifQ <- quantile(inDelsChromossome$diference)

inDelDifOutSup <- inDelsDifQ["75%"] + 1.5 * (inDelsDifQ["75%"] - inDelsDifQ["25%"])
names(inDelDifOutSup) <- NULL
inDelDifOutSup

inDelsDifOutInf <- inDelsDifQ["25%"] - 1.5 * (inDelsDifQ["75%"] - inDelsDifQ["25%"])
names(inDelsDifOutInf) <- NULL
inDelsDifOutInf <- -1/inDelsDifOutInf

jpeg(
  filename = "indel_foldchange_per_Chr.jpeg",
  width = 700,
  height = 700,
  quality = 100
)
par(mar = c(8, 4, 1, 1) + 0.1)
barplot(
  t(data.frame(inDelsChromossome$diferencePlot, row.names = row.names(inDelsChromossome))),
  ylab = "",
  xlab = "",
  ylim = c(-2, 2),
  las = 2,
  cex.names = 0.9,
  cex.axis = 0.9
)
title(xlab="Chromossome", line=6, cex.lab=0.9)
title(ylab="Female/male fold change", line=2.5, cex.lab=0.9)
abline(h = inDelDifOutSup, col = "red")
abline(h = inDelsDifOutInf, col = "blue")
text(40, 1.75, "♀", cex = 2)
text(40, -1.75, "♂", cex = 2)
dev.off()
# Outliers much more on...
# ... Male: chr5_random
# ... Female: chr3_random, chr11_random, chr13_random, chr16_random, chr17_random

save.image()



##########



chromossomes <- row.names(inDelsChromossome)

join <- data.frame(a = 1:100)
for (i in 1:length(chromossomes)) {
  cat(paste("\n ...", chromossomes[i]))
  join <- cbind(join, chrHeatmap(chromossomes[i], inDelsChromossome$chrLength[i], femaleIndels, maleIndels, 100), Na = rep (0, 100))
}

join$a <- NULL

jpeg(
  filename = "indel_chr_heatmaps_scale_none.jpeg",
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
  filename = "indel_chr_heatmaps_scale.jpeg",
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

max(join)
View(join[, 22:ncol(join)])

##########

join

chr <- 8
start <- 18
end <- 20
paste(round(inDelsChromossome$chrLength[chr]/100 * start), round(inDelsChromossome$chrLength[chr]/100 * end), sep="..")



##########



inDelsImpact <- cbind(data.frame(summary(femaleIndels$IMPACT[femaleIndels$FILTER == "PASS"])), data.frame(summary(maleIndels$IMPACT[maleIndels$FILTER == "PASS"])))
names(inDelsImpact) <- c("Female", "Male")
inDelsImpact

barplot(
  t(inDelsImpact),
  beside = T,
  main = "Number os InDels according the impact",
  ylab = "Number of ocurrences",
  xlab = "Impact",
  col = c("Red", "Blue")
)
# It seems that there are no differences between male and female
# but it is super clear that MODIFIER is the most relevant class.



##########



inDelsEffect <- cbind(data.frame(summary(femaleIndels$EFFECT[femaleIndels$FILTER == "PASS"])), data.frame(summary(maleIndels$EFFECT[maleIndels$FILTER == "PASS"])))
names(inDelsEffect) <- c("Female", "Male")
row.names(inDelsEffect)
row.names(inDelsEffect) <- c("Codon change + deletion", "Codon change + insertion",
  "Codon deletion", "Codon insertion", "Downstream", "Exon deleted", "Frame shift",
  "Intergenic", "Intragenic", "Intron", "Splice site acceptor", "Splice site donor",
  "Start lost", "Stop gained", "Stop lost", "Upstream", "UTR (3')", "UTR (5')", "")
inDelsEffect[1:18,]
# No significative differences are found between sexes

jpeg(
  filename = "inDels_per_effect.jpeg",
  width = 700,
  height = 700,
  quality = 100
)
par(mar = c(4, 10, 1, 2) + 0.1)
barplot(
  t(inDelsEffect[1:18,]),
  beside = T,
  horiz = T,
  xlim = c(0,500000),
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



##########



genotype <- data.frame(Female = summary(femaleIndels$Female.GENOTYPE[femaleIndels$FILTER == "PASS"]), Male = summary(maleIndels$Male.GENOTYPE[maleIndels$FILTER == "PASS"]))
genotype

#    Female   Male
#0/1 701974 752511
#1/1 416953 389688
#1/2  33070  35719

jpeg(
  filename = "indels_homozigotic_levels.jpeg",
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



##########



library(VennDiagram)

femaleInDelsVenn <- paste(femaleIndels$CHROMOSOME[femaleIndels$FILTER == "PASS"], femaleIndels$POSITION[femaleIndels$FILTER == "PASS"], sep = "_")
maleInDelsVenn <- paste(maleIndels$CHROMOSOME[maleIndels$FILTER == "PASS"], maleIndels$POSITION[maleIndels$FILTER == "PASS"], sep = "_")

overlap <- calculate.overlap(
  x = list(femaleInDelsVenn, maleInDelsVenn)
)

length(overlap$a1)
length(femaleInDelsVenn)

length(overlap$a2)
length(maleInDelsVenn)

length(overlap$a3)

jpeg(
  filename = "inDels_venn_ocurrences.jpeg",
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