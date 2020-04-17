library(biomaRt)
library(ggplot2)
library(GenomicRanges)
library(pheatmap)
library(tidyr)
library(dplyr)
library(rstan)
library(ggjoy)
library(ggthemes)
library(SummarizedExperiment)
library(funciVar)
setwd("/Users/hazelettd/Projects/GWAS/OVARIAN_CANCER/ASIAN_STUDY/")

fg <- read.delim(file="bfdp05.list", header = TRUE) # foreground list
fg <- mutate(fg, chrom_label=paste("chr", Chromosome, sep=""))
fg$chrom_label[fg$chrom_label=="chr23"] <- "chrX"
fg.gr <- GRanges(seqnames = fg$chrom_label,
                 ranges = IRanges(start = fg$Position,
                                  end = fg$Position),
                 strand = "*")
fg.gr$blockStarts <- start(fg.gr)
fg.gr$blockStarts <- start(fg.gr)
fg.gr$blockSizes <- end(fg.gr)-start(fg.gr)+1
## for bdfp < 0.5 SNPs "effect"=allele1, "baseline"=allele2
values(fg.gr) <- DataFrame(label=fg$SNP, linkage_group=fg$Chromosome, 
                           allele1=fg$Effect, 
                           allele2=fg$Baseline,
                           blockStarts=start(fg.gr),
                           blockSizes=end(fg.gr)-start(fg.gr))

bg <- read.delim(file = "Selected_SNP_regions.txt", header = TRUE) # background list
bg <- mutate(bg, chrom_label=paste("chr", Chromosome, sep = ""))
bg$chrom_label[bg$chrom_label=="chr23"] <- "chrX"
bg.gr <- GRanges(seqnames = bg$chrom_label,
                 ranges = IRanges(start = bg$Position,
                                  end = bg$Position),
                 strand = "*")
values(bg.gr) <- DataFrame(label=bg$SNP, linkage_group=bg$Chromosome, 
                           allele1 = bg$A1, 
                           allele2 = bg$A2, 
                           blockStarts=start(bg.gr),
                           blockSizes=end(bg.gr)-start(bg.gr))

background.bed <- list(seqnames = seqnames(bg.gr), starts = start(bg.gr), ends = end(bg.gr), names = mcols(bg.gr)$label)
cat("browser position chr1:157878955-158078955\ntrack name='ASIAN BKG SNPS' description='Background SNPs for Ovarian Cancer study in Asian population'\n", file = "asian_oc_background.hg19.bed")
write.table(background.bed, file = "asian_oc_background.hg19.bed", sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)

aoc.variants <- list(fg=fg.gr, bg=bg.gr)

annofiles <- list.files("biofeatures", full.names = TRUE)
annos <- data.frame(read.delim("Asian_samplenames.csv", sep = ',', row.names = 1))

# store the indices of each feature type for subsets analysis
h3k4me1  <- c(3, 11, 13, 29, 32, 35, 38, 41, 75)
h3k27ac  <- c(2, 4:8, 10, 15, 17, 18, 20, 22:25, 26, 28, 31, 34, 37, 40, 49, 51:56, 60:65, 68:72, 74)
h3k4me3  <- c(12, 16, 21)
dhsfaire <- c(1, 9, 14, 27, 30, 33, 36, 39, 73)
ctcf     <- c(19, 48, 58, 59, 66, 67)
pax8     <- c(42:47, 50, 57)

biofeatures <- GetBioFeatures(files = annofiles[h3k27ac], genome = "hg19") # hang on! this will take a while

enrichment <- CalculateEnrichment(aoc.variants, biofeatures, feature.type="biofeatures", CI = 0.95, return.overlaps = TRUE)

global.plot <- PlotEnrichment(enrichment$enrichment, value = "difference", color.by = NULL) +
  ggtitle("Global Enrichment of H3K27ac")
ggsave(filename = "plots/h3k27ac_global_enrichment_plot.pdf", plot = global.plot, device = "pdf", width = 18, height = 12, units = "in")

overlap.mat <- mcols(enrichment$overlaps$foreground.overlaps)
write.csv(overlap.mat, file="foreground_overlaps.csv", quote = FALSE, row.names = FALSE)
fg.track <- list(chrom=seqnames(fg.gr), start = start(fg.gr), end = end(fg.gr), label = mcols(fg.gr)$label)
cat("browser position chr6:3571487-3588186\ntrack name='OVC BDFP 0.5' description='Ovarian cancer genetic associations (Bayes False Discovery Prob < 0.5)' visibility=pack\n", file = "foreground_snps.bed")
write.table(fg.track, file = "foreground_snps.bed", sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE, col.names = FALSE)

locus <- with(read.delim(file="regions.txt", sep="\t"), GRanges(seqnames = seqnames, ranges = IRanges(start = start, end = end), region = region))
enrichment$enrichment$sample


#define a function to create an enrichment analysis for all SNPs that overlap an interval in the locus table
length(h3k27ac)
enrichment.locus <- list(probability = data.frame(matrix(NA, nrow = 73, ncol = 30)),
                         difference  = data.frame(matrix(NA, nrow = 73, ncol = 30)))
colnames(enrichment.locus$probability) <- c("global", as.character(locus$region))
colnames(enrichment.locus$difference)  <- c("global", as.character(locus$region))
rownames(enrichment.locus$probability) <- rownames(enrichment$enrichment)
rownames(enrichment.locus$difference) <- rownames(enrichment$enrichment)
enrichment.locus$probability$global <- enrichment$enrichment$probability
enrichment.locus$difference$global  <- enrichment$enrichment$difference
enrichment.locus



for (a in 1:length(locus)) {
  print(paste(as.character(locus[a]$region)[1], "Enrichment", sep=" "))
  locus.variants <- list(fg=fg.gr[subjectHits(findOverlaps(locus[a], fg.gr))],
                         bg=bg.gr[subjectHits(findOverlaps(locus[a], bg.gr))])
  a.locus <- CalculateEnrichment(locus.variants, biofeatures, feature.type="biofeatures", CI = 0.95, return.overlaps = FALSE)
  a.locus$probability
  enrichment.locus$probability[rownames(a.locus),a+1] <- a.locus$probability
  enrichment.locus$difference[rownames(a.locus),a+1] <- a.locus$difference
  #a.plot <- PlotEnrichment(a.locus, value = "difference", color.by = NULL) + ggtitle(paste(locus[a]$region[1], "enrichment", sep=" ")) 
  #ggsave(filename = paste("plots/", locus[a]$region[1], ".pdf", sep=""), plot = a.plot, device = "pdf", width = 18, height = 12, units = "in")
  #ggsave(filename = paste("plots/", locus[a]$region[1], ".png", sep=""), plot = a.plot, device = "png", dpi = 92, width = 18, height = 12, units = "in")
}

load("asiandata.Rda")
save(enrichment.locus, file="asiandata.Rda")
dim(enrichment.locus$probability)
dim(annos)
png("plots/locus_heat_difference.png", res = 300, width = 12, height = 13, units = "in")
hm1 <- pheatmap(enrichment.locus$difference, annotation_row = annos[,c("Classification", "Feature", "Tissue.CellLine")])
pheatmap(enrichment.locus$difference, annotation_row = annos[,c("Classification", "Feature", "Tissue.CellLine")], labels_row = as.character(annos$SampleName)[hm1$tree_row$order])
dev.off()
png("plots/locus_heat_probability.png", res = 300, width = 12, height = 13, units = "in")
svg("plots/locus_heat_probability.svg", height = 13, width = 12)
hm2 <- pheatmap(enrichment.locus$probability, annotation_row = annos[,c("Classification", "Feature", "Tissue.CellLine")])
pheatmap(enrichment.locus$probability, annotation_row = annos[,c("Classification", "Feature", "Tissue.CellLine")], labels_row = as.character(annos$SampleName)[hm2$tree_row$order])
dev.off()

 
 
 
 
 rs6902488snp <- snps.from.rsid("rs6902488", dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh37, search.genome = BSgenome.Hsapiens.UCSC.hg19)
 
 rs6902488 <- motifbreakR(snpList = rs6902488snp, filterp = T,
                          pwmList = encodemotif,
                          threshold = 5e-5,
                          method = "ic",
                          bkg = c(A=0.28, C=0.22, G=0.22, T=0.28),
                          BPPARAM = BiocParallel::bpparam())
 
 rs6902488 <- calculatePvalue(rs6902488, background = c(A=0.28, C=0.22, G=0.22, T=0.28))
 ?pdf
 pdf(file = "plots/rs6902488_MBplot.pdf")
 plotMB(rs6902488, rsid = "rs6902488")
 dev.off()
 