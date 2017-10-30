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

snp <- read.delim(file="snps.txt", header = TRUE, sep = '\t')
snp$Chromosome[snp$Chromosome==23] <- "X"
snp <- mutate(snp, seqnames = paste("chr", Chromosome, sep = ''))
labels <- t(matrix(unlist(strsplit(as.character(snp$SNP), ":")), nrow = 4))
snp$Chromosome <- NULL
snps <- GRanges(seqnames = snp$seqnames,
                ranges = IRanges(start = snp$Position,
                                 end = snp$Position),
                name = labels[,1],
                baseline = snp$Baseline,
                effect = snp$Effect,
                genome = "hg19")

my.files <- list.files("biofeatures", full.names = TRUE)
biofeatures <- GetBioFeatures(files = my.files, genome = "hg19")

# AOC "asion ovarian cancer"
enrichAOC <- CalculateEnrichment(variants = snps,
                                 features = biofeatures,
                                 feature.type = 
                                 CI = 0.95,
                                 return.overlaps = TRUE)
