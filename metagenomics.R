getwd()
setwd("/home/brenda/eanbit_cohort4/programming/Residential_training/metagenomics/")
install.packages("vegan")
library(vegan)
install.packages("BiocManager")
library(phyloseq)
BiocManager::install("phyloseq")
library(phyloseq)
BiocManager::install("Biostrings")
library(Biostrings)
BiocManager::install("DECIPHER")
library(DECIPHER)
BiocManager::install("dada2")
library(dada2)
library(ggplot2)
getwd()
## confirm that you downloaded the sequences
path <- "/home/brenda/eanbit_cohort4/programming/Residential_training/metagenomics/miseqsopdata/MiSeq_SOP/"
list.files(path)

##Quality trimming and filtering
###Sort to get matched lists of forward and reverse files
####forward reads
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnFs

####reverse reads
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
fnRs

###Extract sample names from the first element for the file names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

plotQualityProfile(fnRs[7:10])

##creating subdirectories
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz")) 
filtFs

filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
filtRs

names(filtFs) <- sample.names
names(filtFs)
names(filtRs) <- sample.names
names(filtRs)

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE) # On Windows, multithread=FALSE head(out)
out

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

##Dereplicate to get unique sequences and save on computation time downstream
derep_forward <- derepFastq(filtFs, verbose=TRUE)
derep_reverse <- derepFastq(filtRs, verbose=TRUE)

names(derep_forward) <- sample.names
names(derep_reverse) <- sample.names

##Apply the core sequence-variant inference algorithm to the dereplicated sequences.
dadaFs <- dada(derep_forward, err=errF, multithread=TRUE)
dadaFs
dadaRs <- dada(derep_reverse, err=errR, multithread=TRUE)
dadaRs
##How many sequence variants were inferred
dadaFs[[1]]

##merging forward and reverse reads
merged_reads <- mergePairs(dadaFs, derep_forward, dadaRs, derep_reverse, verbose=TRUE)
head(merged_reads)

#construct sequence table
seqtab <- makeSequenceTable(merged_reads)
View(seqtab)
dim(seqtab)

table(nchar(getSequences(seqtab)))

##removing chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
View(seqtab.nochim)
dim(seqtab.nochim)

##What percentage of our reads did we keep after removing chimeras?
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merged_reads, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

rownames(track) <- sample.names

head(track)

##Assigning taxonomy
getwd()
taxa <- assignTaxonomy(seqtab.nochim, "~/eanbit_cohort4/programming/Residential_training/metagenomics/miseqsopdata/MiSeq_SOP/silva_nr_v138_train_set.fa", multithread=TRUE)
View(taxa)
taxa <- addSpecies(taxa, "./miseqsopdata/MiSeq_SOP/silva_species_assignment_v138.fa")
View(taxa)
taxa.print <- taxa # Removing sequence rownames for display only 
taxa.print
rownames(taxa.print) <- NULL
head(taxa.print)

##Analysis and Visualization
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, 'D'), `[`, 1)
gender <- substr(subject, 1, 1)
subject <- substr(subject, 2, 999)
day <- as.integer(sapply(strsplit(samples.out, 'D'), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- 'Late'
rownames(samdf) <- samples.out

##Construct a phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), sample_data(samdf), tax_table(taxa))

ps <- prune_samples(sample_names(ps) !='Mock', ps)
ps

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
names(dna) <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

##rarefaction curve
rarecurve(otu_table(ps), step = 100, ylab = "Observed ASVs", xlab = "Number of reads", main = "Rarefaction curves", label = FALSE)

##Alpha diversity
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")

##Beta diversity
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))

ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")

##Barplot
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")
