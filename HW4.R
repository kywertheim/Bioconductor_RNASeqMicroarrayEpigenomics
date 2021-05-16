#The yeastRNASeq experiment data package contains FASTQ files from an RNA seq experiment in yeast.
library(yeastRNASeq)

#Q1. What fraction of reads in the file 'wt_1_f.fastq.gz' has an A nucleotide in the 5th base of the read?

#Obtain the specified file.
fastqFilePath <- system.file("reads", "wt_1_f.fastq.gz", package = "yeastRNASeq")

#Extract the sequences of the reads in the file.
library(ShortRead)
library(Biostrings)
reads <- readFastq(fastqFilePath)
reads_seq <- sread(reads)

#Calculate the statistics.
cM_all <- consensusMatrix(reads_seq, as.prob=TRUE, baseOnly=FALSE)
cM_bases <- consensusMatrix(reads_seq, as.prob=TRUE, baseOnly=TRUE)
cM_all['A', 5]
cM_bases['A', 5]

#Q2. What is the average numeric quality value of the 5th base of these reads?
mean(as(quality(reads), "matrix")[,5])

#The leeBamViews experiment data package contains aligned BAM files from an RNA seq experiment in yeast (the same experiment as in Questions 1 and 2, but that is not pertinent to the question).
#These reads are short reads (36bp) and have been aligned to the genome using a standard aligner, so potential junctions have been ignored (this makes some sense as yeast has very few junctions and the reads are very short).
library(leeBamViews)

#We will consider one of these files.
bamFilePath <- system.file("bam", "isowt5_13e.bam", package="leeBamViews")

#We will focus on the interval from 800,000 to 801,000 on yeast chromosome 13.
#A read duplicated by position is a read where at least one more read shares the same position.
#Q3. In this interval, how many reads are duplicated by position?

#Obtain the aligned reads in the file.
library(Rsamtools)
reads_aligned <- BamFile(bamFilePath)

#Retain the aligned reads in the interval.
interval <- GRanges(seqnames = "Scchr13", ranges = IRanges(start = c(800000), end = c(801000)))
parameters <- ScanBamParam(which = interval, what = scanBamWhat())
reads_alignedANDscanned <- scanBam(reads_aligned, param = parameters)

#Count the duplicated reads.
sum(table(reads_alignedANDscanned[[1]]$pos))-sum(table(reads_alignedANDscanned[[1]]$pos)==1)

#The leeBamViews experiment data package contains 8 BAM files in total, representing 8 different samples from 4 groups.
#Obtain the files.
bpaths <- list.files(system.file("bam", package="leeBamViews"), pattern = "bam$", full=TRUE)
reads_aligned_8 <- BamViews(bpaths)

#An objective of the original paper was the discovery of novel transcribed regions in yeast. One such region is Scchr13:807762-808068.
interval_noveltr <- GRanges(seqnames="Scchr13", ranges=IRanges(start = c(807762), end = c(808068)))

#Q4. What is the average number of reads across the 8 samples falling in this interval?

#Retain the reads in the interval containing the novel transcribed region.
bamRanges(reads_aligned_8) <- interval_noveltr
reads_alignedANDscanned_8 <- scanBam(reads_aligned_8)
dummy <- unlist(reads_alignedANDscanned_8)

#Calculate the mean number of reads across the samples in this interval.
(length(dummy[[1]]$pos)+length(dummy[[2]]$pos)+length(dummy[[3]]$pos)+length(dummy[[4]]$pos)+length(dummy[[5]]$pos)+length(dummy[[6]]$pos)+length(dummy[[7]]$pos)+length(dummy[[8]]$pos))/8

#In the lecture on the oligo package an ExpressionSet with 18 samples is constructed, representing normalized data from an Affymetrix gene expression microarray. The samples are divided into two groups given by the 'group' variable.
#Q5. What is the average expression across samples in the control group for the "8149273" probeset (this is a character identifier, not a row number)?

#Obtain the dataset.
library(oligo)
library(GEOquery)
getGEOSuppFiles("GSE38792")
untar("GSE38792/GSE38792_RAW.tar", exdir = "GSE38792/CEL")
filenames <- list.files("GSE38792/CEL",full=TRUE)
Affymetrix_GEmarray <- read.celfiles(filenames)

#Split the samples into the control group and the test group.
samplenames <- sampleNames(Affymetrix_GEmarray)
pData(Affymetrix_GEmarray)$sample <- samplenames
samplenames <- sub(".*_", "", samplenames)
samplenames <- sub(".CEL.gz$", "", samplenames)
sampleNames(Affymetrix_GEmarray) <- samplenames
pData(Affymetrix_GEmarray)$group <- ifelse(grepl("^OSA", sampleNames(Affymetrix_GEmarray)), "test", "control")

#Normalise the expression levels.
Affymetrix_GEmarray_normalised <- rma(Affymetrix_GEmarray)

#Extract the expression levels for the 8149273 probeset.
probeset <- match("8149273", rownames(Affymetrix_GEmarray_normalised))

#Average the expression levels for the 8149273 probeset across the control samples.
mean(exprs(Affymetrix_GEmarray_normalised[probeset,])[1:8])

#Use the limma package to fit a two group comparison between the control group and the OSA group, and borrow strength across the genes using eBayes(). Include all 18 samples in the model fit.
library(limma)

#Set up the design matrix.
Affymetrix_GEmarray_normalised$group <- factor(Affymetrix_GEmarray_normalised$group)
designmatrix <- model.matrix(~Affymetrix_GEmarray_normalised$group)

#Fit the two group comparison and use eBayes().
twogcomparison <- lmFit(Affymetrix_GEmarray_normalised, designmatrix)
twogcomparison <- eBayes(twogcomparison)

#Q6. What is the absolute value of the log fold change of the gene with the lowest P value?
abs(topTable(twogcomparison)$logFC[1])

#Q7. How many genes are differentially expressed between the two groups at an adj.P.value cutoff of 0.05?
twogcomparison_toptable <- topTable(twogcomparison)
ans <- subset(twogcomparison_toptable, adj.P.Val < 0.05)

#An example 450k dataset is contained in the minfiData package. This dataset contains 6 samples; 3 cancer and 3 normals. Cancer has been shown to be globally hypo-methylated (less methylated) compared to normal tissue of the same kind.
#Take the RGsetEx dataset in this package and preprocess it with the preprocessFunnorm function. For each sample, compute the average Beta value (percent methylation) across so-called OpenSea loci.

#Obtain and transform the data.
library(minfi)
library(minfiData)
RGsetEx_transformed <- preprocessFunnorm(RGsetEx)

#Obtain the beta values in the data.
RGsetEx_beta <- getBeta(RGsetEx_transformed)
beta_control <- RGsetEx_beta[, c(1,2,5)]
beta_cancer <- RGsetEx_beta[, c(3,4,6)]

#Retain the beta values in the OpenSea loci.
RGsetEx_is <- getIslandStatus(RGsetEx_transformed)
OpenSea_beta_control <- beta_control[RGsetEx_is == 'OpenSea', ]
OpenSea_beta_cancer <- beta_cancer[RGsetEx_is == 'OpenSea', ]

#Q8. What is the mean difference in beta values between the 3 normal samples and the 3 cancer samples, across OpenSea CpGs?
mean(OpenSea_beta_control) - mean(OpenSea_beta_cancer)

#The Caco2 cell line is a colon cancer cell line profiled by ENCODE. Obtain the narrowPeak DNase hyper sensitive sites computed by the analysis working group (AWG).
library(AnnotationHub)
ah <- AnnotationHub()
ah <- subset(ah, species=="Homo sapiens")
ah_Caco2 <- query(ah, c("Caco2", "AWG"))
narrowPeak_DNaseHypersensitive <- ah_Caco2[["AH22442"]]

#Q9. How many of these DNase hypersensitive sites contain one or more CpGs on the 450k array?
CpG <- granges(RGsetEx_transformed)
reduce(CpG)
unique(findOverlaps(CpG, narrowPeak_DNaseHypersensitive, type="within"))

#The zebrafishRNASeq package contains summarized data from an RNA-seq experiment in zebrafish in the form of a dataframe called zfGenes. The experiment compared 3 control samples to 3 treatment samples.
#Obtain the data.
library(DESeq2)
library(zebrafishRNASeq)
data("zfGenes")

#Each row is a transcript; the dataframe contains 92 rows with spikein transcripts; these have a rowname starting with "ERCC". Exclude these rows from the analysis.
zfGenes_minusERCC <- zfGenes[grep("^ERCC", rownames(zfGenes), invert = T), ]

#Use DESeq2 to perform a differential expression analysis between control and treatment. Do not discard (filter) genes and use the padj results output as the p-value.
zfGenes_minusERCC <- as.matrix(zfGenes_minusERCC)
labels <- DataFrame(sampleID = colnames(zfGenes_minusERCC), group = as.factor(c("control", "control", "control", "treatment", "treatment", "treatment")))
DEanalysis <- DESeqDataSetFromMatrix(zfGenes_minusERCC, labels, design = ~ group)
DEanalysis <- DESeq(DEanalysis)
DEanalysis_results <- results(DEanalysis)

#Q10. How many features are differentially expressed between control and treatment (padj <= 0.05)?
subset(DEanalysis_results, padj <= 0.05)