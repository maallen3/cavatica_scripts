#install.packages("renv")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

renv::init(force=TRUE, bioconductor = TRUE)
#renv::install("GenomicFeatures")
#renv::install("DESeq2")

#renv::install("BiocVersion")
#renv::install("tximportData")
#renv::install("GenomicFeatures")
#renv::install("tximport")
#renv::install(dplyr)
#renv::install("txdbmaker")

library(GenomicFeatures)
library(tximport)
library(dplyr)
library(txdbmaker)
library(DESeq2)

mainindir="/sbgenomics/project-files/"
countdir="/sbgenomics/project-files/kallisto/"
metadatadir="/sbgenomics/project-files/metadata/"
outdir="/sbgenomics/output-files/"


#numberoflines=200401
numberoflines=244998
df <-read.csv(paste(mainindir,"kallisto_file_sizes.csv", sep=""))
df200401 <- df %>% filter(rows==numberoflines)
files <- file.path(countdir, df200401$filename)
names(files) <- df200401$filename

length(files)

#https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#kallisto_with_TSV_files
#ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.primary_assembly.annotation.gtf.gz
#gtf = paste(mainindir,"gencode.v27.annotation.gtf", sep="") #for 200401 line kallisto files
gtf = paste(mainindir,"gencode.v39.primary_assembly.annotation.gtf", sep="") # for 244998 line kallisto files
gencode_txdb <- makeTxDbFromGFF(gtf)

k <- keys(gencode_txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(gencode_txdb, k, "GENEID", "TXNAME")

txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)


head(txi.kallisto.tsv$counts)

sampleTable <- read.csv(paste(metadatadir,"include_familyManifest_20241029T155031Z.tsv", sep=""), sep="\t")

sampleTable<- sampleTable %>% select(File.Name, Participant.ID, Down.Syndrome.Status) %>% filter(File.Name %in% colnames(txi.kallisto.tsv$counts))

# Get the column names of the count matrix (sample names)
count_colnames <- colnames(txi.kallisto.tsv$counts)

# Reorder sampleTable to match the order of the count columns
sampleTable <- sampleTable[match(count_colnames, sampleTable$File.Name), ]


dds <- DESeqDataSetFromTximport(txi=txi.kallisto.tsv, sampleTable, ~Down.Syndrome.Status)

dds <- DESeq(dds)

normcounts <- counts(dds, normalized=TRUE)
rawcounts <- counts(dds, normalized=FALSE)

write.csv(normcounts,paste(outdir,"kallisto_to_DESeq_normcounts_",numberoflines,".csv", sep=""))
write.csv(rawcounts,paste(outdir,"kallisto_to_DESeq_rawcounts_",numberoflines,".csv", sep=""))
