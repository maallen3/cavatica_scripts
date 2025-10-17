#install.packages("renv")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

renv::init(force=TRUE, bioconductor = TRUE)
#renv::install("readr")
#renv::install("readxl")
#renv::install("GenomicFeatures")
#renv::install("DESeq2")

#renv::install("BiocVersion")
#renv::install("tximportData")
#renv::install("GenomicFeatures")
#renv::install("tximport")
#renv::install(dplyr)
#renv::install("txdbmaker")
#renv::install("rrcov")

library(GenomicFeatures)
library(tximport)
library(dplyr)
library(txdbmaker)
library(DESeq2)
library(readr)
library(readxl)

mainindir="/sbgenomics/project-files/"
countdir="/sbgenomics/project-files/kallisto/"
metadatadir="/sbgenomics/project-files/metadata/"
outdir="/sbgenomics/output-files/"

biospecimans <- paste(metadatadir,"include_biospecimenData_20241028T223221Z.xlsx",sep="")
filemanifest <- paste(metadatadir,"include_familyManifest_20241029T155031Z.tsv",sep="")
kallistolinenumbers <- paste(mainindir,"kallisto_file_sizes.csv", sep="")

#this is the number of lines in the Kalistofile. You can get the file for this by running Kallisto_file_linenumbers
#numberoflines=200401
numberoflines=244998

#https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#kallisto_with_TSV_files
#ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.primary_assembly.annotation.gtf.gz
#for 200401 line kallisto files use gencode.v27.annotation.gtf
#gtf = paste(mainindir,"gencode.v27.annotation.gtf", sep="") 
#for 244998 line kallisto files use gencode.v39.primary_assembly.annotation.gtf
gtf = paste(mainindir,"gencode.v39.primary_assembly.annotation.gtf", sep="") # for 244998 line kallisto files


#read in the metadata
sampleTable <- read.csv(filemanifest, sep="\t")
biospecimendf <- read_excel(biospecimans, sheet = "Biospecimens")
biospecimendf <- biospecimendf %>%
  rename_with(~ gsub(" ", ".", .))
dflinenumbers <-read.csv(kallistolinenumbers)


unique(biospecimendf$Study.Code)
unique(biospecimendf$Collection.Sample.Type)
unique(biospecimendf$External.Sample.ID)

#filter for HTP and not "White" since that's a separate type of RNA-seq
#thisstudy="HTP"
thisstudy="X01-Hakonarson"
onestudy <- sampleTable %>% right_join(biospecimendf, by=c("Participant.ID", "External.Sample.ID","Sample.ID" ),relationship ="many-to-many") 
onestudy<- onestudy %>% filter(Study.Code==thisstudy) %>%
  filter(!grepl("White", External.Sample.ID))

#filter for the correct number of lines


dfonefilesize <- dflinenumbers %>% filter(rows==numberoflines)
onestudy <- onestudy %>% filter(File.Name %in% dfonefilesize$filename)

#now make a list of files
files <- file.path(countdir, onestudy$File.Name)
names(files) <- onestudy$File.Name

length(files)

#set up the name conversion from Transcript name to gene name
gencode_txdb <- makeTxDbFromGFF(gtf)
k <- keys(gencode_txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(gencode_txdb, k, "GENEID", "TXNAME")

#read into a txiimport file for Deseq2
txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)


head(txi.kallisto.tsv$counts)


sampleTable<- sampleTable %>% select(File.Name, Participant.ID, Down.Syndrome.Status) %>% filter(File.Name %in% colnames(txi.kallisto.tsv$counts))

# Get the column names of the count matrix (sample names)
count_colnames <- colnames(txi.kallisto.tsv$counts)

# Reorder sampleTable to match the order of the count columns
sampleTable <- sampleTable[match(count_colnames, sampleTable$File.Name), ]


dds <- DESeqDataSetFromTximport(txi=txi.kallisto.tsv, sampleTable, ~Down.Syndrome.Status)

dds <- DESeq(dds)

normcounts <- counts(dds, normalized=TRUE)
rawcounts <- counts(dds, normalized=FALSE)


write.csv(normcounts,paste(outdir,"kallisto_to_DESeq_normcounts_",numberoflines,"_",thisstudy,"_notwhite",".csv", sep=""))
write.csv(rawcounts,paste(outdir,"kallisto_to_DESeq_rawcounts_",numberoflines,"_",thisstudy,"_notwhite",".csv", sep=""))
