
library("dada2"); packageVersion("dada2")
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")


path <- "/Users/nqhuynh/Documents/rdna"
list.files(path)
##Filter and Trim
# Forward and reverse fastq filenames have format: SAMPLENAME_L001_R1.fastq and SAMPLENAME_L001_R2.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- make.unique(sapply(strsplit(basename(fnFs), "_"), `[`, 1))
sample.names
sample.namesfull <- sapply(strsplit(basename(fnFs), "_16S"), `[`, 1)

##Examine quality profiles of forward and reverse reads
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

##Perform filtering and trimming
#Assign the filenames for the filtered fastq.gz files.
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#Filter the forward and reverse reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220, 220),
                     maxN=0, truncQ=2, maxEE = c(2,2), rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
head(out)



##Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

##Dereplication
#Dereplicate the filtered fastq files
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

##Sample Inference
#Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#Inspecting the dada-class object returned by dada
dadaFs[[1]]
dadaRs[[1]]

##Merge paired reads
#Merge the denoised forward and reverse reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

##Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

##Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab) # 1= 0% seq from chimeras removed, 0,9=1% seq from chimeras removed...etc

##Track reads through the pipeline
getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

##Assign taxonomy
taxaRC <- assignTaxonomy(seqtab.nochim, "/Users/nqhuynh/Documents/Trainset/silva_nr_v123_train_set.fa.gz", tryRC=TRUE)
taxaSp <- addSpecies(taxaRC, "/Users/nqhuynh/Documents/Trainset/silva_species_assignment_v123.fa.gz")


#taxonomic assignments
taxa.print <- taxaSp # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


#Environmental data---------------------------
#prepare it following the script in https://github.com/MicroB3-IS/osd-analysis/wiki/OSD-2014-environmental-data-csv-documentation


envdata <- read.table("/Users/nqhuynh/Documents/rdna/3/data1.txt", header=TRUE, fill=TRUE, row.names=1, sep="\t", dec=",")
envdata
###Mapping sample names
env_Full_names <- rownames(envdata1)
env_Short_names <- make.unique(sapply(strsplit(rownames(envdata1), "_"), `[`,1))
env_Names <- cbind(Short_Names = env_Short_names, env_Full_names)

Sequencing_names <- cbind(Short_Names = sample.names, S_Full_Names = sample.namesfull)
Name_list <- merge(Sequencing_names, env_Names, all = T, all.x = TRUE)
View(Name_list)
write.table(Name_list, "/Users/nqhuynh/Documents/Name_list.txt", sep="\t")



####
rownames(envdata) <- make.unique(sapply(strsplit(rownames(envdata), "_"), `[`, 1))
rownames(envdata)
DAT <- sample_data(envdata)
DAT

#fichier Physeq
ps123 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxaSp))
ps123 <- merge_phyloseq(ps123, DAT)
ps123

save(ps1, "OSD_phyloseq.Rdata")
#change nom seq/OTU
taxa_names(ps1)
taxa_names(ps1) <- paste0("OTU", seq(ntaxa(ps1)))


# Extract abundance matrix from the phyloseq object
OTU1 = as(otu_table(ps1), "matrix")
# transpose if necessary
if(taxa_are_rows(ps1)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)

# Extract abundance matrix from the phyloseq object
TAX = as(tax_table(ps1), "matrix")
# transpose if necessary
if(taxa_are_rows(ps1)){TAX <- t(TAX)}
# Coerce to data.frame
TAXdf = as.data.frame(TAX)

# Extract abundance matrix from the phyloseq object
DAT1 = as(sample_data(ps1), "matrix")
# transpose if necessary
if(taxa_are_rows(ps1)){DAT1 <- t(DAT1)}
# Coerce to data.frame
DAT1df = as.data.frame(DAT1)
write.table(OTU1, "/Users/nqhuynh/Documents/1_OTU1.txt", sep="\t")
write.table(TAX, "/Users/nqhuynh/Documents/1_TAX.txt", sep="\t")
write.table(DAT1, "/Users/nqhuynh/Documents/1_DAT.txt", sep="\t")


save(ps123, file = "OSD_phyloseq_123.Rdata")
