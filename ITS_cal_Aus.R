##  ITS California Australia ##
library(dada2); packageVersion("dada2")
pathF <- "/Users/liaolab/Documents/koko_test/SO_origional_data_Ko_ITS/cut_R1/trimmo" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
fns <- sort(list.files(pathF, pattern="trimo.fastq.gz"))
out<-filterAndTrim(file.path(pathF,fns), file.path(filtpathF,fns), 
                   maxEE=3, truncQ=12, minLen=50,
                   compress=TRUE, verbose=TRUE, multithread=TRUE)

filtpath <- "/Users/liaolab/Documents/koko_test/SO_origional_data_Ko_ITS/cut_R1/trimmo/filtered" # CHANGE ME to the directory containing your filtered fastq files
filts <- list.files(filtpath, pattern="fastq.gz", full.names=TRUE) # CHANGE if different file extensions
sample.names <- sapply(strsplit(basename(filts), "_"), `[`, 1) # Assumes filename = E40_S1_L001_R1_001_trimo.fastq.gz
names(filts) <- sample.names
# Learn error rates
set.seed(100)
err <- learnErrors(filts, nbases = 1e8, multithread=TRUE, randomize=TRUE)
# Infer sequence variants
dds <- vector("list", length(sample.names))
names(dds) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filts[[sam]])
  dds[[sam]] <- dada(derep, err=err, multithread=TRUE)
}

# Construct sequence table and write to disk
seqtab <- makeSequenceTable(dds)
saveRDS(seqtab, "/Users/liaolab/Documents/koko_test/asilomar_redo/ITS_calaus/seqtab.rds")

dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)
saveRDS(seqtab.nochim, "/Users/liaolab/Documents/koko_test/asilomar_redo/ITS_calaus/seqtabnochim.rds")

getN <- function(x) sum(getUniques(x))


track <- cbind(out, rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "nonchim")
rownames(track) <- sample.names
track
write.table(track,"/Users/liaolab/Documents/koko_test/asilomar_redo/ITS_calaus/track_single.txt",quote=F, sep="\t")



tax <- assignTaxonomy(seqtab, "/Users/liaolab/Documents/koko_test/asilomar_redo/ITS_calaus/sh_general_release_dynamic_s_02.02.2019.fasta", multithread=TRUE)
saveRDS(tax, "/Users/liaolab/Documents/koko_test/asilomar_redo/ITS_calaus/tax_assign.rds")

tax.print <- tax # Removing sequence rownames for display only
rownames(tax.print) <- NULL
head(tax.print)