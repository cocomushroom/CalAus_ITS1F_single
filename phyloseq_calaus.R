setwd("/Users/kc178/Documents/Florida_projects/CSP504376_AusCal/SO_analysis")
library("ape")
library("vegan")
library("tidyr")
library("metacoder")
library(phyloseq)
library("ggplot2")

## Load taxa produced by dada2##
tax30<-readRDS(file = "/Users/kc178/Documents/Florida_projects/CSP504376_AusCal/SO_analysis/tax30_assign.rds")
tax30_R<-gsub(".__", "",tax30$tax) ## If having problem with tax rank, skip this



## Load seqtab.nochim produced by dada2##
seqtab.nochim<-readRDS(file = "/Users/kc178/Documents/Florida_projects/CSP504376_AusCal/SO_analysis/seqtabnochim.rds")

## Not a QIIME2 format
meta <-read.delim("/Users/kc178/Documents/Florida_projects/CSP504376_AusCal/SO_analysis/metadata_prep/meta_SO.txt",row.names = 1)

##### Start Phyloseq#####
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(meta), 
               tax_table(tax30_R))
ps


#If want to double check taxonomy assignment
#library(microbiome)
#write_phyloseq(ps, 'TAXONOMY')

#If want to double check taxonomy assignment
#library(microbiome)
#write_phyloseq(pseq, 'TAXONOMY')

# Keep fungi only
psf <- subset_taxa(ps, Kingdom=="Fungi")
saveRDS(psf, "/Users/kc178/Documents/Florida_projects/CSP504376_AusCal/SO_analysis/psf.rds")


# check read distribution
readsumsdf = data.frame(nreads = sort(taxa_sums(psf), TRUE), sorted = 1:ntaxa(psf), 
                        type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(psf), 
                                                        TRUE), sorted = 1:nsamples(psf), type = "Samples"))
png("read_distribution.png", width=2000, height=2500,res=250)
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
dev.off()

nreads = sort(sample_sums(psf)) # if need the lowest read number for a sample to determine the rarification number
nreads

##GP.chl = prune_samples(sample_sums(GP.chl)>=20, GP.chl) ## might consider removing samples with low reads
psftrim = prune_samples(sample_sums(psf)>=4000, psf)

#### Create Rarified table ### 
set.seed(1111)
psf_trimR = rarefy_even_depth(psftrim, sample.size = nreads[[1]])
saveRDS(psf_trimR, "/Users/kc178/Documents/Florida_projects/CSP504376_AusCal/SO_analysis/psf_trimR.rds")


psf_trim_soil <- subset_samples(psftrim, Sample_type=="Soil")
psf_trim_root <- subset_samples(psftrim, Sample_type=="Root")


png("barplot_soilvsroot.png", width=3000, height=3000,res=350)
title = "all_soil root"
plot_bar(psf_trimR, "Sample_type", fill = "Order", title = title) + coord_flip() + 
  ylab("Percentage of Sequences")
dev.off()

#### 
psf_trimR_soil <- subset_samples(psf_trimR, Sample_type=="Soil")
psf_trimR_root <- subset_samples(psf_trimR, Sample_type=="Root")


### Barplot ###
png("soil_country.png", width=3000, height=3000,res=350)
title = "soil_country"
plot_bar(psf_trimR_soil, "Country", fill = "Order", title = title) + coord_flip() + 
  ylab("Percentage of Sequences")
dev.off()

png("root_country.png", width=3000, height=3000,res=350)
title = "root_country"
plot_bar(psf_trimR_root, "Country", fill = "Order", title = title) + coord_flip() + 
  ylab("Percentage of Sequences")
dev.off()


## PCoA ##
png("MDS_all.png", width=3000, height=3000,res=350)
logt  = transform_sample_counts(psf, function(x) log(1 + x) )
out.pcoa.logt <- ordinate(logt, method = "MDS", distance = "bray")
evals <- out.pcoa.logt$values$Eigenvalues
plot_ordination(logt, out.pcoa.logt, type = "samples") + labs(col = "Slash pile number") +
  coord_fixed(sqrt(evals[2] / evals[1]))
dev.off()

png("Pcoa_all.png", width=3000, height=3000,res=350)
logt  = transform_sample_counts(psftrim, function(x) log(1 + x) )
out.pcoa.logt <- ordinate(logt, method = "PCoA", distance = "bray")
evals <- out.pcoa.logt$values$Eigenvalues
plot_ordination(logt, out.pcoa.logt, type = "samples", 
                color = "Forest_type", shape = "Sample_type") + labs(col = "Slash pile number") +
  coord_fixed(sqrt(evals[2] / evals[1]))
dev.off()



png("Pcoa_soil.png", width=3000, height=3000,res=350)
logt  = transform_sample_counts(psf_trim_soil, function(x) log(1 + x) )
out.pcoa.logt <- ordinate(logt, method = "PCoA", distance = "bray")
evals <- out.pcoa.logt$values$Eigenvalues
plot_ordination(logt, out.pcoa.logt, type = "samples", 
                color = "Forest_type", shape = "Country") + labs(col = "Slash pile number") +
  coord_fixed(sqrt(evals[2] / evals[1]))
dev.off()


png("Pcoa_root.png", width=3000, height=3000,res=350)
title = "PCoA root"
psf_trimR_root_ord = ordinate(psf_trimR_root, "PCoA", "bray")
p = plot_ordination(psf_trimR_root, psf_trimR_root_ord, color = "Forest_type", shape = "Country")
p = p + geom_point(size = 6, alpha = 0.7) + ggtitle(title)
p
dev.off()


png("network.png", width=2000, height=2500,res=250)
ig = make_network(psf, type = "samples", distance = "bray", max.dist = 0.95)
plot_network(ig, psf, color = "Forest_type", shape = "Sample_type", line_weight = 0.4, 
             label = NULL)
dev.off()


plot_richness(psf, x = "Forest_type") + geom_boxplot()


