### ENVIRONMENT ###

library("phyloseq")
library("tidyr")
library("ggplot2")
library("viridis")
library("DESeq2")
library("dendextend")
library("ggrepel")
library("vegan")

### READING IN DATA ###

otu_counts_tab <- read.table("combined_otu_counts_tab.txt", header=TRUE, row.names=1, sep="\t")
head(otu_counts_tab)
dim(otu_counts_tab)
colnames(otu_counts_tab)
head(rownames(otu_counts_tab))

taxonomy_tab <- read.table("combined_taxonomy_tab.txt", header=TRUE, row.names=1, sep="\t", na.strings="<NA>")
head(taxonomy_tab)     
dim(taxonomy_tab)
colnames(taxonomy_tab)

sample_info_tab <- read.table("sample_info_tab.txt", header=T, check.names=F, row.names=1, sep="\t")
head(sample_info_tab)
dim(sample_info_tab)
colnames(sample_info_tab)

### TAXONOMIC SUMMARIES ###

counts_phy <- otu_table(otu_counts_tab, taxa_are_rows=TRUE)
tax_phy <- tax_table(as.matrix(taxonomy_tab))
sample_phy <- sample_data(sample_info_tab)
aqua_phy <- phyloseq(counts_phy, tax_phy, sample_phy)







