### This code analyzes the 12S full run data, including 1 full MiSeq run and 1 spike-in
### Inputs: ASV table, taxonomy, and phylogeny from qiime2
### Outputs: phyloseq objects (1. pre-processed data; 2. rarefied dataset), figures
### Last updated 6 Oct 2026 by Kara Andres (kjandr2@ilstu.edu)

### Clear work environment and load packages
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") 
library(ggplot2)
library(phyloseq) ## NOTE: must be installed BiocManager::install("phyloseq"); 'igraph' must be installed first
library(qiime2R) ## devtools::install_github("jbisanz/qiime2R")
library(reshape2)
library(decontam) ## BiocManager::install("decontam")
library(Biostrings)
library(stringr)
library(phylosmith) ## devtools::install_github("schuyler-smith/phylosmith")
library(RColorBrewer)
library(MiscMetabar)
library(colorspace)
library(vegan)
library(forcats)
library(microbiome) ## BiocManager::install("microbiome")
library(dplyr)

setwd("12S_full/")

### Produce phyloseq object
ASVs <- read_qza("12S_merged_table.qza")
taxonomy <- read_qza("taxonomy.qza")
metadata <- read_q2metadata("meta_data_R.txt")
rownames(metadata) <- metadata$SampleID
metadata$SampleID <- as.factor(metadata$SampleID)
metadata$fish_ID <- as.factor(metadata$fish_ID)
metadata$species <- as.factor(metadata$species)
metadata$family <- as.factor(metadata$family)
metadata$sp_abb <- as.factor(metadata$sp_abb)
metadata$diet_cat <- as.factor(metadata$diet_cat)
metadata$habitat_cat <- as.factor(metadata$habitat_cat)
sampledata <- sample_data(metadata)
rownames(sampledata) <- rownames(metadata)
R_tree <- read_qza("rooted-tree.qza")
phyloseq_dat <- phyloseq(otu_table(ASVs$data,taxa_are_rows = TRUE),tax_table(as.matrix(parse_taxonomy(taxonomy$data))),sampledata,phy_tree(R_tree$data))
phyloseq_dat # 240 taxa and 260 samples
sum(sample_sums(phyloseq_dat)) # 18,861,623 reads
saveRDS(phyloseq_dat,file = "phyloseq_dat.rds")

### Remove ASVs that are bacteria or mammal (human, pig, cow) 
phyloseq_euks <- phyloseq_dat %>%
  subset_taxa((Kingdom != "d__Bacteria") &
      (is.na(Class) | Class != "Mammalia"))
# 236 taxa and 260 samples

### Inspect negative controls 
phyloseq_blank <- subset_samples(phyloseq_euks, diet_cat=="Blank")
phyloseq_blank <- prune_taxa(taxa_sums(phyloseq_blank) > 0, phyloseq_blank)
phyloseq_blank@otu_table # 1 ASV present in 1 negative control
rowSums(unique(phyloseq_blank@otu_table)) # 3 reads total 

### Identify and remove contaminants by decontam
sample_data(phyloseq_euks)$is.neg <- sample_data(phyloseq_euks)$species == "Blank"
contam <- isContaminant(phyloseq_euks, method="prevalence", neg="is.neg", threshold=0.1)
contaminants <- rownames(contam)[which(contam$contaminant)] # 0 contaminant sequences identified
phyloseq_decontam <- prune_taxa(!rownames(phyloseq_euks@otu_table) %in% contaminants, phyloseq_euks)

### Remove negative controls
phyloseq_clean <- subset_samples(phyloseq_decontam, diet_cat!="Blank")

### Remove species with n=1 (Ameiurus natalis) 
phyloseq_clean <- subset_samples(phyloseq_clean, species!="Ameiurus_natalis")

### Merge resequencing data with initial run
sample_data(phyloseq_clean)$fish_ID <- as.factor(sample_data(phyloseq_clean)$fish_ID) # merge replicates
phyloseq_clean.nonrep <- merge_samples(phyloseq_clean,"fish_ID")

### Metadata coerced to NA -- manually add back sample data
sample_table <- unique(sample_data(phyloseq_clean)[,c("fish_ID","family","species","sp_abb","diet_cat","habitat_cat")])
rownames(sample_table) <- sample_table$fish_ID
sample_data(phyloseq_clean.nonrep) <- sample_table

### Remove ASVs belonging to host Genus
otu <- as(otu_table(phyloseq_clean.nonrep), "matrix")
tax <- as(tax_table(phyloseq_clean.nonrep), "matrix")
meta <- sample_data(phyloseq_clean.nonrep)
# Loop over samples, remove host ASVs
for (sample in sample_names(phyloseq_clean.nonrep)) {
  host_species <- as.character(meta[sample, "species"]$species)
  host_genus <- unlist(strsplit(host_species, "_"))[1]
  host_asvs <- rownames(tax)[!is.na(tax[, "Genus"]) & tax[, "Genus"] == host_genus] # Identify ASVs matching the host genus
  if (length(host_asvs) > 0) {
    otu[sample,host_asvs] <- 0 # Set host ASV counts to zero in this sample only
  }
}

# Rebuild new phyloseq object
otu_new <- otu_table(otu, taxa_are_rows = FALSE)
phyloseq_clean.nonrep.nohost <- phyloseq(otu_new, tax_table(phyloseq_clean.nonrep), sample_data(phyloseq_clean.nonrep))
phyloseq_nonzero <- prune_taxa(taxa_sums(phyloseq_clean.nonrep.nohost) > 0, phyloseq_clean.nonrep.nohost) # Remove ASVs that are now zero everywhere
phyloseq_nonzero <- prune_samples(sample_sums(phyloseq_nonzero) > 0, phyloseq_nonzero) # Remove samples that are now zero everywhere
phyloseq_nonzero # 136 taxa and 152 samples 

### Inspect sequencing depth per sample
hist(sample_sums(phyloseq_nonzero@otu_table), xlab="Number of reads", main="Histogram of number of reads per sample", breaks=50)
read_depth <- data.frame(fish_ID=sample_data(phyloseq_nonzero)$fish_ID,
                         family=sample_data(phyloseq_nonzero)$family,
                         species=sample_data(phyloseq_nonzero)$species,
                         sp_abb=sample_data(phyloseq_nonzero)$sp_abb,
                         read_depth=sample_sums(phyloseq_nonzero))
p <- ggplot(read_depth, aes(x=species, y=read_depth, fill=species)) + 
  geom_boxplot(outlier.shape=NA, outlier.colour=NA) + 
  geom_point(shape=16, size=2, color="gray") + 
  geom_hline(yintercept=1000, linetype="dashed") +
  labs(x="Species", y="Read Depth") + 
  theme_bw() +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        legend.position = "none")
p
ggsave("Figures/read_depth_all_samples.pdf", plot=p)

### Determine rarefaction depth with ASV accumulation curves by species
accum_curves <- accu_plot(phyloseq_nonzero, "species", add_nb_seq=TRUE, step=1000) + 
  # facet_wrap(~fact) + 
  theme_bw()
accum_curves
ggsave("Figures/12S_accum_curves.pdf", plot=accum_curves)

### Most samples plateau by 1,000 reads -- target rarefaction depth ###

### Identify samples with reads < 1,000
sample_sums(prune_samples(sample_sums(phyloseq_nonzero)<1000, phyloseq_nonzero))
length(sort(sample_sums(prune_samples(sample_sums(phyloseq_nonzero)<1000, phyloseq_nonzero)))) # 117 samples below 1,000 reads
SampleFilter1 <- names(sample_sums(phyloseq_nonzero))[sample_sums(phyloseq_nonzero) < 1000]
phyloseq.1000 <- subset_samples(phyloseq_nonzero,!sample_names(phyloseq_nonzero) %in% SampleFilter1)

### Produce final phyloseq after rarefaction
set.seed(1000)
phyloseq_rarefied <- rarefy_even_depth(phyloseq.1000, sample.size=min(sample_sums(phyloseq.1000)))
phyloseq_rarefied <- prune_taxa(taxa_sums(phyloseq_rarefied)>0,phyloseq_rarefied)
phyloseq_rarefied # 82 taxa and 35 samples
saveRDS(phyloseq_rarefied,file = "phyloseq_rarefied.rds")
sum(sample_sums(phyloseq_rarefied)) # 60,760 reads

### Sample size per species 
sample_size <- as.data.frame(with(sample_data(phyloseq_nonzero), table(species)))
clean_n <- as.data.frame(with(sample_data(phyloseq_rarefied), table(species)))
merge(sample_size, clean_n, by="species", all=TRUE)

### Transform final dataset to relative abundance
phyloseq_ReAb <- transform_sample_counts(phyloseq_rarefied, function(x) x / sum(x))

### Pie chart of Family: RelAb
ps_family <- tax_glom(phyloseq_rarefied, taxrank = "Family")
my.cols <- c(brewer.pal(8,"Dark2"), "#C46254","bisque","darkorange", rainbow_hcl(12), rainbow(13), "azure3")
data <- ps_family %>%
  psmelt() %>%
  as_tibble() %>%
  group_by(Family) %>%
  summarise(Reads = sum(Abundance), .groups = "drop") %>%
  mutate(RelAb = Reads / sum(Reads)) %>%
  mutate(Family = ifelse(RelAb < 0.01, "Other", Family)) %>% # Combine families < 1% into "Other"
  group_by(Family) %>%
  summarise(RelAb = sum(RelAb), .groups = "drop") %>%
  arrange(desc(RelAb)) %>%
  mutate(Family = factor(Family, levels = Family))
pie <- ggplot(data, aes(x = "", y = RelAb, fill = Family)) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Family", values = my.cols) +
  coord_polar("y") +
  theme_void() +
  theme(legend.title = element_text(size = 18),legend.text = element_text(size = 16),
  legend.box.spacing = unit(0, "cm"))
pie
ggsave("Figures/pie_chart_family.pdf", plot=pie)

### Closer look at taxonomic proportions (Family) per species
familyabundance <- phyloseq_rarefied %>%
  tax_glom(taxrank = "Family") %>%
  prune_samples(sample_sums(.) > 0, .) %>%  # Remove samples with zero Family-assigned reads
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>% # Melt to long format
  arrange(Family)
familyabundance_avg <- familyabundance %>%
  select(Family, species, sp_abb, diet_cat, habitat_cat, Abundance) %>%
  group_by(Family, species, sp_abb, diet_cat, habitat_cat) %>%
  summarize(avg_abundance = mean(Abundance), .groups = "drop") %>%
  mutate(Family = ifelse(avg_abundance < 0.01, "Other", Family)) %>% # Combine families < 1% into "Other"
  group_by(Family, species, sp_abb, diet_cat, habitat_cat) %>% # Sum families that were combined into "Other"
  summarize(avg_abundance = sum(avg_abundance), .groups = "drop")
p <- ggplot(familyabundance_avg) +
  geom_col(aes(x=species, y=avg_abundance, fill=fct_reorder(Family, avg_abundance, .desc=TRUE)), show.legend=TRUE) +
  facet_grid(.~diet_cat, scales="free", space="free") +
  ylab("Proportion of diet sequences") +
  scale_fill_manual(name="Family", values=my.cols) +
  xlab(NULL) +
  theme_bw() +
  theme(axis.text.y.left=element_text(size=15), axis.text.x=element_text(size=10, angle=45, hjust=1),
    axis.title.y=element_text(size=15))
p
ggsave("Figures/family_prop_by_species.pdf", plot=p)

### Closer look at taxonomic proportions (Genus) per species
genusabundance <- phyloseq_rarefied %>%
  tax_glom(taxrank="Genus") %>%
  prune_samples(sample_sums(.) > 0, .) %>%  # Remove samples with zero Genus-assigned reads
  transform_sample_counts(function(x) x / sum(x)) %>%
  psmelt() %>%
  arrange(Genus)
genusabundance_avg <- genusabundance %>%
  select(Genus, species, sp_abb, diet_cat, habitat_cat, Abundance) %>%
  group_by(Genus, species, sp_abb, diet_cat, habitat_cat) %>%
  summarize(avg_abundance = mean(Abundance), .groups = "drop") %>%
  mutate(Genus=ifelse(avg_abundance < 0.01, "Other", Genus)) %>% # Combine families < 1% into "Other"
  group_by(Genus, species, sp_abb, diet_cat, habitat_cat) %>% # Sum families that were combined into "Other"
  summarize(avg_abundance = sum(avg_abundance), .groups = "drop")
p <- ggplot(genusabundance_avg) +
  geom_col(aes(x=species, y=avg_abundance, fill=fct_reorder(Genus, avg_abundance, .desc=TRUE)), show.legend=TRUE) +
  facet_grid(.~diet_cat, scales="free", space="free") +
  ylab("Proportion of diet sequences") +
  scale_fill_manual(name="Genus", values=my.cols) +
  xlab(NULL) +
  theme_bw() +
  theme(axis.text.y.left=element_text(size=15), axis.text.x=element_text(size=10, angle=45, hjust=1),
        axis.title.y=element_text(size=15))
p
ggsave("Figures/genus_prop_by_species.pdf", plot=p)

### Shannon diversity per species
my.cols <- c(brewer.pal(8,"Dark2"), rainbow_hcl(10))
plot_richness(phyloseq_rarefied, x="species", measures="shannon", color="species") +
  geom_boxplot(alpha=0.6) + 
  scale_color_manual(values=my.cols) +
  theme_bw() + ylab("Shannon Diversity") +
  theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=10))

### Diversity stats per sample
alpha_ASV_meta <- data.frame(estimate_richness(physeq = phyloseq_rarefied, measures = c("Shannon", "Observed")),
                             Sample_ID=sample_data(phyloseq_rarefied)$fish_ID,
                             fish_species=sample_data(phyloseq_rarefied)$species,
                             diet_cat=sample_data(phyloseq_rarefied)$diet_cat,
                             habitat_cat=sample_data(phyloseq_rarefied)$habitat_cat)

### ANOVA to test for differences in Shannon diversity
richness.aov <- aov(Shannon ~ fish_species, data=alpha_ASV_meta)
summary(richness.aov)

### PCOA of species
ord.PCoA.jaccard <- ordinate(phyloseq_ReAb, method="PCoA", distance="jaccard")
pcoa.cols <- c(brewer.pal(8,"Dark2"), rainbow_hcl(10))
plot_ordination(phyloseq_ReAb, ord.PCoA.jaccard, type="samples", color="species") +
  stat_ellipse(type="norm") +
  geom_point(size=3, stroke=1) + 
  scale_color_manual(name="Host Species", values=pcoa.cols) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20), legend.title=element_text(size=20), legend.text=element_text(size=14))

### PERMANOVA
meta <- microbiome::meta(phyloseq_ReAb)
abund <- microbiome::abundances(phyloseq_ReAb)
dist_jaccard <- phyloseq::distance(phyloseq_ReAb, method = "jaccard")
adonis2(formula = dist_jaccard ~ species, data=meta, permutations=999, method="jaccard")

### PERMDISP
dist_adund <- phyloseq::distance(phyloseq_ReAb, method="jaccard")
beta <- betadisper(dist_adund, meta$species)
permutest(beta)

### Alpha diversity values for comparison to gut microbiome
alpha_div <- estimate_richness(physeq = phyloseq_rarefied, measures = c("Shannon", "Observed"))
metadata <- sample_data(object = phyloseq_rarefied) %>% 
  data.frame(.) # convert to a data.frame
alpha_div_metadata <- cbind(alpha_div, metadata)
write.csv(alpha_div_metadata, "12S_alpha_div_metadata.csv")