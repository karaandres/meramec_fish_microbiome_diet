### This code analyzes the 18S full run data, including 2 full MiSeq runs and 1 spike-in
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
library(tidyr)

setwd("18S_full/")

### Produce phyloseq object
ASVs <- read_qza("18S_merged_table.qza")
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
phyloseq_dat # 7184 taxa and 501 samples
sum(sample_sums(phyloseq_dat)) # 17,024,349 reads
saveRDS(phyloseq_dat,file = "phyloseq_dat.rds")

### Remove ASVs that are bacteria or fish (fishes will be captured by 12S) 
phyloseq_euks <- phyloseq_dat %>%
  subset_taxa((Kingdom != "d__Bacteria") &
                (is.na(Class) | Class != "Actinopteri"))
# 7069 taxa and 501 samples

### Inspect negative controls 
phyloseq_blank <- subset_samples(phyloseq_euks, diet_cat=="Blank")
phyloseq_blank <- prune_taxa(taxa_sums(phyloseq_blank) > 0, phyloseq_blank)
phyloseq_blank@otu_table # 2 ASVs present in 2 negative controls
rowSums(unique(phyloseq_blank@otu_table)) # 6 reads total 

### Identify and remove contaminants by decontam
sample_data(phyloseq_euks)$is.neg <- sample_data(phyloseq_euks)$species == "Blank"
contam <- isContaminant(phyloseq_euks, method="prevalence", neg="is.neg", threshold=0.1)
contaminants <- rownames(contam)[which(contam$contaminant)] # 0 contaminant sequences identified
phyloseq_decontam <- prune_taxa(!rownames(phyloseq_euks@otu_table) %in% contaminants, phyloseq_euks)

### Remove negative controls
phyloseq_clean <- subset_samples(phyloseq_decontam, diet_cat!="Blank")

### Remove species with n=1 (Ameiurus natalis) 
phyloseq_clean <- subset_samples(phyloseq_clean, species!="Ameiurus_natalis")

### Inspect and remove samples with low reads (if one of the two replicates has no reads)
lowreads <- subset_samples(phyloseq_clean, sample_sums(phyloseq_clean) == 0)
SampleFilter <- lowreads@sam_data$fish_ID
phyloseq_nonzero <- subset_samples(phyloseq_clean,!fish_ID %in% SampleFilter)

### Transform to relative abundance and calculate correlation coefficient
phyloseq_ReAb <- transform_sample_counts(phyloseq_nonzero, function(OTU) OTU*100/sum(OTU))
samples_to_keep <- !grepl("reseq", sample_names(phyloseq_nonzero)) # remove any resequenced samples
phyloseq_ReAb_filtered <- prune_samples(samples_to_keep, phyloseq_nonzero)
Pearson_cor_table <- as.data.frame(cor(phyloseq_ReAb_filtered@otu_table), method="pearson") # calculate correlation coefficient
keep <- sort(colnames(Pearson_cor_table))
Pearson_cor_table_srt <- Pearson_cor_table[,keep]
Pearson_cor_table_srt <- Pearson_cor_table_srt[order(row.names(Pearson_cor_table_srt)), ]
df <- data.frame(Samplename=NA, Pearson_values=NA, Reads=NA)
number_index <- seq(1,452, by=2)
for(i in number_index){
  df[i,1] <- rownames(Pearson_cor_table_srt)[i]
  df[i,2] <- Pearson_cor_table_srt[i, i+1]
  df[i,3] <- colSums(phyloseq_clean@otu_table)[i]+colSums(phyloseq_clean@otu_table)[i+1]
}
pearson_table <- na.omit(df)
pearson_table$species <- str_split_fixed(pearson_table$Samplename, "_", 2)[,1]
mean(pearson_table$Pearson_values) # mean pearson correlation 0.987
quantile(pearson_table$Pearson_values) # 0%: 0.764, 25%:0.991, 50%: 0.997, 75%: 0.999, 100%: 0.999 
p <- ggplot(pearson_table, aes(x=species, y=Pearson_values, fill=species)) + 
  geom_boxplot(outlier.shape=NA, outlier.colour=NA) + 
  geom_point(shape=16, size=2, color="gray") + 
  labs(x="Species", y="Pearson correlation") + 
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20), legend.position="none")
p
ggsave("Figures/pearson_correlation_technical_reps.pdf", plot=p)

### Merge technical replicates and resequenced samples 
sample_data(phyloseq_clean)$fish_ID <- as.factor(sample_data(phyloseq_clean)$fish_ID) # merge replicates
phyloseq_clean.nonrep <- merge_samples(phyloseq_clean,"fish_ID")

### Metadata coerced to NA -- manually add back sample data
sample_table <- unique(sample_data(phyloseq_clean)[,c("fish_ID","family","species","sp_abb","diet_cat","habitat_cat")])
rownames(sample_table) <- sample_table$fish_ID
sample_data(phyloseq_clean.nonrep) <- sample_table

### Inspect sequencing depth per sample
hist(sample_sums(phyloseq_clean.nonrep@otu_table), xlab="Number of reads", main="Histogram of number of reads per sample", breaks=50)
read_depth <- data.frame(fish_ID=sample_data(phyloseq_clean.nonrep)$fish_ID,
                         family=sample_data(phyloseq_clean.nonrep)$family,
                         species=sample_data(phyloseq_clean.nonrep)$species,
                         sp_abb=sample_data(phyloseq_clean.nonrep)$sp_abb,
                         read_depth=sample_sums(phyloseq_clean.nonrep))
p <- ggplot(read_depth, aes(x=species, y=read_depth, fill=species)) + 
  geom_boxplot(outlier.shape=NA, outlier.colour=NA) + 
  geom_point(shape=16, size=2, color="gray") + 
  geom_hline(yintercept=5000, linetype="dashed") +
  labs(x="Species", y="Read Depth", size=12) + 
  theme_bw() +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        legend.position = "none")
p
ggsave("Figures/read_depth_all_samples.pdf", plot=p)

### Inspect and remove samples with no reads
lowreads <- subset_samples(phyloseq_clean.nonrep, sample_sums(phyloseq_clean.nonrep) == 0)
SampleFilter <- lowreads@sam_data$fish_ID
phyloseq_nonzero <- subset_samples(phyloseq_clean.nonrep,!fish_ID %in% SampleFilter)

### Determine rarefaction depth with ASV accumulation curves by species
accum_curves <- accu_plot(phyloseq_nonzero, "species", add_nb_seq=TRUE, step=1000) + 
  # facet_wrap(~fact) + 
  theme_bw()
accum_curves
ggsave("Figures/18S_accum_curves.pdf", plot=accum_curves)

### Most samples plateau by 5,000 reads -- target rarefaction depth ###

### Identify samples with reads < 5,000
sample_sums(phyloseq_nonzero@otu_table)
length(sort(sample_sums(prune_samples(sample_sums(phyloseq_nonzero)<5000, phyloseq_nonzero)))) # 31 samples below 5000 reads
SampleFilter1 <- names(sample_sums(phyloseq_nonzero))[sample_sums(phyloseq_nonzero) < 5000]
phyloseq.5000 <- subset_samples(phyloseq_nonzero,!sample_names(phyloseq_nonzero) %in% SampleFilter1)

### Produce final phyloseq after rarefaction
set.seed(5000)
phyloseq_rarefied <- rarefy_even_depth(phyloseq.5000, sample.size=min(sample_sums(phyloseq.5000)))
phyloseq_rarefied <- prune_taxa(taxa_sums(phyloseq_rarefied)>0,phyloseq_rarefied)
phyloseq_rarefied # 4878 taxa and 200 samples
saveRDS(phyloseq_rarefied,file = "phyloseq_rarefied.rds")
sum(sample_sums(phyloseq_rarefied)) # 1,032,000 reads

### Sample size per species 
sample_size <- as.data.frame(with(sample_data(phyloseq_clean.nonrep), table(species)))
clean_n <- as.data.frame(with(sample_data(phyloseq_rarefied), table(species)))
merge(sample_size, clean_n, by="species", all=TRUE)

### Transform final dataset to relative abundance
phyloseq_ReAb <- transform_sample_counts(phyloseq_rarefied, function(x) x/sum(x))

### Pie chart of Phyla: RelAb
ps_phylum <- tax_glom(phyloseq_rarefied, taxrank = "Phylum")
my.cols <- c(brewer.pal(8,"Dark2"), "#C46254","bisque","darkorange", rainbow_hcl(12), rainbow(13), "azure3")
data <- ps_phylum %>%
  psmelt() %>%
  as_tibble() %>%
  group_by(Phylum) %>%
  summarise(Reads = sum(Abundance), .groups = "drop") %>%
  mutate(RelAb = Reads / sum(Reads)) %>%
  mutate(Phylum = ifelse(RelAb < 0.01, "Other", Phylum)) %>% # Combine phyla < 1% into "Other"
  group_by(Phylum) %>%
  summarise(RelAb = sum(RelAb), .groups = "drop") %>%
  arrange(desc(RelAb)) %>%
  mutate(Phylum = factor(Phylum, levels = Phylum))
pie <- ggplot(data, aes(x = "", y = RelAb, fill = Phylum)) +
  geom_col(width = 1) +
  scale_fill_manual(name = "Phylum", values = my.cols) +
  coord_polar("y") +
  theme_void() +
  theme(legend.title = element_text(size = 18),legend.text = element_text(size = 16),
        legend.box.spacing = unit(0, "cm"))
pie
ggsave("Figures/pie_chart_phylum.pdf", plot=pie)

# Phylum-level breakdown: plot top phyla
names_unknown <- taxa_names(phyloseq_ReAb@tax_table)[is.na(phyloseq_ReAb@tax_table[,"Phylum"])]
phyloseq_ReAb@tax_table[is.na(phyloseq_ReAb@tax_table[,"Phylum"]),"Phylum"] <- c("unclassified Phyla")
phyloseq_ReAb.phylum <- tax_glom(phyloseq_ReAb, NArm=FALSE, taxrank="Phylum")
base.phylum <- phyloseq_ReAb.phylum %>% psmelt()
base.phylum$Phylum <- factor(base.phylum$Phylum, levels=unique(base.phylum$Phylum))
base.phylum$species <- factor(base.phylum$species)
base.phylum <- base.phylum[base.phylum$Abundance>0.1,]
unique(base.phylum$Phylum) # 34 unique phyla w/ abundance > 0.1
p <- ggplot(base.phylum, aes(x=fish_ID, y=Abundance, fill=Phylum)) + 
  geom_col(width=0.9) +
  scale_fill_manual(name="Phylum", values=my.cols) +
  scale_color_manual(name="Phylum", values=my.cols) +
  ylab("Relative Abundance (%)") + xlab("Sample ID") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text=element_text(size=8), axis.title=element_text(size=14),
        legend.key.size = unit(0.2, 'cm'),
        legend.text = element_text(size=8))
p

### Closer look at taxonomic proportions (Phyla) per species
phylumabundance <- phyloseq_rarefied %>%
  tax_glom(taxrank = "Phylum") %>%
  prune_samples(sample_sums(.) > 0, .) %>%  # Remove samples with zero Phylum-assigned reads
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>%  # Melt to long format
  arrange(Phylum)
phylumabundance_avg <- phylumabundance %>%
  select(Phylum, species, sp_abb, diet_cat, habitat_cat, Abundance) %>%
  group_by(Phylum, species, sp_abb, diet_cat, habitat_cat) %>%
  summarize(avg_abundance=mean(Abundance), .groups="drop") %>%
  mutate(Phylum = ifelse(avg_abundance < 0.01, "Other", Phylum)) %>% # Combine phyla < 1% into "Other"
  group_by(Phylum, species, sp_abb, diet_cat, habitat_cat) %>% # Sum phyla that were combined into "Other"
  summarize(avg_abundance=sum(avg_abundance), .groups="drop")
p <- ggplot(phylumabundance_avg) +
  geom_col(aes(x=species, y=avg_abundance, fill=fct_reorder(Phylum, avg_abundance, .desc=TRUE)), show.legend=TRUE) +
  facet_grid(.~diet_cat, scales="free", space="free") +
  ylab("Proportion of diet sequences") +
  scale_fill_manual(name="Phylum", values=my.cols) +
  xlab(NULL) + theme_bw() +
  theme(axis.text.y.left=element_text(size=15),
        axis.text.x=element_text(size=10, angle=45, hjust=1),
        axis.title.y=element_text(size=15))
p
ggsave("Figures/phyla_prop_by_species.pdf", plot=p)

### Closer look at taxonomic proportions (Class) per species
classabundance <- phyloseq_rarefied %>%
  tax_glom(taxrank = "Class") %>%
  prune_samples(sample_sums(.) > 0, .) %>%  # Remove samples with zero Class-assigned reads
  transform_sample_counts(function(x) x/sum(x)) %>%
  psmelt() %>%  # Melt to long format
  arrange(Class)
classabundance_avg <- classabundance %>%
  select(Class, species, sp_abb, diet_cat, habitat_cat, Abundance) %>%
  group_by(Class, species, sp_abb, diet_cat, habitat_cat) %>%
  summarize(avg_abundance=mean(Abundance), .groups="drop") %>%
  mutate(Class = ifelse(avg_abundance < 0.01, "Other", Class)) %>% # Combine phyla < 1% into "Other"
  group_by(Class, species, sp_abb, diet_cat, habitat_cat) %>% # Sum phyla that were combined into "Other"
  summarize(avg_abundance=sum(avg_abundance), .groups="drop")
p <- ggplot(classabundance_avg) +
  geom_col(aes(x=species, y=avg_abundance, fill=fct_reorder(Class, avg_abundance, .desc=TRUE)), show.legend=TRUE) +
  facet_grid(.~diet_cat, scales="free", space="free") +
  ylab("Proportion of diet sequences") +
  scale_fill_manual(name="Class", values=my.cols) +
  xlab(NULL) + theme_bw() +
  theme(axis.text.y.left=element_text(size=15),
        axis.text.x=element_text(size=10, angle=45, hjust=1),
        axis.title.y=element_text(size=15))
p
ggsave("Figures/class_prop_by_species.pdf", plot=p)

### Shannon diversity per species
my.cols <- c(brewer.pal(8,"Dark2"), rainbow_hcl(10))
plot_richness(phyloseq_rarefied, x="species", measures="Shannon", color="species") +
  geom_boxplot(alpha=0.6) + 
  scale_color_manual(values=my.cols) +
  theme_bw() + ylab("Shannon Diversity") +
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=10))

### Diversity stats per sample
alpha_ASV_meta <- data.frame(estimate_richness(physeq = phyloseq_rarefied, measures = c("Shannon", "Observed")),
                             Sample_ID=sample_data(phyloseq_rarefied)$fish_ID,
                             fish_species=sample_data(phyloseq_rarefied)$species,
                             diet_cat=sample_data(phyloseq_rarefied)$diet_cat,
                             habitat_cat=sample_data(phyloseq_rarefied)$habitat_cat)

### ANOVA to test for differences in Shannon diversity
richness.aov <- aov(Shannon ~ fish_species, data=alpha_ASV_meta)
summary(richness.aov)

### PCOA of species & diet category
ord.PCoA.jaccard <- ordinate(phyloseq_ReAb, method="PCoA", distance="jaccard")
pcoa.cols <- c(brewer.pal(8,"Dark2"), rainbow_hcl(10))
plot_ordination(phyloseq_ReAb, ord.PCoA.jaccard, type="samples", color="species") +
  stat_ellipse(type="norm") +
  geom_point(size=3, stroke=1) + 
  scale_color_manual(name="Host Species", values=pcoa.cols) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20), legend.title=element_text(size=20), legend.text=element_text(size=14))
plot_ordination(phyloseq_ReAb, ord.PCoA.jaccard, type="samples", color="diet_cat", shape="diet_cat") +
  stat_ellipse(type="norm") +
  geom_point(size=3, stroke=1) + 
  scale_color_manual(name="Host Diet", values=pcoa.cols) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20), legend.title=element_text(size=20), legend.text=element_text(size=20))

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
write.csv(alpha_div_metadata, "18S_alpha_div_metadata.csv")