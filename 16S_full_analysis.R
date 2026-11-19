### This code analyzes the 16S full run data (1 full Miseq + 2 spike-ins)
### This includes both environmental and gut-associated microbes 
### Inputs: ASV table, taxonomy, and phylogeny from qiime2
### Outputs: phyloseq objects (1. pre-processed data; 2. rarefied dataset), figures & analyses
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
library(colorspace)
library(forcats)
library(vegan)
library(ggpubr)
library(rstatix)
library(pander)
library(dplyr)
library(lefser) ## BiocManager::install("lefser")
library(dada2) ## BiocManager::install("dada2")
library(MiscMetabar)
library(microbiome) ## BiocManager::install("microbiome")
library(ggVennDiagram)
library(devtools)
library(pairwiseAdonis) ## install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(ggpmisc)
library(ANCOMBC)
library(ggtext)
library(patchwork)
library(lme4)
library(lmerTest)
library(MuMIn)
library(emmeans)
library(UpSetR)
library(tidyr)

setwd("16S_full/")

###############################
### Part I: data processing ###
###############################

### Produce phyloseq object
ASVs <- read_qza("16S_merged_table.qza")
taxonomy <- read_qza("taxonomy.qza")
metadata <- read_q2metadata("meta_data_R_updated.txt")
rownames(metadata) <- metadata$SampleID
metadata$SampleID <- as.factor(metadata$SampleID)
metadata$fish_ID <- as.factor(metadata$fish_ID)
metadata$species <- as.factor(metadata$species)
metadata$family <- as.factor(metadata$family)
metadata$sp_abb <- as.factor(metadata$sp_abb)
metadata$diet_cat <- as.factor(metadata$diet_cat)
metadata$sample_type <- as.factor(metadata$sample_type)
metadata$habitat_cat <- as.factor(metadata$habitat_cat)
sampledata <- sample_data(metadata)
rownames(sampledata) <- rownames(metadata)
R_tree <- read_qza("rooted-tree.qza")
phyloseq_dat <- phyloseq(otu_table(ASVs$data, taxa_are_rows = TRUE), tax_table(as.matrix(parse_taxonomy(taxonomy$data))), sampledata,phy_tree(R_tree$data))
phyloseq_dat # 25,090 taxa and 633 samples
sum(sample_sums(phyloseq_dat)) # 11,524,674 reads
saveRDS(phyloseq_dat,file = "phyloseq_dat.rds")
length(unique(phyloseq_dat@sam_data$fish_ID)) # 283 unique samples

### Remove ASVs that are not bacteria 
phyloseq_bact <- phyloseq_dat %>%
  subset_taxa(Kingdom == "d__Bacteria" &                       # only bacteria
      (is.na(Family) | Family != "Mitochondria") &   # filter out mitochondria
      (is.na(Genus) | Genus != "Chloroplast"))       # filter out chloroplasts

### Inspect negative controls 
phyloseq_blank <- subset_samples(phyloseq_bact, species=="Blank")
phyloseq_blank <- prune_taxa(taxa_sums(phyloseq_blank) > 0, phyloseq_blank)
nrow(phyloseq_blank@otu_table) # 31 ASVs present in blanks
phyloseq_blank_tab <- phyloseq_blank %>% psmelt()
my.cols <- c(brewer.pal(8,"Dark2"), "#C46254","bisque","darkorange", rainbow_hcl(10), rainbow(10), "azure3")
p <- ggplot(phyloseq_blank_tab, aes(x=fish_ID, y=Abundance, fill=Genus)) + 
  geom_col(width=0.9) +
  scale_fill_manual(name="Genus", values=my.cols) +
  scale_color_manual(name="Genus", values=my.cols) +
  ylab("Read Count") + xlab("Sample ID") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text=element_text(size=8), axis.title=element_text(size=14),
        legend.key.size = unit(0.2, 'cm'),
        legend.text = element_text(size=8))
p
ggsave("Figures/negative_controls.pdf", plot=p)

### Identify and remove contaminants by decontam
sample_data(phyloseq_bact)$is.neg <- sample_data(phyloseq_bact)$species == "Blank"
contam <- isContaminant(phyloseq_bact, method="prevalence", neg="is.neg", threshold=0.1)
contaminants <- rownames(contam)[which(contam$contaminant)] # 10 contaminants found
phyloseq_decontam <- prune_taxa(!rownames(phyloseq_bact@otu_table) %in% contaminants, phyloseq_bact)

### Remove negative controls
phyloseq_clean <- subset_samples(phyloseq_decontam, species!="Blank")

### Remove species with n=1 (Ameiurus natalis) 
phyloseq_clean <- subset_samples(phyloseq_clean, species!="Ameiurus_natalis")

### Inspect and remove samples with low reads (if one of the two replicates has low reads)
lowreads <- subset_samples(phyloseq_clean, sample_sums(phyloseq_clean) < 300) # 26 samples w/ <300 reads
SampleFilter <- lowreads@sam_data$fish_ID
phyloseq_nonzero <- subset_samples(phyloseq_clean,!fish_ID %in% SampleFilter)

### Transform to relative abundance and calculate correlation coefficient
phyloseq_ReAb <- transform_sample_counts(phyloseq_nonzero, function(OTU) OTU*100/sum(OTU))
samples_to_keep <- !grepl("reseq", sample_names(phyloseq_nonzero)) # remove any resequenced samples
phyloseq_ReAb_filtered <- prune_samples(samples_to_keep, phyloseq_nonzero)
Pearson_cor_table <- as.data.frame(cor(phyloseq_ReAb_filtered@otu_table), method="pearson") # calculate correlation coefficient
df <- data.frame(Samplename=NA, Pearson_values=NA, Reads=NA, Species=NA)
number_index <- seq(1,524, by=2)
for(i in number_index){
  df[i,1] <- rownames(Pearson_cor_table)[i]
  df[i,2] <- Pearson_cor_table[i, i+1]
  df[i,3] <- colSums(phyloseq_clean@otu_table)[i]+colSums(phyloseq_clean@otu_table)[i+1]
  df[i,4] <- as.vector(phyloseq_ReAb@sam_data$species)[i]
}
pearson_table <- na.omit(df)
mean(pearson_table$Pearson_values) # mean pearson correlation 0.984
quantile(pearson_table$Pearson_values)  # 0%: 0.774, 25%: 0.979, 50%: 0.998, 75%: 0.999, 100%: 0.999
p <- ggplot(pearson_table, aes(x=Species, y=Pearson_values, fill=Species, group=Species)) + 
  geom_boxplot(outlier.shape=NA, outlier.colour=NA) + 
  geom_point(shape=16, size=2, color="gray") + 
  labs(x="Species", y="Pearson correlation") + 
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20), legend.position="none", axis.text.x=element_text(angle=45, vjust=1, hjust=1))
p
ggsave("Figures/pearson_correlation_technical_reps.pdf", plot=p)

### Merge technical replicates and resequenced samples 
sample_data(phyloseq_clean)$fish_ID <- as.factor(sample_data(phyloseq_clean)$fish_ID) # merge replicates
phyloseq_clean.nonrep <- merge_samples(phyloseq_clean,"fish_ID")

### Metadata coerced to NA -- manually add back sample data
sample_table <- unique(sample_data(phyloseq_clean)[,c("fish_ID","sample_type", "family","species","sp_abb","diet_cat","habitat_cat")])
rownames(sample_table) <- sample_table$fish_ID
sample_data(phyloseq_clean.nonrep) <- sample_table
as.data.frame(with(sample_data(phyloseq_clean.nonrep), table(species)))

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
  geom_hline(yintercept=10000, linetype="dashed") +
  labs(x="Species", y="Read Depth") + 
  theme_bw() +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        legend.position = "none")
p
ggsave("Figures/read_depth_all_samples.pdf", plot=p)

### Determine rarefaction depth with ASV accumulation curves by species
accum_curves <- accu_plot(phyloseq_clean.nonrep, "species", add_nb_seq=TRUE, step=1000) + 
  # facet_wrap(~fact) + 
  theme_bw()
accum_curves
ggsave("Figures/16S_accum_curves.pdf", plot=accum_curves)

### Most samples plateau by 10,000 reads -- target rarefaction depth ###

### Identify samples with reads < 10,000
sample_sums(phyloseq_clean.nonrep@otu_table)
length(sort(sample_sums(prune_samples(sample_sums(phyloseq_clean.nonrep)<10000, phyloseq_clean.nonrep)))) # 14 samples below 10,000 reads
SampleFilter1 <- rownames(phyloseq_clean.nonrep@otu_table)[sample_sums(phyloseq_clean.nonrep@otu_table)<10000]
phyloseq.10000 <- subset_samples(phyloseq_clean.nonrep,!sample_names(phyloseq_clean.nonrep) %in% SampleFilter1)

### Produce final phyloseq after rarefaction
set.seed(10000)
phyloseq_rarefied <- rarefy_even_depth(phyloseq.10000, sample.size=min(sample_sums(phyloseq.10000)))
phyloseq_rarefied <- prune_taxa(taxa_sums(phyloseq_rarefied)>0,phyloseq_rarefied)
phyloseq_rarefied # 19,221 taxa and 258 samples
saveRDS(phyloseq_rarefied,file = "phyloseq_rarefied.rds")
sum(sample_sums(phyloseq_rarefied)) # 3,352,710 reads

### Sample size per species: 222 total gut microbiome samples retained, 36 environmental
as.data.frame(with(sample_data(phyloseq_rarefied), table(species)))

### Transform final dataset to relative abundance
phyloseq_ReAb <- transform_sample_counts(phyloseq_rarefied, function(x) x/sum(x))

###############################################################
### Part II: taxonomic overview by sample type/host species ###
###############################################################

# Phylum-level breakdown: plot top phyla
my.cols <- c(brewer.pal(8,"Dark2"), rainbow_hcl(10), rainbow(15), "azure3")
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
  tax_glom(taxrank = "Phylum") %>% # Aggregate at phylum level
  prune_samples(sample_sums(.) > 0, .) %>%  # Remove samples with zero Phylum-assigned reads
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%  # Melt to long format
  arrange(Phylum) 
phylumabundance_avg <- phylumabundance %>%
  select(Phylum, species, sample_type, diet_cat, habitat_cat, Abundance) %>%
  group_by(Phylum, species, sample_type, diet_cat, habitat_cat) %>%
  summarize(avg_abundance = mean(Abundance),  .groups = "drop") %>%
  mutate(Phylum = ifelse(avg_abundance < 0.01, "Other", Phylum)) %>% # Combine phyla <1% into "Other"
  group_by(Phylum, species, sample_type, diet_cat, habitat_cat) %>%
  summarize(avg_abundance=sum(avg_abundance), .groups="drop")
p <- ggplot(phylumabundance_avg) +
  geom_col(aes(x=species, y=avg_abundance, fill=fct_reorder(Phylum, avg_abundance, .desc=TRUE))) +
  facet_grid(.~diet_cat, scales="free", space="free") +
  ylab("Proportion of microbiome sequences") +
  scale_fill_manual(name="Phylum", values = my.cols) +
  xlab(NULL) + theme_bw() +
  theme(axis.text.y.left = element_text(size=16),
        axis.text.x = element_text(size=14, angle=45, hjust=1),
        axis.title.y = element_text(size=16))
p
ggsave("Figures/phylum_per_species.pdf", plot=p)

### Closer look at taxonomic proportions (Genus) per species
genusabundance <- phyloseq_rarefied %>%
  tax_glom(taxrank = "Genus") %>% # agglomerate at genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()  # Melt to long format
genusabundance_avg <- genusabundance %>%
  select(Genus, species, sp_abb, diet_cat, habitat_cat, Abundance) %>%
  group_by(Genus, species, sp_abb, diet_cat, habitat_cat) %>%
  summarize(avg_abundance = mean(Abundance)) %>%
  filter(avg_abundance > 0.05)
p <- ggplot(genusabundance_avg) +
  geom_col(mapping = aes(x=species, y=avg_abundance, fill=fct_reorder(Genus, avg_abundance, .desc=TRUE)), position="fill", show.legend=TRUE) +
  facet_grid(.~diet_cat, scales="free", space="free") +
  ylab("Proportion of microbiome sequences") +
  scale_fill_manual(name="Genus",values = my.cols) +
  xlab(NULL) + theme_bw()+
  theme(axis.text.y.left = element_text(size=15),
        axis.text.x = element_text(size=10, angle=45, hjust=1),
        axis.title.y = element_text(size=15))
p
ggsave("Figures/genus_per_species.pdf", plot=p)

##################################################################
### Part III: alpha/beta diversity by sample type/host species ###
##################################################################

### Shannon diversity per sample type
p <- plot_richness(phyloseq_rarefied, x="sample_type", measures="Shannon", color="sample_type") +
  geom_boxplot(alpha=0.6) + 
  scale_color_manual(values=c("#1B9E77", "#E6AB02", "#5DADE2")) +
  theme_bw() + ylab("Shannon diversity of microbial sequences") +
  theme(legend.position="none", axis.title=element_text(size=18), axis.text=element_text(size=16), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
ggsave("Figures/Shannon_diversity_sample_type.pdf", plot=p, width=6, height=5, units="in")

# ASV overlap between environmental and microbiome samples
p <- ggvenn_pq(phyloseq_rarefied, fact = "sample_type") + 
  ggplot2::scale_fill_distiller(palette = "BuPu", direction = 1)
p
ggsave("Figures/venn_diag_env_microbiome.pdf", plot=p)

### PCOA of species & environmental samples
ord.PCoA.bray <- ordinate(phyloseq_ReAb, method="PCoA", distance="bray")
pcoa.cols <- c(brewer.pal(8,"Dark2"), rainbow_hcl(10))
p <- plot_ordination(phyloseq_ReAb, ord.PCoA.bray, type="samples", color="sample_type") +
  geom_point(size=3, stroke=1) + 
  stat_ellipse(type="norm") +
  scale_color_manual(name="", values=c("#1B9E77", "#E6AB02", "#5DADE2")) +
  theme_bw() +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=20), legend.text=element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
ggsave("Figures/PCoA_gut_and_environment_microbiome.pdf", plot=p, width=8, height=5, units="in")

### Shannon diversity per species/sample type
all.samples.cols <- c("#D95F02","#1B9E77","#66A61E","#7570B3","#E7298A","#E6AB02","#A6761D","#666666",
                      "#E495A5","#D89F7F","#BDAB66","#96B56C","#65BC8C","#91ACE1","#C29DDE","#5DADE2")
p <- plot_richness(phyloseq_rarefied, x="species", measures="Shannon", color="species") +
  geom_boxplot(alpha=0.6) + 
  facet_grid(.~diet_cat, scales="free", space="free") +
  scale_color_manual(values=all.samples.cols) +
  theme_bw() + ylab("Shannon diversity of microbial sequences") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16), 
        legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=10),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
ggsave("Figures/Shannon_diversity_per_species.pdf", plot=p, width=8, height=6, units="in")

### Kruskall-Walllis test to assess differences in Shannon Diversity among sample types 
shannon_df <- estimate_richness(phyloseq_rarefied, measures = "Shannon")
shannon_df$sample_type <- sample_data(phyloseq_rarefied)$sample_type # Add metadata back in
kruskal.test(Shannon ~ sample_type, data = shannon_df)
library(FSA)
pairwise.wilcox.test(shannon_df$Shannon, shannon_df$sample_type, p.adjust.method = "bonferroni")

### Beta diversity differences in gut and environmental samples 
### PERMANOVA
meta <- microbiome::meta(phyloseq_ReAb)
dist_bray <- phyloseq::distance(phyloseq_ReAb, method = "bray")
adonis2(dist_bray ~ sample_type, data=meta, permutations = 999, by="terms")

### PERMDISP
dist_adund <- phyloseq::distance(phyloseq_ReAb, method="bray")
beta <- betadisper(dist_adund, meta$sample_type)
permutest(beta)

######################################################################
### Part IV: alpha/beta diversity by host trophic category/habitat ###
######################################################################

### Subset to just species
phyloseq_rarefied_fish <- subset_samples(phyloseq_rarefied, sample_type=="Microbiome")
phyloseq_ReAb_fish <- subset_samples(phyloseq_ReAb, sample_type=="Microbiome")

### Diversity stats per sample
alpha_ASV_meta <- data.frame(Shannon_ASV=diversity(phyloseq_rarefied_fish@otu_table,index="shannon"),
                             Observed_ASV=specnumber(phyloseq_rarefied_fish@otu_table),
                             Sample_ID=sample_data(phyloseq_rarefied_fish)$fish_ID,
                             fish_family=sample_data(phyloseq_rarefied_fish)$family,
                             fish_species=sample_data(phyloseq_rarefied_fish)$species,
                             diet_cat=sample_data(phyloseq_rarefied_fish)$diet_cat,
                             habitat_cat=sample_data(phyloseq_rarefied_fish)$habitat_cat)
mean(alpha_ASV_meta$Observed_ASV)
sd(alpha_ASV_meta$Observed_ASV)
max(alpha_ASV_meta$Observed_ASV)

### Top phyla in gut microbiome 
gut_phylum <- tax_glom(phyloseq_ReAb_fish, taxrank = "Phylum") # 47 unique phyla
gut_df <- psmelt(gut_phylum)
phylum_props <- gut_df %>% # Calculate mean proportion per phylum
  group_by(Phylum) %>%
  summarize(mean_abundance = mean(Abundance)) %>%
  arrange(desc(mean_abundance))
phylum_props

### Shannon diversity by trophic category
phyloseq_rarefied_fish <- subset_samples(phyloseq_rarefied, sample_type=="Microbiome")
p <- plot_richness(phyloseq_rarefied_fish, x="diet_cat", measures="Shannon", color="diet_cat") +
  geom_point(size=3) +
  geom_boxplot(alpha=0.6) + 
  scale_color_manual(values=c("#3058A8","#0EAA4B","#D3B42A","#DB559F")) +
  theme_bw() + ylab("Shannon diversity of gut microbiome") + xlab("Trophic category") +
  theme(legend.position="none", axis.text=element_text(size=16), axis.title = element_text(size=18), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
ggsave("Figures/Shannon_diversity_diet_cat.pdf", plot=p, width=5, height=5, units="in")

### Shannon diversity by habitat category
p <- plot_richness(phyloseq_rarefied_fish, x="habitat_cat", measures="Shannon", color="habitat_cat") +
  geom_point(size=3) +
  geom_boxplot(alpha=0.6) + 
  scale_color_manual(values=c("#712B90","#539997","#E75725")) +
  theme_bw() + ylab("Shannon diversity of gut microbiome") + xlab("Habitat category") +
  theme(legend.position="none", axis.text=element_text(size=16), axis.title = element_text(size=18), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
ggsave("Figures/Shannon_diversity_habitat_cat.pdf", plot=p)

### Mixed models to test for differences in Shannon diversity among diet and habitat
lmer_model_diet <- lmer(shannon ~ diet_cat + (1 | fish_species), data = alpha_ASV_meta)
summary(lmer_model_diet)
r.squaredGLMM(lmer_model_diet)
anova(lmer_model_diet, type = 3)
emmeans(lmer_model_diet, pairwise ~ diet_cat)

lmer_model_habitat <- lmer(shannon ~ habitat_cat + (1 | fish_species), data = alpha_ASV_meta)
summary(lmer_model_habitat)
r.squaredGLMM(lmer_model_habitat)
anova(lmer_model_habitat, type = 3)
emmeans(lmer_model_habitat, pairwise ~ habitat_cat)

### PCOA of species colored by diet, habitat, and species
ord.PCoA.bray <- ordinate(phyloseq_ReAb_fish, method="PCoA", distance="bray")
p <- plot_ordination(phyloseq_ReAb_fish, ord.PCoA.bray, type="samples", color="diet_cat") +
  stat_ellipse(type="norm") +
  geom_point(size=5, stroke=1) + 
  scale_color_manual(name="Trophic category", values=c("#3058A8","#0EAA4B","#D3B42A","#DB559F")) +
  theme_bw() +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20), legend.title=element_text(size=20), legend.text=element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
ggsave("Figures/PCoA_diet.pdf", plot=p, width=8, height=5, units="in")

p <- plot_ordination(phyloseq_ReAb_fish, ord.PCoA.bray, type="samples", color="habitat_cat") +
  stat_ellipse(type="norm") +
  geom_point(size=5, stroke=1) + 
  scale_color_manual(name="Habitat category", values=c("#712B90","#539997","#E75725")) +
  theme_bw() +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20), legend.title=element_text(size=20), legend.text=element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
ggsave("Figures/PCoA_habitat.pdf", plot=p, width=8, height=5, units="in")

fish.samples.cols <- c("#D95F02","#66A61E","#7570B3","#E7298A","#E6AB02","#A6761D","#666666",
                       "#E495A5","#D89F7F","#BDAB66","#96B56C","#65BC8C","#91ACE1","#C29DDE")

p <- plot_ordination(phyloseq_ReAb_fish, ord.PCoA.bray, color = "species") +
  geom_point(size=5, stroke=1) + 
  stat_ellipse(type="norm") +
  scale_color_manual(name="Host Species", values=fish.samples.cols) +
  theme_bw() +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20), legend.title=element_text(size=20), legend.text=element_text(size=16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
ggsave("Figures/PCoA_species.pdf", plot=p, width=8, height=5, units="in")

### PERMANOVA
meta <- microbiome::meta(phyloseq_ReAb_fish)
dist_bray <- phyloseq::distance(phyloseq_ReAb_fish, method = "bray")
adonis2(dist_bray ~ habitat_cat + diet_cat + species, data=meta, permutations = 999, by="terms")

### PERMDISP
dist_adund <- phyloseq::distance(phyloseq_ReAb_fish, method="bray")
beta <- betadisper(dist_adund, meta$habitat_cat)
permutest(beta)
beta <- betadisper(dist_adund, meta$diet_cat)
permutest(beta)
beta <- betadisper(dist_adund, meta$species)
permutest(beta)

#######################################################################
### Part V: top taxa contributing to differences among sample types ###
#######################################################################

### Generate PCoA with loadings for top contributing taxa
ordination <- ordinate(phyloseq_ReAb_fish, method = "PCoA", distance = "bray")
# Extract OTU table and PCoA sample scores
otu_matrix <- otu_table(phyloseq_ReAb_fish)
ordination_scores <- as.data.frame(ordination$vectors[, 1:2])  # First two axes (PC1 and PC2)
ordination_scores$SampleID <- rownames(ordination_scores)  # Add sample IDs
# Correlate taxa with ordination axes (axis 1 and axis 2)
taxa_correlation <- cor(otu_matrix, ordination_scores[, 1:2], use = "pairwise.complete.obs")
# Add taxonomic information to the correlations
taxa_info <- as.data.frame(tax_table(phyloseq_ReAb_fish))
taxa_correlation <- cbind(taxa_info, taxa_correlation)
# Inspect the top contributors for axis 1 and axis 2
top_taxa_PC1 <- taxa_correlation[order(abs(taxa_correlation[, "Axis.1"]), decreasing = TRUE), ][1:10, ]
top_taxa_PC2 <- taxa_correlation[order(abs(taxa_correlation[, "Axis.2"]), decreasing = TRUE), ][1:10, ]
# Subset top contributors for visualization
top_taxa <- unique(c(rownames(top_taxa_PC1), rownames(top_taxa_PC2)))
top_taxa_matrix <- taxa_correlation[top_taxa, ]
taxa_arrows <- as.data.frame(top_taxa_matrix[, c("Phylum", "Class", "Family", "Axis.1", "Axis.2")])
# Plot the ordination with taxa contributions
ordination_plot <- plot_ordination(phyloseq_ReAb_fish, ordination, color = "diet_cat") +
  geom_point(size=5, stroke=1) + 
  stat_ellipse(type="norm") +
  scale_color_manual(name="Diet category", values=c("#3058A8","#0EAA4B","#D3B42A","#DB559F")) +
  geom_segment(data=taxa_arrows, aes(x=0, y=0, xend=Axis.1, yend=Axis.2),
               arrow=arrow(length=unit(0.2, "cm")), color="lightgray", linewidth=0.5) +
  geom_text(data=taxa_arrows, aes(x=Axis.1, y=Axis.2, label=Class),
            color="black", size=5, vjust=1.5, hjust=0.4) +
  theme_bw() +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16), legend.title=element_text(size=20), legend.text=element_text(size=16),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ordination_plot
ggsave("Figures/PCOA_top_taxa_microbiome.pdf", plot=ordination_plot, width=10, height=7, units="in")

### ANCOMBC analysis to identify top taxa (at Genus level) differentiating trophic categories ###
phyloseq_genus <- tax_glom(phyloseq_rarefied_fish, taxrank="Genus", NArm=TRUE)
phyloseq_genus <- prune_taxa(taxa_sums(phyloseq_genus) > 0, phyloseq_genus) # Remove taxa with no counts
# Calculate prevalence = proportion of samples containing each genus
otu <- as(otu_table(phyloseq_genus), "matrix")
if (!taxa_are_rows(phyloseq_genus)) otu <- t(otu)
prevalence <- rowMeans(otu > 0)
prevalence_df <- data.frame(Taxon = names(prevalence), Prevalence=prevalence)
ggplot(prevalence_df, aes(x=Prevalence)) + 
  geom_histogram(binwidth=0.05, boundary=0) +
  geom_vline(xintercept=0.05, linetype="dashed") +
  scale_x_continuous(breaks=seq(0, 1, 0.1), limits=c(0, 1)) +
  labs(x="Prevalence across samples", y="Number of genera") +
  theme_bw() + theme(axis.title=element_text(size = 16), axis.text=element_text(size=14))
# Use 5% prevalence as cutoff in analysis for differentially abundant taxa
# Run ANCOM-BC2
ancom_diet <- ancombc2(data=phyloseq_genus, tax_level=NULL, fix_formula="diet_cat",
                       group="diet_cat", p_adj_method="BH", prv_cut=0.05, lib_cut=0, struc_zero=TRUE,
                       neg_lb=TRUE, global=TRUE, pairwise=TRUE)
top_taxa <- ancom_diet$res_global %>% filter(q_val < 0.05)
nrow(top_taxa) # 166 taxa w/ q-values <0.05
# Subset overall taxa list to just those with the largest effect sizes in the pairwise comparisons
pairwise_lfc <- c("lfc_diet_catHerbivore", "lfc_diet_catInvertivore", "lfc_diet_catOmnivore", 
                  "lfc_diet_catInvertivore_diet_catHerbivore", "lfc_diet_catOmnivore_diet_catHerbivore", "lfc_diet_catOmnivore_diet_catInvertivore")
effect_sizes <- ancom_diet$res_pair %>% rowwise() %>%
  mutate(max_abs_lfc=max(abs(c_across(all_of(pairwise_lfc))), na.rm=TRUE)) %>%
  ungroup()
top_taxa <- ancom_diet$res_global %>%
  filter(q_val < 0.05) %>%
  select(taxon, q_val) %>%
  left_join(effect_sizes %>% select(taxon, max_abs_lfc), by="taxon") %>%
  arrange(desc(max_abs_lfc)) %>%
  head(20) # Select top 20 taxa
phyloseq_genus_RelAb <- transform_sample_counts(phyloseq_genus, function(x) x / sum(x)) # Transform to rel. abund.
genus_RelAb <- psmelt(phyloseq_genus_RelAb) # Melt to long format
top_taxa_ids <- top_taxa$taxon
heatmap_dat <- genus_RelAb %>% # Filter to top taxa
  filter(OTU %in% top_taxa_ids) %>%
  group_by(OTU, Phylum, Genus, diet_cat) %>%
  summarize(mean_abundance = mean(Abundance), .groups = "drop") %>%
  mutate(Taxon_label = paste0(Phylum," (<i>", Genus, "</i>)")) %>%
  left_join(top_taxa %>% select(taxon, max_abs_lfc), by=c("OTU"="taxon"))
ancombc_plot_microbiome <- ggplot(heatmap_dat, aes(x = diet_cat, y=reorder(Taxon_label, max_abs_lfc), fill = mean_abundance)) +
  geom_tile(color = "white") + 
  scale_fill_viridis_c(option="viridis", name="Mean relative\nabundance") +
  labs(x = "Trophic category", y = NULL) + theme_bw() +
  theme(axis.text.x = element_text(size = 12),
    axis.text.y = element_markdown(size = 8), axis.title.x = element_text(size = 16), 
    legend.title = element_text(size = 12), legend.text = element_text(size = 11))
ancombc_plot_microbiome
ggsave("Figures/ancombc_plot_microbiome.pdf", plot=ancombc_plot_microbiome, width=7, height=4, units="in")

##################################################
### Part V: core gut microbiome across hosts ###
##################################################

### Plot abundance-prevalence relationship for fish gut microbiome
fish_counts <- psmelt(phyloseq_rarefied_fish)
fish_RelAb <- psmelt(phyloseq_ReAb_fish)
# Calculate prevalence across fish samples
prevalence_dat <- fish_counts %>%
  group_by(OTU) %>%
  summarize(prevalence=mean(Abundance > 0), .groups="drop")
# Calculate mean relative abundance across fish samples
abundance_dat <- fish_RelAb %>%
  group_by(OTU, Phylum, Class, Family, Genus) %>%
  summarize(mean_abundance=mean(Abundance), .groups="drop")
# Combine prevalence and abundance
core_plot_dat <- left_join(abundance_dat, prevalence_dat, by="OTU") %>%
  mutate(Core=ifelse(prevalence >= 0.80, "Core", "Non-core"))
p <- ggplot(core_plot_dat, aes(x=mean_abundance, y=prevalence, color=Core)) +
  geom_point(size=3, alpha=0.7) +
  geom_hline(yintercept=0.80, linetype="dashed") +
  scale_x_log10() +
  labs(x="log (Mean relative abundance)", y="Prevalence") +
  theme_bw() +
  theme(axis.title=element_text(size=16), axis.text=element_text(size=14),
        legend.title=element_text(size=14), legend.text=element_text(size=12))
p
ggsave("Figures/prevalence_vs_abundance_gut_microbiome.pdf", plot=p, width=8, height=6, units="in")

### UpSet plot to show core taxa sharing across trophic categories ###
# Convert taxa to presence/absence
otu <- as(otu_table(phyloseq_rarefied_fish), "matrix")
if (taxa_are_rows(phyloseq_rarefied_fish)) otu <- t(otu)
otu_pa <- otu > 0
metadata <- data.frame(sample_data(phyloseq_rarefied_fish))
# Identify overall core taxa (>=80% prevalence)
overall_prevalence <- colMeans(otu_pa)
core_taxa_ids <- names(overall_prevalence[overall_prevalence >= 0.80])
# Calculate prevalence of overall core taxa within each trophic category
diet_core <- lapply(unique(metadata$diet_cat), function(diet) {
  samples <- rownames(metadata)[metadata$diet_cat == diet]
  data.frame(OTU=core_taxa_ids,
             diet_cat=diet,
             prevalence=colMeans(otu_pa[samples, core_taxa_ids, drop=FALSE]))}) %>% 
  bind_rows() %>% mutate(core_within_diet=prevalence>0.80)
# Convert to wide 0/1 format for UpSetR
diet_core_wide <- diet_core %>%
  select(OTU, diet_cat, core_within_diet) %>%
  pivot_wider(names_from = diet_cat, values_from = core_within_diet, values_fill = FALSE)
# Define trophic category order
diet_names <- c("Generalist predator", "Herbivore","Invertivore", "Omnivore")
# Colors for trophic categories
diet_colors <- c("Generalist predator"="#3058A8", "Herbivore"="#0EAA4B", "Invertivore"="#D3B42A", "Omnivore"="#DB559F")
# Format metadata for UpSetR
diet_metadata <- data.frame(sets=diet_names, diet_cat=diet_names)
# Plot
pdf("Figures/core_microbiome_overlap_diet.pdf", width = 8, height = 5)
upset(diet_upset,
      sets = diet_names, keep.order = TRUE, order.by = "freq",
      set.metadata=list(data=diet_metadata, plots=list(list(type="matrix_rows", column="diet_cat", colors=diet_colors, alpha=0.5))))
dev.off()

############################################################
### Part VI: core gut microbiome compared to environment ###
############################################################

otu_mat <- as(otu_table(phyloseq_rarefied), "matrix")
if (!taxa_are_rows(phyloseq_rarefied)) otu_mat <- t(otu_mat)
# Extract sample metadata
metadata <- data.frame(sample_data(phyloseq_rarefied))
metadata <- metadata[colnames(otu_mat), , drop=FALSE]
# Convert to presence/absence
otu_pa <- otu_mat > 0
# Identify core gut taxa (>=70% prevalence)
gut_samples <- rownames(metadata)[metadata$sample_type == "Microbiome"]
gut_prevalence <- rowMeans(otu_pa[, gut_samples, drop=FALSE])
core_taxa <- names(gut_prevalence[gut_prevalence >= 0.80])

# Calculate prevalence of core taxa within each sample type
prevalence_dat <- lapply(unique(metadata$sample_type), function(stype) {
  samples <- rownames(metadata)[metadata$sample_type == stype]
  data.frame(Taxon=core_taxa,
             prevalence=rowMeans(otu_pa[core_taxa, samples, drop=FALSE]),
             sample_type=stype)}) %>%
  bind_rows()
# Calculate mean relative abundance of core taxa in gut microbiome only
gut_rel <- sweep(otu_mat[, gut_samples, drop=FALSE], 2,
                 colSums(otu_mat[, gut_samples, drop=FALSE]), "/")
abundance_dat <- data.frame(Taxon=core_taxa, 
                            mean_abundance=rowMeans(gut_rel[core_taxa, , drop=FALSE]))
# Add taxonomy
tax_df <- as.data.frame(tax_table(phyloseq_rarefied))
tax_df$Taxon <- rownames(tax_df)
prevalence_dat <- prevalence_dat %>%
  left_join(tax_df, by="Taxon") %>%
  mutate(Phylum=ifelse(is.na(Phylum), Kingdom, Phylum), 
         Genus=ifelse(is.na(Genus), "Unclassified", Genus), 
         taxon_label=paste0(Phylum, " (", Genus, ")"))
# Order taxa by gut prevalence
taxon_order <- prevalence_dat %>%
  filter(sample_type == "Microbiome") %>%
  arrange(prevalence) %>%
  pull(Taxon)
prevalence_dat$Taxon <- factor(prevalence_dat$Taxon, levels=taxon_order)
abundance_dat$Taxon <- factor(abundance_dat$Taxon, levels=taxon_order)
# Set sample type order and colors
prevalence_dat$sample_type <- factor(prevalence_dat$sample_type, levels=c("Microbiome", "Biofilm", "Water"))
sample_colors <- c("Microbiome"="#E6AB02", "Biofilm"="#1B9E77", "Water"="#5DADE2")

# Prevalence plot
p1 <- ggplot(prevalence_dat, aes(x=prevalence, y=Taxon, group=Taxon, color=sample_type)) +
  geom_line(color="darkgray", linewidth=0.4) +
  geom_point(size=5) +
  scale_y_discrete(labels=setNames(prevalence_dat$taxon_label, prevalence_dat$Taxon)) +
  scale_color_manual(values=sample_colors) +
  labs(y="Phylum (Genus)", x="Prevalence (proportion of samples)", color="Sample Type") +
  theme_bw() +
  theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=16),
        legend.title=element_text(size=20), legend.text=element_text(size=20),
        panel.grid.minor = element_blank())

# Mean relative abundance in gut microbiome
p2 <- ggplot(abundance_dat, aes(x=mean_abundance, y=Taxon)) +
  geom_col(fill="black", width=0.7) +
  labs(y=NULL, x="Mean relative abundance\n(gut microbiome)") +
  theme_bw() +
  theme(axis.text.x=element_text(size=12), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title=element_text(size=16), 
        panel.grid.minor = element_blank())

# Combine
p <- p1 + p2 + plot_layout(widths=c(1, 0.5), guides="collect") & theme(legend.position="bottom")
p
ggsave("Figures/core_gut_microbiome.pdf", plot=p, width=10, height=8, units="in")

# Calculate average prevalence of these taxa 
prevalence_dat %>%
  group_by(sample_type) %>%
  summarise(mean_prevalence = mean(prevalence)) %>%
  arrange(sample_type)

# Look at phyla of core microbiome 
prevalence_dat %>%
  distinct(Taxon, Phylum) %>%
  dplyr::count(Phylum, sort = TRUE)
