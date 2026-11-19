### This code combines outputs from separate loci and completes several analyses 
### Inputs: 16S, 18S, COI, and 12S rarified phyloseq objects, phylogenetic distances
### Outputs: statistics, models, figures
### Last updated 22 Oct 2026 by Kara Andres (kjandr2@ilstu.edu)

### Clear work environment and load packages
library(phyloseq)
library(vegan)
library(ggplot2)
library(plyr)
library(dplyr)
library(colorspace)
library(RColorBrewer)
library(forcats)
library(ecodist)
library(ape)
library(cluster)
library(lme4)
library(lmerTest)
library(MuMIn)
library(emmeans)
library(ANCOMBC)
library(ggtext)

# Load the phyloseq objects for all loci
phy16S <- readRDS("16S_full/phyloseq_rarefied.rds")
phy18S <- readRDS("18S_full/phyloseq_rarefied.rds")
phyCOI <- readRDS("COI_full/phyloseq_rarefied.rds")
phy12S <- readRDS("12S_full/phyloseq_rarefied.rds")

# Remove environmental samples from 16S
microbiome_phyloseq <- subset_samples(phy16S, sample_type=="Microbiome")
microbiome_phyloseq_rel <- transform_sample_counts(microbiome_phyloseq, function(OTU) OTU*100/sum(OTU))

# Remove the phylogenetic tree from diet phyloseq objects
phy18S_no_tree <- phyloseq(otu_table(phy18S), sample_data(phy18S), tax_table(phy18S))
phy12S_no_tree <- phyloseq(otu_table(phy12S), sample_data(phy12S), tax_table(phy12S))
phyCOI_no_tree <- phyloseq(otu_table(phyCOI), sample_data(phyCOI), tax_table(phyCOI))

# Merge the diet metabarcoding phyloseq objects
diet_phyloseq <- merge_phyloseq(phy18S_no_tree, phy12S_no_tree, phyCOI_no_tree)
diet_phyloseq
# 5796 taxa and 208 samples

#############################################################################
### Part I: Diet alpha diversity across species, diet, and habitat categories
#############################################################################

### Diet ASV richness across species
fish.samples.cols <- c("#D95F02","#66A61E","#7570B3","#E7298A","#E6AB02","#A6761D","#666666",
                       "#E495A5","#D89F7F","#BDAB66","#96B56C","#65BC8C","#91ACE1","#C29DDE")
p <- plot_richness(diet_phyloseq, x="species", measures="Observed", color="species") +
  geom_point(size=3) +
  geom_boxplot(alpha=0.6) + 
  scale_color_manual(values=fish.samples.cols) +
  facet_grid(.~diet_cat, scales="free", space="free") +
  theme_bw() + ylab("Diet ASV richness") +
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=10))
p
ggsave("Figures/Diet_ASV_richness_per_species.pdf", plot=p, width=8, height=6, units="in")

### Diet ASV richness by trophic category
p <- plot_richness(diet_phyloseq, x="diet_cat", measures="Observed", color="diet_cat") +
  geom_point(size=3) +
  geom_boxplot(alpha=0.6) + 
  scale_color_manual(values=c("#3058A8","#0EAA4B","#D3B42A","#DB559F")) +
  theme_bw() + ylab("ASV richness of diet") + xlab("Trophic Category") +
  theme(legend.position="none", axis.text=element_text(size=16), axis.title = element_text(size=18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
ggsave("Figures/Diet_ASV_richness_diet_cat.pdf", plot=p, width=5, height=5, units="in")

### Diet ASV richness by habitat category
p <- plot_richness(diet_phyloseq, x="habitat_cat", measures="Observed", color="habitat_cat") +
  geom_point(size=3) +
  geom_boxplot(alpha=0.6)+ 
  scale_color_manual(values=c("#712B90","#539997","#E75725")) +
  theme_bw() + ylab("Observed ASV diversity") + xlab("Habitat Category") +
  theme(legend.position="none", axis.text=element_text(size=16), axis.title = element_text(size=18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
ggsave("Figures/Diet_ASV_richness_habitat_cat.pdf", plot=p, width=5, height=5, units="in")

### Diversity stats per sample
alpha_ASV_meta <- data.frame(Observed_ASV=specnumber(diet_phyloseq@otu_table),
                             Sample_ID=sample_data(diet_phyloseq)$fish_ID,
                             fish_species=sample_data(diet_phyloseq)$species,
                             diet_cat=sample_data(diet_phyloseq)$diet_cat,
                             habitat_cat=sample_data(diet_phyloseq)$habitat_cat)

# Statistical test for differences in ASV richness across samples 
lmer_model_diet <- lmer(Observed_ASV ~ diet_cat + (1 | fish_species), data = alpha_ASV_meta)
summary(lmer_model_diet)
r.squaredGLMM(lmer_model_diet)
anova(lmer_model_diet, type = 3)
emmeans(lmer_model_diet, pairwise ~ diet_cat)

lmer_model_habitat <- lmer(Observed_ASV ~ habitat_cat + (1 | fish_species), data = alpha_ASV_meta)
summary(lmer_model_habitat)
r.squaredGLMM(lmer_model_habitat)
anova(lmer_model_habitat, type = 3)
emmeans(lmer_model_habitat, pairwise ~ habitat_cat)

#############################################################################
### Part II: Diet beta diversity across species, diet, and habitat categories
#############################################################################

### PCOA by diet category
dist_jaccard <- phyloseq::distance(diet_phyloseq, method="jaccard", binary=TRUE)
diet_ord <- ordinate(diet_phyloseq, method="PCoA", distance=dist_jaccard)
p <- plot_ordination(diet_phyloseq, diet_ord, type="samples", color="diet_cat") +
  stat_ellipse(type="norm") +
  geom_point(size=5, stroke=1) + 
  scale_color_manual(name="Diet category", values=c("#3058A8","#0EAA4B","#D3B42A","#DB559F")) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20), legend.title=element_text(size=20), legend.text=element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
ggsave("Figures/dietmet_PCoA_diet.pdf", plot=p, width=8, height=5, units="in")

### PCOA by habitat category
p <- plot_ordination(diet_phyloseq, diet_ord, type="samples", color="habitat_cat") +
  stat_ellipse(type="norm") +
  geom_point(size=5, stroke=1) + 
  scale_color_manual(name="Habitat category", values=c("#712B90","#539997","#E75725")) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=20), legend.title=element_text(size=20), legend.text=element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
ggsave("Figures/dietmet_PCoA_habitat.pdf", plot=p, width=8, height=5, units="in")

### PCOA by species
fish.samples.cols <- c("#D95F02","#66A61E","#7570B3","#E7298A","#E6AB02","#A6761D","#666666",
                       "#E495A5","#D89F7F","#BDAB66","#96B56C","#65BC8C","#91ACE1","#C29DDE")
p <- plot_ordination(diet_phyloseq, diet_ord, color = "species") +
  geom_point(size=5, stroke=1) + 
  stat_ellipse(type="norm") +
  scale_color_manual(name="Host Species", values=fish.samples.cols) +
  theme_bw() +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16), legend.title=element_text(size=20), legend.text=element_text(size=16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
ggsave("Figures/dietmet_PCoA_species.pdf", plot=p, width=8, height=5, units="in")

### Generate PCoA with loadings for top contributing taxa
phy_rel <- transform_sample_counts(diet_phyloseq, function(x) x / sum(x)) # Normalize to relative abundance
# Create PCoA ordination using Jaccard distance
ordination <- ordinate(diet_phyloseq, method="PCoA", distance="jaccard")
# Extract OTU table and PCoA sample scores
otu_matrix <- otu_table(diet_phyloseq)
ordination_scores <- as.data.frame(ordination$vectors[, 1:2])  # First two axes (PC1 and PC2)
ordination_scores$SampleID <- rownames(ordination_scores)  # Add sample IDs
# Correlate taxa with ordination axes (axis 1 and axis 2)
taxa_correlation <- cor(otu_matrix, ordination_scores[, 1:2], use = "pairwise.complete.obs")
# Add taxonomic information to the correlations
taxa_info <- as.data.frame(tax_table(diet_phyloseq))
taxa_correlation <- cbind(taxa_info, taxa_correlation)
# Inspect the top contributors for axis 1 and axis 2
top_taxa_PC1 <- taxa_correlation[order(abs(taxa_correlation[, "Axis.1"]), decreasing = TRUE), ][1:10, ]
top_taxa_PC2 <- taxa_correlation[order(abs(taxa_correlation[, "Axis.2"]), decreasing = TRUE), ][1:10, ]
# Subset top contributors for visualization
top_taxa <- unique(c(rownames(top_taxa_PC1), rownames(top_taxa_PC2)))
top_taxa_matrix <- taxa_correlation[top_taxa, ]
taxa_arrows <- as.data.frame(top_taxa_matrix[, c("Phylum", "Class", "Family", "Axis.1", "Axis.2")])
# Plot the ordination with taxa contributions
ordination_plot <- plot_ordination(diet_phyloseq, ordination, color = "diet_cat") +
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
ggsave("Figures/PCOA_top_taxa_diet.pdf", plot=ordination_plot, width=10, height=7, units="in")

### ANCOMBC analysis to identify top taxa (at Genus level) differentiating trophic categories ###
phyloseq_genus <- tax_glom(diet_phyloseq, taxrank="Genus", NArm=TRUE)
phyloseq_genus <- prune_taxa(taxa_sums(phyloseq_genus) > 0, phyloseq_genus) # Remove taxa with no counts
phyloseq_genus <- prune_samples(sample_sums(phyloseq_genus) > 0, phyloseq_genus) # Remove samples with no counts
# Calculate prevalence = proportion of samples containing each genus
otu <- as(otu_table(phyloseq_genus), "matrix")
prevalence <- colSums(otu > 0) / nrow(otu)
prevalence_df <- data.frame(Taxon=names(prevalence), Prevalence=prevalence)
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
nrow(top_taxa) # 39 taxa w/ q-values <0.05
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
ancombc_plot_diet <- ggplot(heatmap_dat, aes(x=diet_cat, y=reorder(Taxon_label, max_abs_lfc), fill=mean_abundance)) +
  geom_tile(color="white") + 
  scale_fill_viridis_c(option="viridis", name="Mean relative\nabundance") +
  labs(x = "Trophic category", y = NULL) + theme_bw() +
  theme(axis.text.x=element_text(size = 12),
        axis.text.y=element_markdown(size=8), axis.title.x=element_text(size=16), 
        legend.title=element_text(size=12), legend.text=element_text(size=11))
ancombc_plot_diet
ggsave("Figures/ancombc_plot_diet.pdf", plot=ancombc_plot_diet, width=6.5, height=5, units="in")

### PERMANOVA
meta <- microbiome::meta(diet_phyloseq)
dist_jaccard <- phyloseq::distance(diet_phyloseq, method = "jaccard", binary=TRUE)
adonis2(dist_jaccard ~ habitat_cat + diet_cat + species, data=meta, permutations = 999, by="terms")

### PERMDISP
beta <- betadisper(dist_jaccard, meta$habitat_cat)
permutest(beta)
beta <- betadisper(dist_jaccard, meta$diet_cat)
permutest(beta)
beta <- betadisper(dist_jaccard, meta$species)
permutest(beta)

#########################################################################
### Part III: Regression between microbiome diversity and diet diversity
#########################################################################

# Calculate richness for microbiome
microb_alpha <- estimate_richness(microbiome_phyloseq, measures="Shannon") %>%
  tibble::rownames_to_column("sample_id") %>%
  rename(microb_diversity = Shannon)
# Calculate richness for diet
diet_alpha <- estimate_richness(diet_phyloseq, measures="Observed") %>%
  tibble::rownames_to_column("sample_id") %>%
  rename(diet_richness = Observed)
# Merge by sample_id
alpha_df <- inner_join(microb_alpha, diet_alpha, by="sample_id")
# Add metadata 
meta_df <- data.frame(sample_data(microbiome_phyloseq)) %>%
  tibble::rownames_to_column("sample_id")
alpha_df <- left_join(alpha_df, meta_df, by = "sample_id")
# Linear regression
lm_model <- lm(microb_diversity ~ diet_richness + diet_cat + habitat_cat, data = alpha_df)
library(car)
vif(lm_model)
plot(fitted(lm_model), resid(lm_model))
qqline(resid(lm_model), col = "red")
abline(h = 0, col = "red")
model_summary <- summary(lm_model)
anova(lm_model)
r2 <- round(model_summary$r.squared, 2)
pval <- signif(model_summary$coefficients[2, 4], 3)

p <- ggplot(alpha_df, aes(x=diet_richness, y=microb_diversity, color=diet_cat)) +
  geom_point(size=3) +
  geom_smooth(method="lm", se=TRUE, color="darkgray") +
  labs(x = "Diet ASV richness", y = "Microbiome Shannon diversity") +
  scale_color_manual(name="Diet category", values=c("#3058A8","#0EAA4B","#D3B42A","#DB559F")) +
  annotate("text", x=Inf, y=Inf, hjust=1.1, vjust=2, label = paste0("R² = ", r2, ", p = ", pval), size=8) +
  theme_bw() +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20), legend.title=element_text(size=20), legend.text=element_text(size=20),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank())
p
ggsave("Figures/diet_div_microb_div_regression.pdf", plot=p, width=9, height=6, units="in")

###########################################################################
### Part IV: MRM: relative contributions of diet & phylogeny to microbiome
###########################################################################
# Subset each phyloseq object to only the samples in common between them
microbiome_samples <- sample_names(microbiome_phyloseq)
diet_samples <- sample_names(diet_phyloseq)
common_samples <- intersect(microbiome_samples, diet_samples) # Find common sample names
microbiome_common <- prune_samples(common_samples, microbiome_phyloseq) # Subset both phyloseq objects to only the common samples
diet_common <- prune_samples(common_samples, diet_phyloseq)
microbiome_common # 19221 taxa and 201 samples
diet_common # 5796 taxa and 201 samples

# Calculate Bray-Curtis distance matrix for microbiome and Jaccard for diet data
microbiome_otu <- as(otu_table(microbiome_common), "matrix")
diet_otu <- as(otu_table(diet_common), "matrix")
microbiome_dist <- vegdist(microbiome_otu, method="bray")
diet_dist <- vegdist(diet_otu, method="jaccard", binary=TRUE)

# Get metadata in the same sample order
meta_data <- data.frame(sample_data(microbiome_common), stringsAsFactors=FALSE)
meta_data <- meta_data[rownames(microbiome_otu), ]

# Create individual-level phylogenetic distance matrix
phylo_dist <- read.csv("cytb_distMatrix.csv", row.names=1)
individual_phylo_dist <- outer(meta_data$species, meta_data$species,
                               Vectorize(function(x, y) phylo_dist[x, y]))
rownames(individual_phylo_dist) <- rownames(meta_data)
colnames(individual_phylo_dist) <- rownames(meta_data)
individual_phylo_dist <- as.dist(individual_phylo_dist)
plot(hclust(individual_phylo_dist), cex=0.4)

# Create individual-level habitat category distance matrix 
individual_habitat <- data.frame(habitat_cat=meta_data$habitat_cat)
rownames(individual_habitat) <- rownames(meta_data)
individual_habitat_dist <- cluster::daisy(individual_habitat, metric="gower")

# Create individual-level species identity matrix
species_vec <- meta_data$species
species_dist <- outer(species_vec, species_vec, FUN=function(x, y) ifelse(x == y, 0, 1))
rownames(species_dist) <- rownames(meta_data)
colnames(species_dist) <- rownames(meta_data)
species_dist <- as.dist(species_dist)

# Create individual-level diet distance matrix
diet_vec <- meta_data$diet_cat
diet_cat_dist <- outer(diet_vec, diet_vec, FUN=function(x, y) ifelse(x == y, 0, 1))
rownames(diet_cat_dist) <- rownames(meta_data)
colnames(diet_cat_dist) <- rownames(meta_data)
diet_cat_dist <- as.dist(diet_cat_dist)

# Check collinearity among predictor distance matrices
vegan::mantel(diet_dist, individual_phylo_dist, permutations = 999) # r = 0.2151
vegan::mantel(diet_dist, individual_habitat_dist, permutations = 999) # r = 0.1187
vegan::mantel(individual_phylo_dist, species_dist, permutations = 999) # r = 0.9038
# drop species-level dissimilarity matrix due to high collinearity

# Microbiome dissimilarity explained by diet metabarcoding, phylogeny, and habitat differences
mrm_data <- data.frame(microbiome=as.vector(microbiome_dist),
                       diet=as.vector(diet_dist),
                       phylogeny=as.vector(individual_phylo_dist),
                       habitat=as.vector(individual_habitat_dist))
mrm_model <- ecodist::MRM(microbiome ~ diet + phylogeny + habitat, data=mrm_data, nperm=9999)
mrm_model
predictors <- c("Intercept", "Diet dissimilarity", "Phylogenetic distance", "Habitat dissimilarity")
# Extract coefficients and p-values
coef_vals <- mrm_model$coef[, 1]
p_vals <- mrm_model$coef[, 2]
# Create labeled data frame
mrm_results <- data.frame(Predictor=predictors, Estimate=round(coef_vals, 3), P_value=round(p_vals, 3))
predictors <- c("Habitat dissimilarity", "Phylogenetic distance", "Diet dissimilarity")
mrm_results$Predictor <- factor(mrm_results$Predictor, levels=predictors)
mrm_results$SigLabel <- ifelse(mrm_results$P_value < 0.05, "*", "")
# Bar plot of effect sizes
p <- ggplot(mrm_results[-1, ], aes(x=Predictor, y=Estimate)) +
  geom_col(fill="steelblue", width=0.7) + coord_flip() +
  geom_text(aes(label=SigLabel), vjust=0.7, hjust=-0.3, size=10) +
  labs(x=NULL, y="Effect size") +
  theme_bw() +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=20), legend.title=element_text(size=20), legend.text=element_text(size=16),
        panel.grid.minor=element_blank())
p
ggsave("Figures/MRM_diet_metabarcoding.pdf", plot=p, width=6, height=6, units="in")

# Microbiome dissimilarity explained by diet category, phylogeny, and habitat differences
mrm_data <- data.frame(microbiome=as.vector(microbiome_dist),
                       diet=as.vector(diet_cat_dist),
                       phylogeny=as.vector(individual_phylo_dist),
                       habitat=as.vector(individual_habitat_dist))
mrm_model <- ecodist::MRM(microbiome ~ diet + phylogeny + habitat, data=mrm_data, nperm=9999)
mrm_model
predictors <- c("Intercept", "Diet dissimilarity", "Phylogenetic distance", "Habitat dissimilarity")
# Extract coefficients and p-values
coef_vals <- mrm_model$coef[, 1]
p_vals <- mrm_model$coef[, 2]
# Create labeled data frame
mrm_results <- data.frame(Predictor=predictors, Estimate=round(coef_vals, 3), P_value=round(p_vals, 3))
predictors <- c("Habitat dissimilarity", "Phylogenetic distance", "Diet dissimilarity")
mrm_results$Predictor <- factor(mrm_results$Predictor, levels=predictors)
mrm_results$SigLabel <- ifelse(mrm_results$P_value < 0.05, "*", "")
# Bar plot of effect sizes
p <- ggplot(mrm_results[-1, ], aes(x=Predictor, y=Estimate)) +
  geom_col(fill="steelblue", width=0.7) + coord_flip() +
  geom_text(aes(label=SigLabel), vjust=0.7, hjust=-0.3, size=10) +
  labs(x=NULL, y="Effect size") +
  theme_bw() +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=20), legend.title=element_text(size=20), legend.text=element_text(size=16),
        panel.grid.minor=element_blank())
p
ggsave("Figures/MRM_diet_categorical.pdf", plot=p, width=6, height=6, units="in")
