### This code calculates the phylogenetic similarity between species
### Last updated 26 February 2025 by Kara Andres (kjandr2@ilstu.edu)

### Clear work environment and load packages
rm(list = ls())
library(traits) # devtools::install_github("ropensci/traits")
library(msa) # BiocManager::install("msa")
library(seqinr)
library(ape)
library(ggplot2)
library(ggtree)
library(vegan)
library(phyloseq)
library(reshape2)

# Download sequences from NCBI by accession number
# Full mitogenomes (only cytb for Cottus_carolinae and Etheostoma_erythrozonum)
# seqs <- ncbi_byid(ids=c("JN389795.2","AP012082.1","KY660677.1","OR552032.1","OR552056.1","AP006776.1",
#                         "KY660678.1","DQ536421.1","OR552040.1","AB378749.1","OR002151.1",
#                         "OR002153.1","AB242169.1","MF621735.1","HQ128171.1", verbose=TRUE))
seqs <- ncbi_byid(ids=c("NC_035659.1", # Ambloplites rupestris mitochondrion, complete genome
                        "NC_008102.1", # Campostoma anomalum mitochondrion, complete genome
                        "AF549113.1", # Cottus carolinae isolate CB36 cytochrome b (cytb) gene
                        "NC_035660.1", # Etheostoma caeruleum mitochondrion, complete genome
                        "HQ128171.1", # Etheostoma erythrozonum isolate EmerB voucher YPM 20829 cytochrome b (cytb) gene
                        "NC_011380.1", # Fundulus olivaceus mitochondrion, complete genome
                        "NC_008676.1", # Hypentelium nigricans mitochondrion, complete genome
                        "JN389795.2", # Lepomis macrochirus mitochondrion, complete genome
                        "NC_084073.1", # Lepomis megalotis isolate SMCr-10 mitochondrion, complete genome
                        "NC_080909.1", # Luxilus zonatus mitochondrion, complete genome
                        "NC_011361.1", # Micropterus dolomieu mitochondrion, complete genome
                        "NC_084084.1", # Moxostoma duquesnii isolate LP-2 mitochondrion, complete genome
                        "NC_033924.1", # Nocomis biguttatus mitochondrial DNA, complete genome
                        "NC_080908.1", # Noturus exilis mitochondrion, complete genome
                        verbose=TRUE))
# write.csv(seqs, "NCBI_accessions_distMatrix.csv")

# Align sequences with Clustal Omega
myDNAseqs <- DNAStringSet(seqs[,16]) # Just the column of sequences
seqs$taxon <- sub(" ", "_", seqs$taxon)
names(myDNAseqs) <- seqs$taxon
myAlignment <- msa(myDNAseqs, "ClustalOmega")
aligned_DNA <- as(myAlignment, "DNAStringSet")

# Identify and subset all sequences to cytB region
cytb_seq <- aligned_DNA["Cottus_carolinae"] # species with shortest cytb sequence available
cytb_chars <- unlist(strsplit(as.character(cytb_seq), split = ""))
cytb_positions <- which(cytb_chars %in% c("A", "T", "G", "C"))
cytb_start <- min(cytb_positions) # get start and end positions
cytb_end <- max(cytb_positions)
subset_cytb <- subseq(aligned_DNA, start = cytb_start, end = cytb_end)

# Subset aligned sequences to cytB Region
cytb_bin <- as.DNAbin(subset_cytb)

# Calculate pairwise distance 
cytb_bin <- as.DNAbin(subset_cytb)
dist_cytb <- dist.dna(cytb_bin, model = "TN93", pairwise.deletion = TRUE)
write.csv(as.matrix(dist_cytb), "cytb_distMatrix.csv")

# Build and plot tree
tree_cytb <- nj(dist_cytb)
pdf("phylo_tree.pdf", width=8, height=8)
plot(tree_cytb)
dev.off()
