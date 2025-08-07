# 18S-tree

### 1. Install.packages ####
library(phyloseq)
library(dplyr)
library(ggplot2)
library(vegan)
library(ape)
library(Biostrings)

setwd("D:/Marine_Iguanas_Project/MARINE_IGUANAS/SECOND PAPER 2025")
list.files()

### convert csv with sequences into fasta if not done already
library(tidyverse)
library(readr)
list.files()
csv = read_csv("18S_OTUs_Seqs.csv")
writeLines(paste0(">", csv$zOTU, "\n", csv$sequence), "18S_OTUs_Seqs.fasta")

### 2. Import data ####
otu_mat    <- read.table("OTUtable532samples_18S_GPT_01082025.txt", 
                        header = TRUE, row.names = 1, sep = "\t", check.names=FALSE)
dim(otu_mat)
meta_df   <- read.table("metadata_18S_GPT_01082025.txt", 
                        header = TRUE, row.names = 1, sep = "\t", check.names=FALSE)
dim(meta_df)
tax_df      <- read.table("taxonomy_18S_GPT_01082025.txt",   # your taxonomy file
                          header = TRUE, row.names = 1, sep = "\t", check.names=FALSE)
dim(tax_df)

seqs <- readDNAStringSet("18S_OTUs_Seqs.fasta", format="fasta")

### Convert to phyloseq components ####
OTU <- otu_table(as.matrix(otu_mat), taxa_are_rows = TRUE)
SAM <- sample_data(meta_df)
TAX <- tax_table(as.matrix(tax_df))
REFSEQ <- refseq(seqs)

# Create physeq, first without sequences to merge samples
physeq <- phyloseq(OTU, TAX, SAM)
physeq

### 3.Identify and remove low-replicate Locations ####
# Build a small df of sample → Location
sample_info <- data.frame(
  SampleID = sample_names(physeq),
  Location = sample_data(physeq)$Location,
  stringsAsFactors = FALSE
)
# Count samples per Location
loc_counts <- sample_info %>%
  count(Location, name = "n_samples")

# Which Locations have ≥5 samples?
keep_locs <- loc_counts %>%
  filter(n_samples >= 5) %>%
  pull(Location)

# Prune phyloseq to keep only those samples
physeq_filt <- subset_samples(physeq, Location %in% keep_locs)

# (Optional) drop any taxa that now have zero counts
physeq <- prune_taxa(taxa_sums(physeq_filt) > 0, physeq_filt)
physeq
 
### 3. MERGE ALL SAMPLES FROM EACH ISLAND TOGETHER ####

phy_island <- merge_samples(physeq, "Island", fun = sum)
sample_data(phy_island) <- data.frame(
  Island = sample_names(phy_island),
  row.names = sample_names(phy_island)
)
phy_island

### 4. Align reference sequences #####

# Incorporate sequences in physeq element
physeq1 = merge_phyloseq(phy_island, REFSEQ)
physeq1

# Load DECIPHER
if (!requireNamespace("DECIPHER", quietly=TRUE)) {
  BiocManager::install("DECIPHER")
}
library(DECIPHER)

# Extract sequences and align
seqs <- refseq(phy_island) #changed from: seqs <- refseq(physeq)
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

# Quick check
alignment[1:5]  # print first 5 aligned sequences

### Step 5: infer a phylogenetic tree with phangorn ####
# 5a. Load phangorn
if (!requireNamespace("phangorn", quietly=TRUE)) {
  install.packages("phangorn")
}
library(phangorn)

# 5b. Convert the DECIPHER alignment into a phyDat object
phang.align <- phyDat(as.matrix(alignment), type = "DNA")

# 5c. Estimate a substitution model distance matrix (GTR by default)
dm <- dist.ml(phang.align)

# 5d. Build a neighbor-joining (NJ) tree
treeNJ <- NJ(dm)

# 5e. Optionally midpoint‐root the tree
treeNJ <- midpoint(treeNJ)

# 5f. (Optional) perform a quick ML optimization
fit  <- pml(treeNJ, data = phang.align)
fitGTR <- update(fit, k = 4, inv = 0.2)
fitGTR <- optim.pml(fitGTR,
                    model = "GTR",
                    optInv = TRUE,
                    optGamma = TRUE,
                    rearrangement = "stochastic")

# 5g. Extract the optimized tree
final_tree <- fitGTR$tree

# 5h. Merge into your phyloseq object
physeq_tree <- merge_phyloseq(phy_island, phy_tree(final_tree))

# 5i. Quick plot
plot_tree(physeq_tree,
          ladderize  = "left",
          label.tips = "Phylum",    # or "taxa_names" for OTU IDs
          color = "Island",
          size = "abundance",
          title = "18S OTU Phylogeny (NJ + ML optimized)")

### collapse OTUs at Phylum level ####
phy_phylum <- tax_glom(phy_island, taxrank = "Phylum")
phy_subphy <- tax_glom(phy_island, taxrank = "Subphylum")

phy_phylum <- merge_phyloseq(phy_phylum, phy_tree(phy_tree(physeq_tree)))
plot_tree(phy_phylum,
          ladderize   = "left",
          label.tips  = "Phylum",
          color       = "Island",
          size = "abundance",
          title = "18S OTU Phylum tree with samples by island (NJ)")

phy_subphy <- merge_phyloseq(phy_subphy, phy_tree(phy_tree(physeq_tree)))
plot_tree(phy_subphy,
          ladderize   = "left",
          label.tips  = "Subphylum",
          color       = "Island",
          size = "abundance",
          title = "18S OTU Subphylum tree with samples by island (NJ)")

### Plot ABUNDANCES per taxa per island ####
