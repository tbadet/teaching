# Orthology inference and comparative genomics: uncovering the evolutionary relationships among genes  
> This tutorial helps you assign orthologous relationships across the full set of proteins from different Saccharomyces species and begin exploring the corresponding results.  

Orthology inference seeks to identify genes across different species or individuals that share a common ancestral origin, providing a foundation for understanding genome evolution and functional diversification.  

By comparing gene sequences and reconstructing their evolutionary histories, orthology inference distinguishes orthologs (genes diverged through speciation) from paralogs (genes diverged through duplication). This distinction is critical for comparative genomics, where the goal is to separate conservation from innovation, revealing which genes retain ancestral functions and which have adapted or diversified across lineages.  

In the era of pangenomics, orthology-based analyses allow researchers to compare gene content across many genomes simultaneously, defining the core genome shared by all members of a species and the accessory genome that contributes to phenotypic diversity and ecological adaptation.  

By reconstructing networks of gene relationships, orthology inference enables studies of gene family evolution, horizontal gene transfer, and the emergence of new functions, linking molecular evolution to adaptation and speciation. As a central framework in evolutionary and functional genomics, orthology inference bridges sequence similarity with biological meaning—transforming raw genomic data into evolutionary insights about how genes, pathways, and species evolve.  

----------------------------------------

First, in your shell terminal, activate micromamba and install `orthofinder` with:  
```
micromamba activate
micromamba install orthofinder
```

For more information about the tool, visit this [link](https://github.com/davidemms/OrthoFinder)   

- Next, copy the data from this Github folder:
```
wget https://github.com/tbadet/teaching/raw/refs/heads/main/Bioinformatic_tools/Course_7/yeast_proteins.zip
```

- You first need to unzip the file:  
```
unzip yeast_proteins.zip 
```

**Have a look at what's inside the folder (number of files, what type ..)**  

----------------------------------------

## This next step, not to overload the cluster, decide on a few groups that will run it (you can then share the results)  

- Now, you can run the `orthofinder` pipeline by calling:  
```
orthofinder.py -a 1 -t 1 -f /path/to/yeast_proteins
```

You will see that `orthofinder` produces a global directory named `Orthofinder` inside of which you'll find many output folders and files.

Key Orthofinder output files and folders  

> **1. Citation.txt**   
> Contains references to the original Orthofinder papers and related tools. Useful for acknowledging the methodology when using results in publications.  
> 
> **2. Comparative_Genomics_Statistics**  
> Summary statistics across all genomes analyzed. Includes metrics like total number of genes, number of orthogroups, proportion of single-copy orthologs, and gene duplication events. Gives a high-level overview of gene content and conservation patterns.  
> 
> **3. Gene_Duplication_Events**  
> Lists all duplication events inferred in each gene family. Helps identify gene family expansions, lineage-specific duplications, or potential paralogs.  
> 
> **4. Gene_Trees**
> Contains phylogenetic trees for each orthogroup (typically in Newick format).
> Trees can be visualized to inspect evolutionary relationships among genes and detect orthologs versus paralogs.
> 
> **5. Log.txt**  
> Execution log for the run. Useful for troubleshooting or verifying parameters and runtime details.
> 
> **6. Orthogroups**  
> Contains tables of orthogroups (sets of genes descended from a single ancestral gene).  
> - Includes:  
> Orthogroups.tsv → all genes grouped by orthogroup.  
> Orthogroups_UnassignedGenes.tsv → genes not assigned to any orthogroup.  
> Key file for analyzing core vs. accessory genome content in pangenome studies. 
> 
> **7. Orthogroup_Sequences** 
> FASTA sequences for each orthogroup. Enables downstream analyses like multiple sequence alignment, Ka/Ks calculations, or domain analysis.  
> 
> **8. Orthologues**  
> Lists pairwise orthologous relationships between genes in different species. Useful for targeted comparative studies between two species or for mapping genes to a reference genome.  
> 
> **9. Phylogenetic_Hierarchical_Orthogroups**  
> Hierarchical grouping of orthogroups based on phylogenetic relationships. Useful for tracking gene family evolution across clades.  
> 
> **10. Putative_Xenologs**  
> Genes that may have been acquired via horizontal gene transfer. Helps detect non-vertical inheritance in comparative genomics studies.  
> 
> **11. Resolved_Gene_Trees**  
> Refined gene trees with inferred duplication and speciation events annotated. Important for accurate inference of orthology and paralogy.  
> 
> **12. Single_Copy_Orthologue_Sequences**  
> FASTA sequences of genes present as single copies in all genomes.  
> Often used for: Species tree reconstruction or Comparative analyses that require one-to-one orthologs.  
> 
> **13. Species_Tree**  
> Reconstructed phylogenetic tree of all species included in the analysis. Built using single-copy orthologs and can be used for evolutionary or comparative studies.  
> 
> **14. Phylogenetically_Misplaced_Genes**  
> Genes that do not conform to expected evolutionary relationships. May indicate misannotations, horizontal transfer, or sequencing errors.  

**Explore the resulting files**    

----------------------------------------

- An important output file is in the `Species_Tree` folder:

You can easily visualize it in R using:  
```
library(ggtree)
tree <- read.tree("/path/to/SpeciesTree_rooted.txt")
ggtree(tree) + geom_tiplab()
```

**What can you say?**    

----------------------------------------

#### This next part takes Orthofinder's `Orthogroups.tsv` output and explores it using R to produce:  

- A presence/absence matrix of orthogroups across strains/species.  
- A heatmap of orthogroup copy numbers (capped at 10 for visualization).  
- A summary of species-specific orthogroups.  
- A pairwise sharing table between species.  


##### 1\. Load required R packages
```
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(ggtree)
library(tibble)
```
Explanation:  
These are standard R libraries for data wrangling (dplyr, tidyr, tibble), and visualization (pheatmap, ggplot2, ggtree).

##### 2\. Read and reshape the Orthofinder table
```
df <- read.delim("/path/to/Orthogroups.txt", h = F, sep = " ")
dt <- reshape2::melt(df, id.vars = "V1")
dt$variable <- NULL
dt <- subset(dt, value != "")
dt$strain <- gsub("_.*$", "", dt$value)
dt$V1 <- gsub(":", "", dt$V1)
names(dt) <- c("orthogroup", "protein", "strain")
```
Explanation:   
- Orthogroups.txt is a table where each row = orthogroup, each column = species/strain, and cells = protein IDs.  
- melt() converts the wide table into a long format (tidy data).  
- Empty cells are removed (subset(dt, value != "")).  
- Strain names are simplified by removing anything after an underscore.  
- The colon in orthogroup names (e.g., "OG0001:") is removed for cleaner plotting.  

✅ Resulting data frame:  
> orthogroup | protein | strain  
> OG000001 | XP_12345 | StrainA  
> OG000001 | XP_67890 | StrainB  
> ...  ...  ...  

##### 3\. Assign species names   
```
dt$species <- with(dt, ifelse(strain %in% c("BJ4", "Y55", "HN1", "IMX2600", "S288C", "SX2"), "Scerevisease", strain))
```
Explanation:   
- Some strains belong to the same species (e.g., *Saccharomyces cerevisiae*). This step maps multiple strains to a single species name for species-level analyses.   

##### 4\. Compute orthogroup species representation   
```
ortho_stats <- dt %>%
  group_by(orthogroup) %>%
  summarise(n_species = n_distinct(species))
```
Explanation:  
- Counts how many distinct species are represented in each orthogroup. This helps detect core, shared, and species-specific orthogroups.

✅ Example output:  
> orthogroup | n_species  
> OG000001 | 6  
> OG000002 | 1  
> OG000003 | 3  

##### 5\. Identify species-specific orthogroups   
```
species_specific <- dt %>%
  group_by(orthogroup) %>%
  summarise(n_species = n_distinct(species),
            species = paste(unique(species), collapse = ",")) %>%
  filter(n_species == 1)
```
Explanation:  
- Filters orthogroups present in only one species → potential novel or lineage-specific genes.  

✅ Example output:  
> orthogroup | n_species | species  
> OG004500 | 1 | Scerevisease  
> OG012300 | 1 | Kluyveromyces  

##### 6\. Compute pairwise sharing of orthogroups   
```
shared <- dt %>%
  distinct(orthogroup, species) %>%
  inner_join(., ., by = "orthogroup") %>%
  filter(species.x < species.y) %>%
  count(species.x, species.y, name = "n_shared")
```
Explanation:
- For each orthogroup, all pairs of species that share it are generated.  
- The count of shared orthogroups per pair is computed.  
- This is the pairwise similarity matrix (can be plotted as a heatmap or network).  

✅ Example:  
> species.x | species.y | n_shared  
> Scerevisease | Kluyveromyces | 5243  
> Scerevisease | Lachancea | 4120  

📊 Visual suggestion:  
**Try to visualize it as a heatmap of `n_shared` values.**  

##### 7\. Build a presence/absence matrix
```
presence_matrix <- dt %>%
  group_by(orthogroup, strain) %>%
  summarise(present = 1, .groups = "drop") %>%
  pivot_wider(names_from = strain, values_from = present, values_fill = 0)
```
Explanation:  
- Creates a binary matrix of orthogroup presence (1) or absence (0) per strain.  

✅ Example:  
> orthogroup | StrainA | StrainB | StrainC  
> OG000100 | 1 | 1 | 1  
> OG000200 | 0 | 1 | 0  
> OG000300 | 1 | 0 | 0  

##### 8\. Build a copy-number matrix  
```
count_matrix <- dt %>%
  group_by(orthogroup, strain) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = strain, values_from = count, values_fill = 0)
```
Explanation:  
- Counts how many proteins from each strain belong to each orthogroup. This allows visualization of gene family expansions.  

✅ Example:  
> orthogroup | StrainA | StrainB | StrainC  
> OG000100 | 2 | 1 | 1  
> OG000200 | 0 | 3 | 0  
> OG000300 | 1 | 0 | 0  

##### 9\. Prepare matrices for plotting   
- for the presence / absence matrix  
```
mat1 <- presence_matrix %>%
  tibble::column_to_rownames("orthogroup") %>%
  as.matrix()

mat1[is.na(mat1)] <- 0
```
- for the count matrix
```
mat2 <- count_matrix %>%
  tibble::column_to_rownames("orthogroup") %>%
  as.matrix()

mat2[mat2 > 10] <- 10  # cap high expansions
```
Explanation:  
- Convert data frames to matrices for pheatmap.
- Cap high values at 10 to prevent large gene families from dominating color scales.

##### 10\. Visualize orthogroup counts (copy numbers)
```
col_fun <- colorRampPalette(c("white", "darkgreen", "orange", "red", "darkred", "purple"))(10)
breaks <- seq(0, 10, length.out = 11)

pheatmap(mat2, legend_breaks = breaks, 
         color = col_fun,  
         breaks = breaks,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         border_color = NA,
         main = "Orthogroup counts (capped at 10)")
```
🎨 Interpretation:   
- Each cell = number of gene copies in a strain.  
- Columns (strains) are clustered to show similarity in gene family composition.  

📊 **Explore additional visualization / summary ideas**  

##### 11\. Visualize orthogroup presence/absence
```
pheatmap(mat1, legend_breaks = c(0,1), 
         scale = "none",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         border_color = NA,
         fontsize = 10,
         main = "Presence/Absence of Orthogroups")
```
🎨 Interpretation:   
- 1 = orthogroup present  
- 0 = orthogroup absent  

**Columns (strains/species) clustered by orthogroup composition → gives a pan-genome map.**  

----------------------------------------

