----------------------------------------

#### This script takes Orthofinder's Orthogroups.tsv (or a simplified version like Orthogroups_2.txt) and produces using R:  

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
df <- read.delim("~/Documents/Orthogroups_2.txt", h = F)
dt <- reshape2::melt(df, id.vars = "V1")
dt$variable <- NULL
dt <- subset(dt, value != "")
dt$strain <- gsub("_.*$", "", dt$value)
dt$V1 <- gsub(":", "", dt$V1)
names(dt) <- c("orthogroup", "protein", "strain")
```
Explanation:   
- Orthogroups_2.txt is a table where each row = orthogroup, each column = species/strain, and cells = protein IDs.  
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
- Some strains belong to the same species (e.g., Saccharomyces cerevisiae). This step maps multiple strains to a single species name for species-level analyses.   

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
**Try to visualize it as a heatmap of `n_shared` values or a network graph showing species as nodes, connected by edge width ∝ number of shared orthogroups.**  

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

