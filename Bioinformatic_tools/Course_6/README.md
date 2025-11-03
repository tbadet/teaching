# 🧬 Genome-wide association studies [GWAS]: linking variation in phenotypes to genomic variation     
   
> This tutorial guides you through a **GWAS workflow**, using variant calls and various phenotypes on a large sample of yeast strains.     
   
Genome-Wide Association Studies (GWAS) leverage natural **intraspecific** genetic variation across populations to uncover the genetic basis of phenotypic diversity. By scanning millions of genetic markers across the genome and correlating them with measurable traits or disease states, GWAS identifies loci where variation is statistically associated with phenotype differences. This population-level approach provides a powerful framework to link genotype to phenotype, revealing genes, pathways, and regulatory regions that contribute to complex traits, adaptation, and disease susceptibility. As a cornerstone of modern population genomics and medical genetics, GWAS transforms our understanding of how inherited variation shapes function, health, and evolution.  

---   
   
$\color{Orange}\Huge{\textbf{Major aims}}$   
   
- Understanding the major prerequisites to conduct a GWAS.     
- Appreciate the power amd limitations of such studies.     
- Perform the association and investigate results.
     
---     

- We'll be using the extensive dataset generated on >1k yeast strains collected from different environments. Here is the [link](https://www.nature.com/articles/s41586-018-0030-5) to the study first describing the **genomic dataset** and [here](https://www.cell.com/cell-reports/fulltext/S2211-1247(16)30802-6) the study that describes the **phenotypic variation** in the population.   

From the Terminal, you can copy the genotyping (HapMap format) and the phenotyping tables from my Temp folder using:  
```
cp -r /legserv/Temp/Thomas/GWAS_dataset /path/where/you/want/to/have/it
```

You can have a look at those tables with:   
```
head /path/to/pheno_filtered.tab
head /path/to/geno_random_1kb.hmp.txt
```

What is represented in rows and columns?    

**Note**:  
- The Hapmap file format is a table which consists of 11 columns plus one column for each sample genotyped. The first row contains the header labels of your samples, and each additional row contains all the information associated with a single SNP. You can get a Hapmap file by chromosome or a general file.     
     
The current release consists of attributes, as follows:     
     
> rs#    alleles    chrom     pos     strand    assembly#    center    protLSID    assayLSID    panelLSID    QCcode    sample1 |     
rs# contains the SNP identifier;     
alleles contains SNP alleles according to NCBI database dbSNP;     
chrom contains the chromosome that the SNP was mapped;     
pos contains the respective position of this SNP on chromosome;     
strand contains the orientation of the SNP in the DNA strand. Thus, SNPs could be in the forward (+) or in the reverse (-) orientation relative to the reference genome;     
assembly# contains the version of reference sequence assembly (from NCBI);     
center contains the name of genotyping center that produced the genotypes;     
protLSID contains the identifier for HapMap protocol;     
assayLSID contain the identifier HapMap assay used for genotyping;     
panelLSID contains the identifier for panel of individuals genotyped;     
QCcode contains the quality control for all entries;     
subsequently, the list of sample names.     
Below it is an example of VCF file:     
      
As a starter, have a quick look at those two studies. What are the main informations you get from those?   

---    

#### Performing the association mapping   

You'll perform the association in R using the [GAPIT](https://zzlab.net/GAPIT/) package (Genome Association and Prediction Integrated Tool).   

1️⃣ Load necessary packages and GAPIT functions 
```
library(dplyr)
source("http://zzlab.net/GAPIT/gapit_functions.txt")
```

  Purpose:   
    dplyr: for clean data manipulation and subsetting.     
    The source() command loads all GAPIT functions from the public repository.      

Note: GAPIT implements several GWAS models, including MLM (Mixed Linear Model), CMLM, and FarmCPU, with built-in PCA correction for population structure.   

   
2️⃣ Load phenotype data. 
Input: A tab-delimited file containing phenotypes for all accessions/individuals.   
```
pheno <- data.table::fread("pheno_filtered.tab", header = TRUE)
```

3️⃣ Load genotype data   
Input: A HapMap-format genotype file (.hmp.txt), containing SNP genotypes across all individuals   
```
geno <- data.table::fread("geno_random_1kb.hmp.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
geno <- as.data.frame(geno)
```

  Explanation:    
    Each SNP is a row; each column (after metadata) is a genotype for a given individual.   
    The number of taxa in genotype and phenotype files must match.   

- Now select a phenotype for which you want to perform the association. For an idea of what these phenotypes identifiers represent you can have a look at the supplementaries of the [study](https://www.cell.com/cms/10.1016/j.celrep.2016.06.048/attachment/f7cca740-8ef2-42a6-815d-a0e61e132efb/mmc1.pdf)

##### As a example using the `YPACETATE` phenotype (e.g. yeast growing on plate with acetate as unusual carbon source)   
```
trait="YPACETATE"
trait.index <- which(names(pheno)==trait)
var <- pheno %>% dplyr::select(1, trait.index)
var_df <- as.data.frame(var)
rownames(var_df) <- var_df$V1
names(var_df) <- c("Taxa", trait)

# now you can run the association with:
myGAPIT <- GAPIT(
  Y = var_df, # phenotype
  G = geno,
  PCA.total=3
  )
```

GAPIT produces a a large number of output files, here are the files generated with the `YPACETATE` example above:   

>GAPIT.Phenotype.View.YPACETATE.pdf
GAPIT.Association.Filter_GWAS_results.csv
GAPIT.Association.GWAS_Results.MLM.YPACETATE(NYC).csv
GAPIT.Association.GWAS_StdErr.MLM.YPACETATE(NYC).csv
GAPIT.Association.Manhattans_Symphysic_Traitsnames.csv
GAPIT.Association.Prediction_results.MLM.YPACETATE.csv
GAPIT.Association.PVE.MLM.YPACETATE.csv
GAPIT.Association.Vairance_markers.MLM.YPACETATE.csv
GAPIT.Genotype.Distance.Rsquare.csv
GAPIT.Genotype.Frequency_MAF.csv
GAPIT.Genotype.Kin_Zhang.csv
GAPIT.Genotype.PCA_eigenvalues.csv
GAPIT.Genotype.PCA.csv
GAPIT.Association.Manhattan_Chro.MLM.YPACETATE(NYC).pdf
GAPIT.Association.Manhattan_Geno.MLM.YPACETATE(NYC).pdf
GAPIT.Association.Manhattans_Symphysic_Legend.pdf
GAPIT.Association.Manhattans_Symphysic.pdf
GAPIT.Association.Optimum.MLM.YPACETATE.pdf
GAPIT.Association.QQ.MLM.YPACETATE(NYC).pdf
GAPIT.Association.Significant_SNPs.MLM.YPACETATE.pdf
GAPIT.Genotype.Distance_R_Chro.pdf
GAPIT.Genotype.Frequency.pdf
GAPIT.Genotype.Kin_Zhang.pdf
GAPIT.Genotype.MAF_Heterozosity.pdf
GAPIT.Genotype.PCA_2D.pdf
GAPIT.Genotype.PCA_3D.pdf
GAPIT.Genotype.PCA_eigenValue.pdf
GAPIT.Phenotype.Distribution_Significantmarkers.MLM.YPACETATE.pdf


The full statistics outcome of the GWAS is in the file `GAPIT.Association.GWAS_Results.MLM.YPACETATE(NYC).csv`. This file reports for every SNP in the dataset the association with the analyzed trait.   

Let's use R to identify the SNP showing the most significant association with the phenotype. You can graphically identify this SNP also by opening the file `GAPIT.Association.Manhattan_Geno.MLM.YPACETATE(NYC).pdf` (adjust to your trait).   

#### 🧠 Digging deeper into the GWAS results   

- You also should find a file `GAPIT.Association.PVE.MLM.YPACETATE.csv` that list the top associated variants. It contains important output metrics from the applied GWAS model:   

The table has two R² values:  
GAPIT fits a linear mixed model (LMM) for each SNP:  
> y = Xβ + Zu + e   

where:  
> y = Vector of phenotypic values >> The observed trait measurements for all individuals (e.g., plant height, yield, etc.).   
X = Design matrix for fixed effects	Encodes known effects such as intercepts, population structure (e.g., PCs or covariates), and SNP genotype being tested.   
β = Vector of fixed effect coefficients >> The estimated effects corresponding to X (e.g., SNP effect sizes, covariate coefficients).   
Z = Design matrix for random effects >> Links individuals to their genetic background or kinship structure.   
u = Vector of random genetic effects >> Captures the contribution of background polygenic effects assumed to follow u∼N(0,Kσ_g²), where K is the kinship matrix and σ_g² is the genetic variance.   
e = Vector of residuals (errors)	>> Represents environmental noise and unexplained variance, assumed e∼N(0,Iσ_e²).   

When testing one SNP, it compares two models:   

(i) Without SNP (null model):    
Only population structure (PCA) and kinship (K) terms are included.   
> y = Xβ + Zu + e   

→ gives R²_without SNP   
   
(ii) With SNP (full model):   
Adds the tested SNP as a fixed effect.   
> y = Xβ + g_"SNP" + Zu + e   

where `g_SNP` is the fixed effect of the SNP >> The contribution of the SNP being tested to the phenotype. This is usually additive (0, 1, 2 for genotype counts).   
→ gives R²_with SNP   
   
!!! The difference between these two R² values shows how much more variance the model explains by including that SNP:   
> ΔR² = R²_with SNP − R²_without SNP

- How much of the phenotypic variance is explained by the model and how much is the contribution of your top associated variant?    

