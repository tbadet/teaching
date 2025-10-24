# 🧬 Transcriptomic analysis: estimating differential gene expression      
   
> This tutorial guides you through a **RNA sequencing analysis workflow**, from raw reads to differential gene expression analysis.     
   
RNA sequencing (RNA-seq) and differential gene expression (DEG) analysis provide a comprehensive view of how genes are transcriptionally regulated across conditions, tissues, or developmental stages. By sequencing and quantifying RNA molecules, RNA-seq captures the dynamic landscape of gene activity, enabling comparisons between biological states such as treated versus control, healthy versus diseased, or wild type versus mutant. Through computational statistical analysis we aim to identify genes with significant expression changes, revealing underlying regulatory mechanisms and functional pathways. This approach has become a cornerstone of modern genomics, linking sequence information to cellular function and offering insights into how organisms respond and adapt at the molecular level.       

---   
   
$\color{Orange}\Huge{\textbf{Major aims}}$   
   
- Understanding the challenges inherent to gene expression quantification.     
- Appreciate the importance of the experimental design.     
- Perform differential gene expression analysis and investigate results.
     
---     
   
#### :rotating_light: No report to submit :rotating_light:     
   
**Tips**: be aware of `where` you currently are on the server (PATH), which system (console==R, terminal==shell) and environment (micromamba activate)   

---     

> $\color{Green}{\textbf{The following part illustrates the different steps from raw reads to the reads count table}}$

> $\color{Red}\Huge{\textbf{Do NOT run it}}$  - they are simply here for reference   

### 1. Retrieve metadata and accession lists   
  
You can get the `Yeast_transcriptomic_data.csv` table with the list of RNA-seq experiments and metadata from this github:  
```

```


```
cut -d',' -f1 Yeast_transcriptomic_data.csv | tail -n +2 > srx.txt
sed -i '' 's/"//g' srx.txt
esearch -db sra -query "$QUERY" | efetch -format runinfo | sed 1d | cut -d',' -f1 > srr.txt
cut -d',' -f2 Yeast_transcriptomic_data.csv | tail -n +2 > sra.txt
```





