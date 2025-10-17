# 🧬 DNA Sequencing Data Analysis: _de novo_ genome assembly      
   
> This tutorial guides you through a **DNA sequencing analysis workflow**, from raw reads to _de novo_ genome assembly.     
   
Mapping sequencing reads to a reference genome and calling variants is a central method in modern biology for understanding genetic differences and their functional consequences. By aligning raw reads to a reference, researchers can pinpoint where in the genome each fragment originates, which allows the detection of single nucleotide variants, insertions, deletions, and larger structural changes. This approach is useful in many contexts, from identifying mutations linked to diseases, to studying genetic diversity in populations, to exploring evolutionary relationships between species. Once variants are identified, they can be connected to phenotypes, adaptive traits, or molecular mechanisms, making this pipeline a foundation for both basic research and applied genomics.     
   
---   
   
$\color{Orange}\Huge{\textbf{Major aims}}$   
   
- Understanding the challenges inherent to _de novo_ genome assembly.     
- Appreciate the importance of certain parameters.     
- Generate an assembly and investigate results.     
   
---     
   
#### :rotating_light: No report to submit :rotating_light:     
   
   
- Using the same data you copied during the 3rd course, you'll `_de novo_` assemble this yeast genome using the same clean, paired-end reads used for variant calling.    
   
For that you'll be using `SPAdes`. SPAdes (St. Petersburg genome assembler) is a popular tool for de novo genome assembly from short-read sequencing data. It is particularly well-suited for paired-end Illumina reads and uses a de Bruijn graph–based approach with sophisticated error correction to reconstruct contiguous genome sequences (contigs). For more information about the tools you can visit their github [page](https://ablab.github.io/spades/).     

$\color{Green}\{\textbf{To explore the effect of `_k-mer_` size on assembly output, we'll split the class into 4 groups:}}$  
  
--> **4 groups that will try different k-mer sizes for scaffolding, talk among yourselves to agree on 4 set sof k-mer sizes to try and define them as in**    
```
k1=x
k2=y
k3=z
```
  
You only need your clean, trimmed set of reads generated ahead (the assembly should take ~10min):     
```   
clean_fwd_reads="BBtrim/ERR1309170_clean1.fq"   
clean_rev_reads="BBtrim/ERR1309170_clean2.fq"   
seq_id="ERR1309170"
spades.py --only-assembler -t 16 -m 252 -k $k1,$k2,$k3 --isolate --tmp-dir /scratch/ -o ${seq_id}_SpadesAssembly -1 $clean_fwd_reads -2 $clean_rev_reads   
```   
   
Some of the main output files include:     
- contigs.fasta – FASTA file containing the assembled contigs. **These are contiguous sequences without gaps; often the main file you’ll analyze.**     
- scaffolds.fasta – FASTA file containing scaffolds, which are contigs joined together with estimated gaps (usually represented by Ns). Generally more contiguous than contigs.fasta.     
- spades.log – Detailed log file of the run, including steps performed and parameters used.     
- params.txt – The exact parameters used by SPAdes.     
- dataset.info / input_dataset.yaml – Metadata about your input reads and datasets.     
- K21, K33, K55 – Directories containing intermediate assemblies for different k-mer sizes tested during the run.     
   
**Investigate the resulting files.**     
   
--> **How many scaffolds were assembled during genome assembly?**     
   
Investigate how it compares to the reference genome assembly used for read mapping at previous course.     
   
More or less contigs?     
How does the total genome size compares?  
What impact of k-mer sizes on the contiguity?  
   
> Illustrates the difficulties in assemblies full chromosomes (even for species with relatively small genomes and low repeat content)   
   
**Why assemblies are useful**     
👉 Using the generated scaffolds, you can start looking for genes (i.e. coding sequences), and infer RNA/protein sequences.     
