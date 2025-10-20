# 🧬 _de novo_ genome assembly and annotation       
   
> This tutorial guides you through a **DNA sequencing analysis workflow**, from raw reads to _de novo_ genome assembly.     
   
De novo genome assembly is a cornerstone of genomics that enables researchers to reconstruct an organism’s genome directly from sequencing reads, without relying on a reference. This approach is essential for studying non-model species, exploring genomic novelty, and uncovering structural features that reference-based methods can miss. By overlapping and connecting short or long reads into contiguous sequences (contigs) and scaffolds, assembly algorithms rebuild the genomic landscape from scratch. The resulting assemblies provide the foundation for gene annotation, comparative genomics, and evolutionary analyses, allowing researchers to investigate genome organization, repeat content, and the genetic basis of biological diversity.   
     
---   
   
$\color{Orange}\Huge{\textbf{Major aims}}$   
   
- Understanding the challenges inherent to _de novo_ genome assembly.     
- Appreciate the importance of certain parameters.     
- Generate an assembly and investigate results.
- Understanding the concept of genome annotation and gff/bed file format.    
   
---     
   
#### :rotating_light: No report to submit :rotating_light:     
   
**Tips**: be aware of `where` you currently are on the server (PATH), which system (console==R, terminal==shell) and environment (micromamba activate)   

---     

- Using the same data you copied during the 3rd course, you'll `_de novo_` assemble this yeast genome using the same clean, paired-end reads used for variant calling.    
   
For that you'll be using `SPAdes`. SPAdes (St. Petersburg genome assembler) is a popular tool for de novo genome assembly from short-read sequencing data. It is particularly well-suited for paired-end Illumina reads and uses a de Bruijn graph–based approach with sophisticated error correction to reconstruct contiguous genome sequences (contigs). For more information about the tools you can visit their github [page](https://ablab.github.io/spades/).     

---     

$\color{Green}\{\textbf{To explore the effect of k-mer size on assembly output, we'll split the class into 4 groups:}}$  
  
--> **4 groups that will try different k-mer sizes for scaffolding, talk among yourselves to agree on 4 set sof k-mer sizes to try and define them as in**    
```
k1=x
k2=y
k3=z
```
  
You only need your clean, trimmed set of reads generated ahead (the assembly should take ~10min):     
```   
micromamba activate DNA_env
clean_fwd_reads="BBtrim/ERR1309170_clean1.fq"   
clean_rev_reads="BBtrim/ERR1309170_clean2.fq"   
seq_id="ERR1309170"
spades.py --only-assembler -t 16 -m 252 -k $k1,$k2,$k3 --isolate --tmp-dir /scratch/ -o ${seq_id}_SpadesAssembly -1 $clean_fwd_reads -2 $clean_rev_reads   
```   

---     

Some of the main output files include:     
- contigs.fasta – FASTA file containing the assembled contigs. **These are contiguous sequences without gaps; often the main file you’ll analyze.**     
- scaffolds.fasta – FASTA file containing scaffolds, which are contigs joined together with estimated gaps (usually represented by Ns). Generally more contiguous than contigs.fasta.     
- spades.log – Detailed log file of the run, including steps performed and parameters used.     
- params.txt – The exact parameters used by SPAdes.     
- dataset.info / input_dataset.yaml – Metadata about your input reads and datasets.     
- K21, K33, K55 – Directories containing intermediate assemblies for different k-mer sizes tested during the run.     


---     

**Investigate the resulting files.**     
   
--> **How many scaffolds were assembled during genome assembly?**     
   
Investigate how it compares to the reference genome assembly used for read mapping at previous course.     
   
More or less contigs?     
How does the total genome size compares?  
What impact of k-mer sizes on the contiguity?  
   
> Illustrates the difficulties in assemblies full chromosomes (even for species with relatively small genomes and low repeat content)   
   
**Why assemblies are useful**     
👉 Using the generated scaffolds, you can start looking for genes (i.e. coding sequences), and infer RNA/protein sequences.     

- To get an idea of the types of hints / clues we can use to perform `eukaryote` gene annotation, have a look at this [link](https://neutra.bzh.uni-heidelberg.de/jbrowse/JBrowse-1.12.3/index.html?loc=1%3A4540661..4558420&tracks=DNA%2CWT%2CNC10&highlight=)   

---     

**Annotating a genome assembly**  


> Genome annotation is the process of identifying and describing functional elements within a newly assembled genome, transforming raw sequence data into biologically meaningful information. This involves predicting protein-coding genes, non-coding RNAs, transposable elements, and regulatory regions. Gene annotation combines ab initio prediction—using statistical models trained on known gene structures—with evidence-based approaches that align transcriptomic (RNA-seq) or proteomic data to the genome. In eukaryotes, features such as intron–exon boundaries, open reading frames, and conserved protein domains help refine predictions, while repeat annotation tools identify transposable elements and repetitive DNA that shape genome evolution. Together, these methods yield a comprehensive map of genomic features, providing the basis for understanding gene function, regulation, and evolutionary dynamics.

---

Using the `S288C_reference_genome_Current_Release` folder downloaded previously, search for a file with a `gff.gz` suffix.   


Once you have found it, you can uncompress it with:  
```
gunzip the_actual_file
```

- investigate the GFF file structure, what does it represents?

Now you'll use two aditional tools to extract some biological sequences from this GFF file providing the reference genome sequence:   

- First, remember to activate your micromamba emvironment and then install the required tools  
```
micromamba activate DNA_env
micromamba install gffread
micromamba install emboss
```

Using a GFF parser called `gffread` you'll recover the sequences corresponding to the predicted transcripts with (careful with the provided PATH):   

```
REF=where_your_S288C_reference_sequence_R64-5-1_20240529.fna_file_is_located
GFF=where_your_saccharomyces_cerevisiae_R64-5-1_20240529.gff_file_is_located
output=saccharomyces_cerevisiae_R64-5-1_20240529_transcripts.fna 
gffread -x $output -g $REF $GFF
```

- Look at your output file, how many predicted transcripts?

You can now translate those transcripts into protein sequences using `transeq`:   
```
transeq -sequence saccharomyces_cerevisiae_R64-5-1_20240529_transcripts.fna -outseq saccharomyces_cerevisiae_R64-5-1_20240529_proteins.faa
```

- Look at your output file, how many predicted proteins?


