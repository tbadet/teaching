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
wget https://raw.githubusercontent.com/tbadet/teaching/refs/heads/main/Bioinformatic_tools/Course_5/datasets/Yeast_transcriptomic_data.csv
```

Then you can extract the sequencing run identifiers for download:  
```
cut -d',' -f1 Yeast_transcriptomic_data.csv | tail -n +2 > srx.txt
sed -i '' 's/"//g' srx.txt
esearch -db sra -query "$QUERY" | efetch -format runinfo | sed 1d | cut -d',' -f1 > srr.txt
cut -d',' -f2 Yeast_transcriptomic_data.csv | tail -n +2 > sra.txt
```

**Explanation**      
  
`cut -d',' -f1` extracts the first column of the CSV file (here, SRX accession numbers).   
`tail -n +2` skips the header line.    
`sed -i '' 's/"//g'` removes any quotation marks from the file.   
`esearch` and `efetch` are NCBI Entrez tools to query the SRA database:   
`-db sra` specifies the database.   
`-query "$QUERY"` uses a query term (e.g., an SRX number).   
`efetch -format runinfo` retrieves run metadata (like SRR IDs, platforms, layout, etc.).   
`sed 1d` removes the header line.   
`cut -d',' -f1` extracts the SRR run IDs.  
  
The final cut extracts the second column (SRA IDs) into a separate file.


### 🧪 2. Select the strain and conditions of interest   
```
grep "YPS606" Yeast_transcriptomic_data.csv | egrep -v "mutant|TF" | egrep "unstressed|NaCl" > get_files.txt
```

**Explanation**    
`grep "YPS606"` filters for rows corresponding to strain YPS606.   
`egrep -v "mutant|TF"` excludes any entries related to mutants or transcription factor experiments.   
`egrep "unstressed|NaCl"` selects the two conditions of interest — unstressed (control) and NaCl (salt stress).   

The filtered lines are written to get_files.txt.    
- Purpose: To define a focused experimental comparison, limiting downloads and analyses to specific biological conditions (here: how strain YPS606 responds to salt stress).    

  
### 💾 3. Download raw reads   

```
export TMPDIR="./"
IFS=','
while read line
do
linearray=( $line )
SRA=${linearray[1]}
TYPE=${linearray[2]}
TREATMENT=$( echo $TYPE | perl -pe 's/; Saccharomyces cerevisiae; RNA-Seq//' | perl -pe 's/^.*: //' | tr ' ' '_' )
echo $TREATMENT
fastq-dump --outdir fastq/ --split-files --gzip $SRA
mv fastq/${SRA}_1.fastq.gz fastq/${TREATMENT}_1.fastq.gz 
mv fastq/${SRA}_2.fastq.gz fastq/${TREATMENT}_2.fastq.gz 
done < get_files.txt
```

**Explanation**   
   
`export TMPDIR="./"` sets a temporary directory for downloads.    
`IFS=','` tells Bash that input fields are comma-separated.    
The while read line loop reads each line from `get_files.txt`   
`linearray=( $line )` splits the line into an array.   
`SRA=${linearray[1]}` stores the SRA run accession (e.g., SRR123456).   
`TYPE=${linearray[2]}` stores a description field (e.g., “RNA-Seq; unstressed”).  
`TREATMENT=$(...)` cleans the text, extracting only the condition (e.g., unstressed or NaCl) and replacing spaces with underscores.   
`fastq-dump` downloads sequencing reads from SRA:   
`--outdir fastq/` outputs into the fastq folder.   
`--split-files` splits paired-end reads into _1 and _2.   
`--gzip` compresses the FASTQ files.   
Files are renamed to meaningful names based on the treatment.    
    
Purpose: To download and organize raw FASTQ reads in a structured, condition-specific manner.  
  
### 🧬 4. Download and prepare the reference genome   
```
wget http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_Current_Release.tgz
gunzip S288C_reference_genome_Current_Release.tgz
tar xvf S288C_reference_genome_Current_Release.tar
zcat S288C_reference_genome_R64-5-1_20240529/S288C_reference_sequence_R64-5-1_20240529.fsa.gz | perl -pe 's/>.*chromosome=/>chr/' | tr -d ']' | perl -pe 's/>ref.*/>chrmt/' > S288C_reference_sequence_R64-5-1_20240529.fna
```
  
**Explanation**  
`wget` downloads the current S. cerevisiae S288C genome release.  
`gunzip` and `tar` extract the archive.   
`zcat` decompresses the `.fsa.gz` FASTA file.    

The series of `perl` and `tr` commands reformat the FASTA headers:   
Replace long headers with short ones like >chrI, >chrII, etc.   
Clean up any bracket or reference annotations.   
Rename the mitochondrial chromosome as `chrmt`.   

Purpose: To standardize chromosome names and prepare the genome FASTA file for downstream tools (STAR requires simple, unique headers).   
  
### 🧾 5. Convert genome annotation from GFF3 to GTF
```
gffread saccharomyces_cerevisiae_R64-5-1_20240529.gff -T -o saccharomyces_cerevisiae_R64-5-1_20240529.gtf
```

**Explanation**  
`gffread` (from the Cufflinks suite) converts GFF3 → GTF.  
`-T` specifies GTF output format.  
The output file is written with -o.  
  
Purpose: STAR requires a GTF annotation (not GFF3) for indexing and spliced alignment.   

###  🧩 6. Build STAR genome index   
```
STAR --runThreadN 4 --runMode genomeGenerate \
--genomeSAindexNbases 12 \
--genomeDir S288C_index \
--genomeFastaFiles S288C_reference_sequence_R64-5-1_20240529.fna \
--sjdbGTFfile saccharomyces_cerevisiae_R64-5-1_20240529.gtf \
--sjdbOverhang 74
Explanation:
--runMode genomeGenerate tells STAR to build an index.
--runThreadN 16 uses 16 threads (parallel processing).
--genomeDir specifies output directory for the index.
--genomeFastaFiles provides the reference FASTA.
--sjdbGTFfile adds splice junction annotations from the GTF file.
--sjdbOverhang 74 should equal read length minus 1, to optimize splice junction mapping.
```
  
Purpose: To prepare the genome index needed for fast, splice-aware RNA-seq alignment.

###  🧮 7. Read trimming, alignment, and quantification
```
for i in fastq/*_1.fastq.gz
do
resource_dir="/data2/toto/softs/bbmap/resources"
dir="/data2/toto/softs/bbmap"
gtf="saccharomyces_cerevisiae_R64-5-1_20240529.gff"
i2=$(echo ${i%_1.fastq.gz} | perl -pe 's/^.*\///' )
echo $i2
zcat $i > ${i%_1.fastq.gz}_1.fq
$dir/bbduk.sh threads=20 -Xmx20g in1=${i%_1.fastq.gz}_1.fq out1=BBtrim/${i2}_clean1.fq ref=$resource_dir/polyA.fa.gz,$resource_dir/truseq.fa.gz ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10 minlength=20
rm ${i%_1.fastq.gz}_1.fq
STAR --runMode alignReads \
--genomeDir S288C_index \
--runThreadN 4 \
--readFilesIn BBtrim/${i2}_clean1.fq \
--outFilterMultimapNmax 100  \
--winAnchorMultimapNmax 200  \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI NM MD \
--outFileNamePrefix BAM/${i2}_multi
samtools index -@ 4 BAM/${i2}_multiAligned.sortedByCoord.out.bam
htseq-count BAM/${i2}_multiAligned.sortedByCoord.out.bam $gtf -n 4 -r pos -f bam -s no -t CDS -i Parent -m union > ${i2}_HTSeq.txt
done
```

**Explanation**  

- 🔁 Loop & variables:    
  
The `for` loop iterates over all first-read FASTQ files (*_1.fastq.gz).
`i2` extracts the basename (e.g., unstressed_rep1).
`resource_dir` and `dir` specify the BBMap software and resources (adapters, etc.).
`gtf` defines the gene annotation used for counting.

🧹 Read trimming:    

`zcat` decompresses FASTQ.  
`bbduk.sh` (BBMap toolkit) trims and cleans reads:  
`threads=20 -Xmx20g` use 20 threads, 20 GB RAM.    
`ref=` points to adapter sequences (polyA, Illumina TruSeq).   
`ktrim=r` trim from the right end.    
`k=23 mink=11 hdist=1` k-mer size 23, minimum 11, allowing 1 mismatch.    
`tpe tbo`trim both reads consistently and based on overlap.   
`qtrim=r trimq=10` quality-trim right end, cutoff Q10.   
`minlength=20` discard reads shorter than 20 nt.    
  
 - 🎯 Alignment with STAR:
   
`--genomeDir`: use previously built index.  
`--runThreadN 4`: 4 threads.   
`--outFilterMultimapNmax 100`: allow up to 100 multiple alignments per read (for repetitive genes).   
`--winAnchorMultimapNmax 200`: increases sensitivity in repetitive regions.   
`--outSAMtype BAM SortedByCoordinate`: output sorted BAM.   
`--outSAMattributes NH HI NM MD`: include mapping statistics (multi-mapping count, edit distance, etc.).   
   
Note that output files are named `BAM/<sample>_multiAligned.sortedByCoord.out.bam`.

- 🧾 Quantification:
  
`samtools index` creates BAM index files.   
`htseq-count` counts reads overlapping genes:   
`-n 4`: 4 threads.  
`-r pos`: BAM is sorted by position.   
`-f bam`: input format.   
`-s no`: library is unstranded (reads not strand-specific).   
`-t CDS`: count only reads in coding sequences.   
`-i Parent`: use the Parent attribute in GFF (links CDS to mRNA/gene).   
`-m union`: count a read if it overlaps any CDS of a gene (most common mode).   
    
Purpose: To clean, align, and count reads per gene — generating raw counts for downstream DESeq2 analysis.    




