## 🧬 DNA Sequencing Data Analysis: Read Mapping & Genome Assembly

This tutorial guides you through a **DNA sequencing analysis workflow**, from raw reads to variant calling and de novo genome assembly.

---

**Learning objectives:**

- Understand the structure of sequencing metadata (ERX/ERR, library, instrument).  
- Preprocess reads by trimming adapters and low-quality bases.  
- Map reads to a reference genome using **Bowtie2**.  
- Identify SNPs and small indels using **GATK**.  
- Filter variants to retain high-confidence calls.  
- Perform **de novo genome assembly** with SPAdes.

---

**Key concepts:**  

- **Read trimming** reduces artifacts from sequencing adapters and low-quality bases, improving mapping accuracy.  
- **Reference indexing** allows mapping tools to efficiently locate sequences in the genome.  
- **Variant calling** identifies differences between the sample and reference genome.  
- **Variant filtering** ensures only reliable variants are used for downstream analyses.  
- **De novo assembly** reconstructs the genome without relying on a reference, useful for novel strains.

---

## Dataset Overview

- Metadata: `ERR.csv` → includes raw run IDs (ERR) from DNA sequencing experiments available on NCBI.  
- Reference genome: **S. cerevisiae S288C**  
- Adapter sequences: `truseq.fa.gz` (Illumina TruSeq adapters)




