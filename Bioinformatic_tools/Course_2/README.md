## 🧬 DNA Sequencing Data Analysis: Read Mapping & Genome Assembly

This tutorial guides you through a **DNA sequencing analysis workflow**, from raw reads to variant calling and de novo genome assembly.

---

$\color{Orange}\Huge{\textbf{Major aims}}$

- Navigate sequencing metadata (ERX/ERR, library, instrument).  
- Understand the structure of sequencing data.  
- Generate DNA alignemnts from DNA sequencing datasets.  
- Lear to apprehend different associated file formats.  

---

#### :rotating_light: Report to submit :rotating_light:  

Please compile brief answers to the questions for your report ("Q1", "Q2", etc.).  
Depending on the question, you can add a short text, some code or a graphic. The answers can be in a text file, Word doc, etc.  

You can work alone or in groups. Every student should submit their own report through Moodle though. No copy-pasting, please. Formulate answers in your own words.  

---------------------------------------------------------------------------------------   

:reminder_ribbon: $\color{Orange}\Large{\textbf{General troubleshooting tips}}$ 

Ask yourself whether you are in the correct folder / location. Try to be consistent and make some notes in your script file. Use pwd to check where you are and cd if necessary.  

If you want to check out e.g. the command cp. Googling for "Linux cp" tells you fairly well what to do. Try for yourself how to best Google/chatGPT/etc. your answers.  

Read the (error) message that you get in the Terminal following your commands. You may not understand these in full, but it's a good sign that you should expect something strange (a missing file, etc.).  

Google/chatGPT the error message. Copy the error message and search for it. Helpful? Try to improve your code googling skills during the course.  

If you can't find a file that you have downloaded or created, start by checking with pwd where you are. Then check whether you might have saved the file in a different folder? Was there any error message about a location or file not found?  

Ask chatGPT (or similar) to write code for you. Does it work? For fairly simple tasks, it gets it often right. However, where chatGPT will fail is if it does not know everything about e.g. your files or folders. Tell it specifically what to do then. 

---------------------------------------------------------------------------------------   

#### (1) Install dedicated softwares  

The first step in a bioinformatics project often includes installing software. Past week we learned how to create a dedicated environment using ```micromamba```. Using environments avoids version conflicts and ensures reproducibility across systems.  

Using `wget` you can download a file with a list of softwares needed for today's exercice:   

```
wget https://raw.githubusercontent.com/tbadet/teaching/refs/heads/main/Bioinformatic_tools/Course_2/DNA_env.yaml
```

Try using `cat` or `nano` / `less` to see the programs that will be installed.  

You can then create a new environment using:  

```
micromamba clean --all --yes
micromamba create -f DNA_env.yaml --channel-priority flexible
```

and remermber to activate it:  
```
micromamba activate DNA_env
```

#### (2) Retrieve datasets  

Previous weeks you've seen how to access genome sequences from the NCBI database. Raw sequecing data can also be retrieved from the 


#### Dataset Overview

- Metadata: `ERR.csv` → includes raw run IDs (ERR) from DNA sequencing experiments available on NCBI.  
- Reference genome: **S. cerevisiae S288C**  
- Adapter sequences: `truseq.fa.gz` (Illumina TruSeq adapters)




