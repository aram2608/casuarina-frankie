# casuarina-frankie
This is a pipeline for analyzing 100 base-pair single-end Illumina reads for an experiment involving the actinorhizal plant, *Casuarina glauca*.

# Project Members
- Javier
- Matt

# Software Used 💻

- `FastQC v. 0.12.1`
- `Trimmomatic v. 0.39`
- `HISAT2 v. 2.2.1`
- `samtools v. 1.18`
- `BRAKER2 v. 2.1.6`
- `AUGUSTUS v. 3.4.0`
- `GeneMark-ET v. 4.72_lic`
- `gffread v. 0.12.7`
- `BUSCO v. 5.8.2`
- `featurecounts v. 2.2.3`
- `DESeq2 v. 3.20`

# Pipeline Schematic

<img width="751" alt="Screenshot 2025-04-12 at 11 50 26 AM" src="https://github.com/user-attachments/assets/9fff82d4-d821-49c8-b585-5d195fbdd577" />
<img width="748" alt="Screenshot 2025-04-12 at 11 50 34 AM" src="https://github.com/user-attachments/assets/6ee20e14-a602-47ea-9c75-947cf702edae" />

# Introduction

The data used for this study originated from the French National Research Institute for Sustainable Development (IRD). An experiment was conducted where *Casuarina glauca* roots were exposed to either *Frankia* cells or a purified putative signaling molecule (NINA), that we believe is involved in establishing the actinorhizal symbiosis. The plants were exposed for either 24 o4 48 hours respectively. Following exposure, mRNA was extracted and then sequenced. A control plant with neither *Frankia* or NINA exposure was included as well.

The goal of this study is to determine the genes involved in the early interactions of the actinorhizal symbiosis. Owing to the recent assembly of the *Casuarina glauca* genome, a reference guided approach was utilized to map reads for subsequent differential gene expression analysis. However, given the recent status of the reference genome, no annotation file is currently available on the NCBI website. Thus, an *ab initio* gene prediction approach was used to generate a gene model for differential expression analysis using mapped RNA-seq reads as evidence for training. 

We hope to uncover the key players in the early actinorhizal symbiosis, and to further the understanding of how nodule development proceeds after initial contact between host and symbiont.

# Methods

The 100 base-pair single end Illumina reads were analyzed on the University of New Hampshire's training server, RON, as well a personal laptop for R analysis of the produced count matrix. The reads were first checked for quality using `FastQC`, trimmed using `trimmomatic`, assessed for quality once again, and then mapped to the reference genome using `HISAT2`. The resulting alignments were converted to BAMs using `samtools` and then used as evidence for training a gene model using `BRAKER2` and the resulting model was checked for completeness using `BUSCO` to compare across three distinct lineages. Finally, the produced gene model and aligned reads were fed into `featurecounts` to quantify gene expression and analyzed using the `DESeq2` package and R Studio IDE.

## FastQC

The raw reads were fed into the `fastqc.sh` script producing an HTML containing information on sequence quality. The information ranged from GC content, presence of adapter sequences, sequence duplicates, and many more features. The output also included plots to visualize the results.

## Trimmomatic

The raw reads were then trimmed using `trimmomatic` with the sliding window parameter enabled. This resulted in fastq files containing the trimmed reads which were then re-assessed for quality and trimming success.

## HISAT2/samtools

The trimmed reads were aligned to the reference genome using `HISAT2`. First the reference genome was indexed using the `HISAT2` indexing tool. The genome indexes were then used to align the reads. This resulted in a SAM file which was then converted to a BAM file using `samtools`. This conversion allowed for more compact files for downstream usage. 

## BRAKER2

The aligned reads in BAM file format were fed into `BRAKER2` for gene model prediction. `BRAKER2` uses `AUGUSTUS` and `GeneMark-ET` for *ab initio* predictions and training. The resulting gene model was a file in GTF format compatable with downstream tools as `featurecounts`.

## gffread

In order to input the trained gene model into `BUSCO`, a fasta file must be produced. This was accomplished using `gffread` with the gene model in GTF format and reference genome as input. The resulting fasta file included the converted nucleotide to protein sequences.

## BUSCO

The resulting protein fasta from the `gffread` conversion was then tested against three `BUSCO` lineages, fabales_odb10, eudicots_odb10, and embryophyta_odb10. The resulting quality scores were output as txt files containing the percentages of matching orthologs to each of the three lineages.

## featurecounts

Once the quality of the gene model was confirmed, the GTF file and aligned reads in BAM format were used as an input for `featurecounts`, resulting in a txt file output. This txt file contained the counts for each aligned gene and tabular format, and with some manual editting was formatted for DESeq2 analysis.

## DESeq2

The resulting counts matrix from `featurecounts` could now finally be used for differential gene expression analysis using DESeq2. The output was a table of the calculated differential gene expression which could then be plotted using other R libraries such as `ggplot2` and `ggVennDiagram`.

# Results

![pca_plot](https://github.com/user-attachments/assets/6ebe2502-d004-47d1-a1fe-afca9191c3d1)

![all_DEGs_venn](https://github.com/user-attachments/assets/9097d032-825b-4cdd-ace4-16d262c503eb)

![volc_cells_24](https://github.com/user-attachments/assets/4256c84e-2c46-441a-a5ff-32a658f2f4a1)

![volc_cells_48](https://github.com/user-attachments/assets/8dd29cb4-45d9-4748-aee1-b86387278a12)

![volc_NINA_24](https://github.com/user-attachments/assets/8659d1d7-c5aa-4765-9369-d924ce1433e6)

![volc_NINA_48](https://github.com/user-attachments/assets/588612d1-070c-4b43-9458-c8874d3f9046)

# Citations

