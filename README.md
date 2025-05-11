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

![pca_plot](https://github.com/user-attachments/assets/d6258207-ff28-4d00-840a-0690cb4c438d)

A PCA plot to demonstrate sample similarity. The clustering of sample types demonstrates strong similarity and the high value for PC1 corresponds to a majority of the variance being explained. 

![all_DEGs_venn](https://github.com/user-attachments/assets/9097d032-825b-4cdd-ace4-16d262c503eb)

![volc_cells_24](https://github.com/user-attachments/assets/4256c84e-2c46-441a-a5ff-32a658f2f4a1)

![volc_cells_48](https://github.com/user-attachments/assets/8dd29cb4-45d9-4748-aee1-b86387278a12)

![volc_NINA_24](https://github.com/user-attachments/assets/8659d1d7-c5aa-4765-9369-d924ce1433e6)

![volc_NINA_48](https://github.com/user-attachments/assets/588612d1-070c-4b43-9458-c8874d3f9046)

# Citations

Babraham Bioinformatics—FastQC A Quality Control tool for High Throughput Sequence Data. (n.d.). Retrieved April 29, 2025, from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/![image](https://github.com/user-attachments/assets/a60dcdd6-8aa1-47b6-a9bb-d747c4ac289e)

Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114–2120. https://doi.org/10.1093/bioinformatics/btu170![image](https://github.com/user-attachments/assets/3dcd67c7-34d2-477c-aa96-fcbda255f386)

Lomsadze, A., Burns, P. D., & Borodovsky, M. (2014). Integration of mapped RNA-Seq reads into automatic training of eukaryotic gene finding algorithm. Nucleic Acids Research, 42(15), e119. https://doi.org/10.1093/nar/gku557![image](https://github.com/user-attachments/assets/c9bb9fb9-e058-40be-b845-469fbe77f907)

Manni, M., Berkeley, M. R., Seppey, M., & Zdobnov, E. M. (2021). BUSCO: Assessing Genomic Data Quality and Beyond. Current Protocols, 1(12), e323. https://doi.org/10.1002/cpz1.323![image](https://github.com/user-attachments/assets/89f50514-c7ca-4b33-9b47-8bd1ad6b6858)

Stanke, M., Keller, O., Gunduz, I., Hayes, A., Waack, S., & Morgenstern, B. (2006). AUGUSTUS: Ab initio prediction of alternative transcripts. Nucleic Acids Research, 34(suppl_2), W435–W439. https://doi.org/10.1093/nar/gkl200![image](https://github.com/user-attachments/assets/716f586a-c690-44cc-a43a-e14e10e79623)

Brůna, T., Hoff, K. J., Lomsadze, A., Stanke, M., & Borodovsky, M. (2021). BRAKER2: Automatic eukaryotic genome annotation with GeneMark-EP+ and AUGUSTUS supported by a protein database. NAR Genomics and Bioinformatics, 3(1), lqaa108. https://doi.org/10.1093/nargab/lqaa108![image](https://github.com/user-attachments/assets/3eb14171-b157-418a-b47b-de8ded5eafec)

Kim, D., Paggi, J. M., Park, C., Bennett, C., & Salzberg, S. L. (2019). Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nature Biotechnology, 37(8), 907–915. https://doi.org/10.1038/s41587-019-0201-4![image](https://github.com/user-attachments/assets/dee08bff-ff52-4584-9bf2-14616c9cce5f)

Liao, Y., Smyth, G. K., & Shi, W. (2014). featureCounts: An efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7), 923–930. https://doi.org/10.1093/bioinformatics/btt656![image](https://github.com/user-attachments/assets/ddbefb4d-7fde-4f69-97cc-30bd5fff0979)

Anders, S., & Huber, W. (2010). Differential expression analysis for sequence count data. Genome Biology, 11(10), R106. https://doi.org/10.1186/gb-2010-11-10-r106![image](https://github.com/user-attachments/assets/a7b78df3-31a3-4bac-95e6-0770fd3fe20a)

Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project Data Processing Subgroup. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352![image](https://github.com/user-attachments/assets/afb012ee-368c-49fc-bfe3-f228ae5bb714)

Pertea, G., & Pertea, M. (2020). GFF Utilities: GffRead and GffCompare. F1000Research, 9, ISCB Comm J-304. https://doi.org/10.12688/f1000research.23297.2![Uploading image.png…]()

