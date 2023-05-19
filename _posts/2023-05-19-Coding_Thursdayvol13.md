---
title: Alternative splicing and Snakemake
date: 2023-05-19
categories: ["Coding Thursday"]
background: /assets/theme/images//as.png
---
In this and next few posts, I will log my attempt to replace shell scripting (which I have been using so far for the analysis) with a more 
reliable way to produce scientific worklows, i.e. [Snakemake](https://snakemake.readthedocs.io/en/stable/) [1].

Snakemake is used to create workflows composed of set of steps (rules) that use some input files to create output files. To learn more about Snakemake and start using it, please navigate to: [Snakemake_installation page] (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

A recommended way to install Snakemake is via Conda/Mamba. Again, I used info on [Snakemake_installation page] (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and this [Biostars post](https://www.biostars.org/p/335903/) as a guideline.

I am going to use data from the recently published paper on elongation rate changes in ageing [2]. This paper finds that RNAPII "speeds up" as we age, and life extending interventions, such as dietary restriction or decrease in IGF-1 signaling (knockout of Irs1 gene) can reverse it. Various groups, includin my work have previously showed that changing RNAPII elongation rate affects RNA processing, including splicing. I want to analyze RNAseq data from aged wild type or Irs1 KO mice and look at the Alternative Splicing (AS) changes. 

I downloaded the respective fastq files from the repository. These are single-end reads, three replicates per condition. For all these samples I will perform an alignment to mouse genome, sort the bam file using STAR, I will then analyze AS between these samples using rMATS.

I first created a yml file with list of programs I need to perform my analysis: 

'''name: rna
channels:
 - bioconda
 - conda-forge
 - defaults
dependencies:
 - star=2.7.8a #mapping
 - rMATS=4.1.2  #AS
 - wget
 - samtools=1.6.0 '''

1. Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33.
2. Debès, C., Papadakis, A., Grönke, S. et al. Ageing-associated changes in transcriptional elongation influence longevity. Nature 616, 814–821 (2023).
