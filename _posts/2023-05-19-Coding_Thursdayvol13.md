---
title: Alternative splicing and Snakemake
date: 2023-05-19
categories: ["Coding Thursday"]
background: /assets/theme/images//as.png
---
In this and next few posts, I will log my attempt to replace shell scripting (which I have been using so far for the analysis) with a more 
reliable way to produce scientific worklows, i.e. [Snakemake](https://snakemake.readthedocs.io/en/stable/) [1].

Snakemake is used to create workflows composed of set of steps (rules) that use some input files to create output files. To learn more about Snakemake and start using it, please navigate to: [Snakemake_installation page](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

A recommended way to install Snakemake is via Conda/Mamba. Again, I used info on [Snakemake_installation page](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and this [Biostars post](https://www.biostars.org/p/335903/) as a guideline.

I am going to use data from the recently published paper on elongation rate changes in ageing [2]. This paper finds that RNAPII "speeds up" as we age, and life extending interventions, such as dietary restriction or decrease in IGF-1 signaling (knockout of Irs1 gene) can reverse it. Various groups, includin my work have previously showed that changing RNAPII elongation rate affects RNA processing, including splicing. I want to analyze RNAseq data from aged wild type or Irs1 KO mice and look at the Alternative Splicing (AS) changes. 

I downloaded the respective fastq files from the repository. These are single-end reads, three replicates per condition. For all these samples I will perform an alignment to mouse genome, sort the bam file using STAR, I will then analyze AS between these samples using rMATS.

1. I first created a yml file with list of programs I need to perform my analysis (env.yml). At the top of the file I specified name of the environment I am going to create.  

```
name: rna
channels:
 - bioconda
 - conda-forge
 - defaults
dependencies:
 - star=2.7.8a #mapping
 - rMATS=4.1.2  #AS
 - wget
 - samtools=1.6.0
 ```
 2. I then created conda environment using env.yml file

```bash
conda env create -f env.yml
```
3. As a sanity check, I verified that it has indeed been created: 

```
conda env list
```

4. Next, I activate my environment calls `rna` and installed Snakemake within that environemt. 

```
conda activate rna
```
```
conda install snakemake
```
5. I prepared few files first to run my workflow `config.json` (path to reads, genome files), `wt.txt` and `ko.txt` (list of bam files for splicing analysis)

```bash
configfile:
    "config.json"

SAMPLES, = glob_wildcards(config['data']+"/{id}.fastq.gz")
EVENTS = ["A3SS", "A5SS", "MXE", "RI", "SE"]
JCS = ["JC", "JCEC"]
BAM_FILES = expand("rnaseq/star/{sample}Aligned.sortedByCoord.out.bam", sample = SAMPLES)
AS_FILES = expand("rnaseq/rmats2/{event}.MATS.{jc}.txt", event = EVENTS, jc = JCS)

rule all:
    input:
        AS_FILES,
        BAM_FILES

rule index_genome:
    input:
        fa = 'rnaseq/genome/mm10.fa',
        gtf = 'rnaseq/genome/mm.gtf'
    params:
        outdir = directory('rnaseq/genome/index')
    output:
        "mockfile.txt",
    shell:
        'mkdir -p {params.outdir} && '
        'touch mockfile.txt && '
        'STAR --runThreadN 4 --runMode genomeGenerate --genomeDir {params.outdir} --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --sjdbOverhang 99'  #reads 100bp

rule align:
    input:
        read = config['data']+"/{sample}_1.fastq.gz", 
        genome = directory('rnaseq/genome/index/')
    params:
        prefix = 'rnaseq/star/{sample}_'
    output:
        'rnaseq/star/{sample}_Aligned.sortedByCoord.out.bam',
        'rnaseq/star/{sample}_Log.final.out'
    message:
        'mapping {wildcards.sample} to genome'
    conda:
        "env.yaml"
    shell:
        "mkdir -p {params.outdir}; "
        "cd {params.outdir}; "
        "STAR --runThreadN 4 --genomeDir {input.genome} --readFilesIn {input.read} --readFilesCommand gunzip -c --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard --alignEndsType EndToEnd"

rule analysesplicing:
    input:
        b1 = 'rnaseq/star/wt.txt',
        b2 = 'rnaseq/star/ko.txt',
        gtf = 'rnaseq/genome/mm.gtf'
    output:
        AS_FILES
    params:
        outdir = directory('rnaseq/rmats2'),
        tmp = directory('rnaseq/rmats2/tmp')
    conda:
       "/Users/magdamaslon/Documents/Magdypliki/analysis/env.yml"
    shell:
       "mkdir -p {params.outdir}; "
       "mkdir -p {params.tmp}; "
       "rmats.py --b1 {input.b1} --b2 {input.b2} --gtf {input.gtf} -t single --variable-read-length --readLength 100 --libType fr-firststrand --od {params.outdir} --tmp {params.tmp}"
```
6. check that the rules are OK by doing a dry run:

```
snakemake --np
```
7. run your workflow

```
snakemake --cores 2 --use-conda
```

Please check what we learnt about AS changes in next post. 

References:
1. Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33.
2. Debès, C., Papadakis, A., Grönke, S. et al. Ageing-associated changes in transcriptional elongation influence longevity. Nature 616, 814–821 (2023).
