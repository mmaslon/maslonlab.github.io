---
title: Downloading fastq files.
date: 2022-05-19 
categories: ["Coding Thursday"]
background: /assets/theme/images/sra.png
---
Let's pick it up where we left it last time. SRA explorer provides as with basic bash script that we can use to download fastqfiles. I saved in my bin directory under the name sra_download.sh. 


```bash
#!/bin/bash -l
#SBATCH -N 1
#SBATCH --mem 5000
#SBATCH --time=20:00:00

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/053/SRR13866853/SRR13866853_1.fastq.gz -o SRR13866853_GSM5137583_LRNA_2i_2d_rep1_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/053/SRR13866853/SRR13866853_2.fastq.gz -o SRR13866853_GSM5137583_LRNA_2i_2d_rep1_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/055/SRR13866855/SRR13866855_1.fastq.gz -o SRR13866855_GSM5137585_LRNA_2i_7d_rep1_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/055/SRR13866855/SRR13866855_2.fastq.gz -o SRR13866855_GSM5137585_LRNA_2i_7d_rep1_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/043/SRR13866843/SRR13866843_1.fastq.gz -o SRR13866843_GSM5137573_FRNA_2i_2d_rep1_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/043/SRR13866843/SRR13866843_2.fastq.gz -o SRR13866843_GSM5137573_FRNA_2i_2d_rep1_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/045/SRR13866845/SRR13866845_1.fastq.gz -o SRR13866845_GSM5137575_FRNA_2i_7d_rep1_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/045/SRR13866845/SRR13866845_2.fastq.gz -o SRR13866845_GSM5137575_FRNA_2i_7d_rep1_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/054/SRR13866854/SRR13866854_1.fastq.gz -o SRR13866854_GSM5137584_LRNA_2i_2d_rep2_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/054/SRR13866854/SRR13866854_2.fastq.gz -o SRR13866854_GSM5137584_LRNA_2i_2d_rep2_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/044/SRR13866844/SRR13866844_1.fastq.gz -o SRR13866844_GSM5137574_FRNA_2i_2d_rep2_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/044/SRR13866844/SRR13866844_2.fastq.gz -o SRR13866844_GSM5137574_FRNA_2i_2d_rep2_Mus_musculus_RNA-Seq_2.fastq.gz
```
The few first lines in the scripts tell the script to use bash (`#!/bin/bash`), `#SBATCH` configures slurm.

The script has been executed using:

```bash
sbatch /path/to/script/script.sh
```

