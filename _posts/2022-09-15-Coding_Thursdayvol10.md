---
title: TT-seq analysis - Metagene profiles
date: 2022-09-15
categories: ["Coding Thursday"]
background: /assets/theme/images/profile.png
---

In this post, I will use ngs.plot to visualize the results of TT-seq experiment around functional genomic regions, i.e. transcription start site (TSS) and transcription end site (TES).

NGS plot download instructions are available at: https://github.com/shenlab-sinai/ngsplot

I would like to be able to look at sense and antisense profiles. As data I am working with is paired, I will first use SAMtools to get MATE 1 reads. 
The script is called as follows: ```sbatch scripts.sh bam samplename```

where ```bam``` is the input bam file, and ```samplename``` is what will be use to name new files. 

```bash
#!/usr/bin/bash

##modules
module load samtools/1.3.1
module load python/3.5.1

#### Set working, temporary and results directories.
export WORKDIR=/pathtoworkdir/
TMPDIR="${WORKDIR}tmp/"
PROFDIR="${WORKDIR}metaprofiles/"
mkdir -p $TMPDIR
mkdir -p $PROFDIR

#### Sample information: sample name and BAM file location.
SAMPLE=$2
BAM=${WORKDIR}${1}
MATE1="${WORKDIR}${SAMPLE}.mate1.bam"

#### Restict to first mate reads and index 
# If the BAM file contains paired reads, create a new version containing only the first mate reads.<br>
# This step may be skipped if your reads are not paired.<br>
samtools view --threads $THREADS -h -b -f 64 $BAM -o $MATE1
samtools index $MATE1

#### BAM header compliance.
# It might be necessary to re-header the BAM file so that the chromosome names match those in the ngs.plot database, e.g. standard chromosomes preceeded with "chr" for mm10
MATE1REHEADER="${WORKDIR}${SAMPLE}.mate1.reheader.bam"
samtools view --threads $THREADS -H ${MATE1} | sed -e 's/SN:\([0-9XY]*\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - ${MATE1} > ${MATE1REHEADER}
samtools index ${MATE1REHEADER}
```
Mate 1 can now be used to create sense and anti-sense profiles. 

```bash
#!/bin/bash
#call script as: script.sh samplename

export WORKDIR=/pathtoworkdir/
PROFDIR="${WORKDIR}metaprofiles/"
mkdir -p $PROFDIR

THREADS=1
SAMPLE=$1
BAM="${WORKDIR}${SAMPLE}.bam"

REGION="tss"
for STRAND in both same opposite
do
    OUTPUT="${PROFDIR}${SAMPLE}.${REGION}.${STRAND}"
    ngs.plot.r -G mm10 -R $REGION -C ${BAM} -O $OUTPUT -P $THREADS -SS $STRAND -SE 1 -L 5000 -F chipseq -D ensembl
done

REGION="genebody"
for STRAND in both same opposite
do
    OUTPUT="${PROFDIR}${SAMPLE}.${REGION}.${STRAND}"
    ngs.plot.r -G mm10 -R $REGION -C ${BAM} -O $OUTPUT -P $THREADS -SS $STRAND -SE 1 -L 5000 -F chipseq -D ensembl
done

REGION="tes"
for STRAND in both same opposite
do
    OUTPUT="${PROFDIR}${SAMPLE}.${REGION}.${STRAND}"
    ngs.plot.r -G mm10 -R $REGION -C ${BAM} -O $OUTPUT -P $THREADS -SS $STRAND -SE 1 -L 5000 -F chipseq -D ensembl
done

```

As a result, you will get pdf files with profiles for sense, antisense, and both strands on seperate graphs, e.g.:

<p align="center">
![](https://github.com/mmaslon/maslonlab.github.io/blob/dde46f5613c18a8ad9cb85b953ae9ed52482df19/assets/theme/images/genebody.same.avgprof.png "sense strand profile over gene body")
    
as well as zipped files with the txt files with data used for drawing the profiles. 

These files can be used to re-draw the profiles, e.g. drawing sense and antisense coverage on the same plot (see the upcoming blog post). 





