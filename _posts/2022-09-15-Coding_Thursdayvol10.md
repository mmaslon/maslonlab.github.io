---
title: TT-seq analysis - Metagene profiles
date: 2022-09-15
categories: ["Coding Thursday"]
background: /assets/theme/images/igv.png
---

In this post, I will use ngs.plot to visualize the results of TT-seq experiment around functional genomic regions, i.e. transcription start site (TSS) and transcription end site (TES).

NGS plot download instructions are available at: https://github.com/shenlab-sinai/ngsplot

I would like to be able to look at sense and antisense profiles. As data I am working with is paired, I will first use SAMtools to get MATE 1 reads. 
The script is called as follows: ```sbatch scripts.sh bam.file samplename```

where 

```
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


#### Threads - set to take advantage of multi-threading and speed things up.
THREADS=1

#### Restict to first mate reads.
# If the BAM file contains paired reads, create a new version containing only the first mate reads.<br>
# This step may be skipped if your reads are not paired.<br>
samtools view --threads $THREADS -h -b -f 64 $BAM -o $MATE1
samtools index $MATE1


#### BAM header compliance.
# It might be necessary to re-header the BAM file so that the chromosome names match those in the ngs.plot database, e.g. standard chromosomes preceeded with "chr" for hg38.  This is most easily achieved using SAMtools "reheader" function.<br>
# For human hg38 Ensembl alignments the following steps should do the job.<br>
# This step may be skipped if your header is already compliant.<br>
MATE1REHEADER="${WORKDIR}${SAMPLE}.mate1.reheader.bam"
samtools view --threads $THREADS -H ${MATE1} | sed -e 's/SN:\([0-9XY]*\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - ${MATE1} > ${MATE1REHEADER}
samtools index ${MATE1REHEADER}


BAMFOR="${MATE1REHEADER}.fwd.bam" 
BAMREV="${MATE1REHEADER}.rev.bam"     # BAM file representing reads mapping to reverse strand

# Forward strand.

```
