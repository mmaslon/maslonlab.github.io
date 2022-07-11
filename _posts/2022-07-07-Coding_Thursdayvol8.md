---
title: TT-seq analysis - BigWig files
date: 2022-07-07 
categories: ["Coding Thursday"]
background: /assets/theme/images/ttseq.png
---

In this post,  we will use previously calculated scale factors to create scale and strand-specific BigWig files. These files can be used to display data in the Genome Browser as a graph. 
We will use SAMtools split the BAM files into reads mapping to the forward and reverse strands. We will then apply bamCoverage function with -scaleFactor flag (deepTools) to convert strand-specific BAM files to a scaled BigWig files.
The following script can be called using:

```bash
sbatch script bamfile samplename scalefactor
```

script.sh: 
```bash
#!/bin/bash
#SBATCH --job-name=alignment
#SBATCH -p biology_weekend
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mem=30gb
#SBATCH -e error.txt
#SBATCH -o output.txt


#### Set working, temporary and results directories
export WORKDIR=/path/to/working/dir/
export HOMEDIR=/path/to/home/dir/

##modules
module load samtools/1.3.1

TMPDIR="${WORKDIR}tmp/"
BIGWIGDIR="${WORKDIR}bigwigdir/"
mkdir -p $BIGWIGDIR
mkdir -p $TMPDIR

##assign variables used when calling the script to BAM, SAMPLE and SCALEFACTOR variables 
BAM=${HOMEDIR}${1}
SAMPLE=$2
SCALEFACTOR=$3
#sanity check:
echo "for file" $BAM "scalefactor is" $SCALEFACTOR "and samplename is" $SAMPLE

#### Create temporaty BAM files representing reads mapping to forward or reverse strands 
BAMFOR="${TMPDIR}${SAMPLE}.fwd.bam"     
BAMREV="${TMPDIR}${SAMPLE}.rev.bam"     
BAMFOR1="${TMPDIR}${SAMPLE}.fwd1.bam"
BAMFOR2="${TMPDIR}${SAMPLE}.fwd2.bam"
BAMREV1="${TMPDIR}${SAMPLE}.rev1.bam"
BAMREV2="${TMPDIR}${SAMPLE}.rev2.bam"
	
#### Create bigwig files.
BIGWIG="${BIGWIGDIR}${SAMPLE}.bigwig"           # all reads
BIGWIGFOR="${BIGWIGDIR}${SAMPLE}.for.bigwig"    # reads mapping to forward strand
BIGWIGREV="${BIGWIGDIR}${SAMPLE}.rev.bigwig"    # reads mapping to reverse strand

#### Create bigwig file for all reads.
bamCoverage --scaleFactor $SCALEFACTOR -b ${BAM} -o $BIGWIG

#### Create bigwig file for the forward strand.
#Include reads that are 2nd in a pair (128).  Exclude reads that are mapped to the reverse strand (16)

samtools view -b -f 128 -F 16 $BAM > $BAMFOR1

#Exclude reads that are mapped to the reverse strand (16) and first in a pair (64): 64 + 16 = 80
samtools view -b -f 80 $BAM > $BAMFOR2
samtools merge -f $BAMFOR $BAMFOR1 $BAMFOR2
samtools index $BAMFOR
bamCoverage --scaleFactor $SCALEFACTOR -b $BAMFOR -o $BIGWIGFOR

#### Create bigwig file for the reverse strand.
#Include reads that map to the reverse strand (128) and are second in a pair (16): 128 + 16 = 144

samtools view -b -f 144 $BAM > $BAMREV1

#Include reads that are first in a pair (64), but exclude those ones that map to the reverse strand (16)<br>

samtools view -b -f 64 -F 16 $BAM > $BAMREV2
samtools merge -f $BAMREV $BAMREV1 $BAMREV2
samtools index $BAMREV
bamCoverage --scaleFactor $SCALEFACTOR -b $BAMREV -o $BIGWIGREV

#### Remove temporary files.

rm $BAMFOR $BAMFFOR1 $BAMFFOR2 $BAMREV $BAMREV1 $BAMREV2

#### Copy data to home directory
cp -R $BIGWIGDIR $HOMEDIR






