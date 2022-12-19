---
title: TT-seq analysis -Spliced reads 
date: 2022-12-22
categories: ["Coding Thursday"]
background: /assets/theme/images/metaprofiles2.png
---

In this post, I will continue with reproducing the analysis from Shao et al and compare spliced reads across various 4sU labelled samples. The assumption is that in the short metabolic labelling experiments, the majority of reads should come from unspliced pre-mRNAs.

In order to extract spliced reads from the data, we will extract information in the 6th column of sam file (the CIGAR string) - see code below. The CIGAR string provides information on how your reads align to the reference sequence. Gaps in the alignments (indicative of splicing) are indicated as “N” in the CIGAR string. I refer you to the following blog for a really clear information on the [CIGAR](https://jef.works/blog/2017/03/28/CIGAR-strings-for-dummies/) 

```{bash}
samtools view -h alignment.bam | awk '$6 ~ /N/ || $1 ~ /^@/' | samtools view -bS - > spliced.bam
```

The first part of the statement uses samtools view to convert the bam to sam and include the header (-h). 
The awk command is then used to filter reads that contain N in the 6th column (the CIGAR string) or start with @ (therefore are header).
The second view statement converts the output back to bam format.
