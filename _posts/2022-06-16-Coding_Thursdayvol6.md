---
title: TT-seq analysis - Step 1
date: 2022-06-16 
categories: ["Coding Thursday"]
background: /assets/theme/images/star.png
---

The data I have been working with are the result of TT-seq experiment performed by Shao et al., who set out to unravel how nascent transcription responds to cell state transitions. 
In this method, RNA is pulse labelled with 4sU for few minutes, isolated from cells and fragmented. Fragmented RNA is biotinylated and purified using streptavidin columns. Purified RNA as well as a fraction of not purified RNA (to represent total RNA); LRNA and FRNA, respectively, are subject to RNA sequencing. Data generated in such experiment allows to estimate the rates of RNA synthesis and degradation. I refer you to the original TT-seq paper for further info (Schwalb et al., 2016)
In the next few posts, I will attempt to repeat Shao’s analysis. 

Shao et al. used ERCC spike-in RNAs as a reference for total and labeled RNA samples normalisation. Spike-ins are used to generate a scale-factor to account for differences in library size between multiple samples. 

Therefore, in the first step of analysis I will align paired Illumina sequence reads against a reference genome as well as against spike-in sequences using STAR (Dobin et al., 2013). 
I will mark duplicates (Picard) and index the genome BAM (Li et al., 2009). This is step is not essential, but can provide additional info on sample quality.

```{bash}
#### Set temporary and home (results) directories
export TMPDIR="/path/to/my/working_directory/"
export HOMEDIR="/path/to/my/home_directory/"
cd $TMPDIR
mkdir -p star_final
mkdir -p index #index is stored here (created for spike-ins and Mouse genome)

####Load modules
module load star/2.7.8a
module load samtools/1.3.1
module load picard/2.9.2

#align to the genome/spike-in.  Sort, mark duplicates and index the genome BAM.

for sample in SRR13866843_GSM5137573_FRNA_2i_2d_rep1 SRR13866844_GSM5137574_FRNA_2i_2d_rep2 SRR13866845_GSM5137575_FRNA_2i_7d_rep1 SRR13866853_GSM5137583_LRNA_2i_2d_rep1 SRR13866854_GSM5137584_LRNA_2i_2d_rep2 SRR13866855_GSM5137585_LRNA_2i_7d_rep1
do
	read1=$TMPDIR/fastq/${sample}_Mus_musculus_RNA-Seq_1.fastq.gz
	read2=$TMPDIR/fastq/${sample}_Mus_musculus_RNA-Seq_2.fastq.gz
	STAR --runThreadN 8 --genomeDir $TMPDIR/index \
	--readFilesIn $read1 $read2 --readFilesCommand zcat \
	--quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --outSAMunmapped None \
	--outSAMattrRGline ID:${sample} PU:${sample} SM:${sample} LB:unknown PL:illumina \
	--outFileNamePrefix ${star_final}/${sample}. \
	--outSAMtype BAM Unsorted \
	--outSAMattributes Standard 
	samtools sort --threads 8 -o ${star_final}/${sample}.sorted.bam ${star_final}/${sample}.Aligned.out.bam
	java -jar picard.jar MarkDuplicates INPUT=${star_final}/${sample}.sorted.bam OUTPUT=${star_final}/${sample}.sorted.marked.bam METRICS_FILE=${star_final}/${sample}.sorted.marked.metrics REMOVE_DUPLICATES=false ASSUME_SORTED=true MAX_RECORDS_IN_RAM=2000000 VALIDATION_STRINGENCY=LENIENT TMP_DIR=${tmpdir}/${sample}
	samtools index ${star_final}/${sample}.sorted.marked.bam
	rm ${star_final}/${sample}.sorted.bam
done

```


### References

* Shao, R. et al. Distinct transcription kinetics of pluripotent cell states. Molecular Systems Biology 18:e10407, doi.org/10.15252/msb.202110407 (2022).

* Schwalb, B. et al. TT-seq maps the human transient transcriptome. Science 352(6290):1225-8, doi: 10.1126/science.aad9841 (2016).

* Dobin, A. et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15-21, doi:10.1093/bioinformatics/bts635 (2013).

* Li, H. et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078-2079, doi:10.1093/bioinformatics/btp352 (2009).

* http://broadinstitute.github.io/picard/.
