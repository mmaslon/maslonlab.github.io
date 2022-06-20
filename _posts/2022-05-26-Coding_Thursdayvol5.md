---
title: Reads alignment
date: 2022-05-26 
categories: ["Coding Thursday"]
background: /assets/theme/images/star.png
---

Today, I will start the analysis of the data and align the pair-end reads to mm10 reference sequence using STAR.

1. Start interactive session.
```srun --pty /bin/bash```

2. Download the reference genome. Navigate to the correct genome on the ensembl page, i.e. use primary assembly rather than top_level, both soft-masked and not masked repeats should work with STAR. For mm10 it will be:

```bash
wget http://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz
```

```bash
wget http://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.chr.gtf.gz
```

3. Aligning reads using STAR requires two steps, i.e. creating a genome index and mapping reads to the genome. 

a) generate genome indices.

The flags to create a genome index using STAR are:

`--runThreadN`: number of threads
`--runMode`: genomeGenerate mode
`--genomeDir`: /path/to/store/index
`--genomeFastaFiles`: /path/to/fasta_file
`--sjdbGTFfile`: /path/to/gtf_file
`--sjdbOverhang`: readlength -1 #

The following script is run to create indices:

```bash
#!/bin/bash
#SBATCH --job-name=index
#SBATCH -p biology_night   #this is cluster specific, here it means the job will run for up to 12h
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mem=30gb
#SBATCH -e error_index.txt    #error output file
#SBATCH -o output_index.txt

# source directories
export TMPDIR=/path/to/temporary-dir
export HOMEDIR=/path/to/home-dir

#load modules
module load star/2.7.8a

# copy required files to temporary directory in scratch area
cp -r $HOMEDIR/genomes/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa $TMPDIR/
cp -r $HOMEDIR/genomes/Mus_musculus.GRCm38.99.chr.gtf $TMPDIR/

cd $TMPDIR

#make subdirectory for the indices
mkdir -p index

#generate index  
#in the case of this data the read lentght is 42, for paired-end data it is a sum of mates' lengths

STAR --runThreadN 6 --runMode genomeGenerate --genomeDir $TMPDIR/index --genomeFastaFiles Mus_musculus.GRCm38.dna_sm.primary_assembly.fa \
--sjdbGTFfile Mus_musculus.GRCm38.99.chr.gtf --sjdbOverhang 83

# Clean up TMPDIR
rm $TMPDIR/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa
rm $TMPDIR/Mus_musculus.GRCm38.99.chr.gtf
```

b) mapping 

Some basic options we will use for aligning reads to the genome using STAR are:

`--runThreadN`: number of threads
`--genomeDir`: /path/to/index_folder
`--readFilesIn`: /path/to/fastq_files 
`--outFileNamePrefix`: prefix for all output files
`--readFilesCommand`: can be used to uncompress the files
`--outSAMtype`: output filetype (SAM default)
`--outSAMunmapped`: what to do with unmapped reads

The script is as follows;

```bash
#!/bin/bash
#SBATCH --job-name=mapping
#SBATCH -p biology_night
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mem=30gb
#SBATCH -e error_mapping.txt
#SBATCH -o output_mapping.txt

export TMPDIR=/home/users/mmaslon/grant_312/scratch/mmaslon/
export HOMEDIR=/home/users/mmaslon/grant_312/project_data/mmaslon

#load modules
module load star/2.7.8a

cd $TMPDIR

#MAKE FEW SUBDIRS UNLESS THEY EXIST
mkdir -p star_final

#align

for sample in SRR13866843_GSM5137573_FRNA_2i_2d_rep1 SRR13866844_GSM5137574_FRNA_2i_2d_rep2 SRR13866845_GSM5137575_FRNA_2i_7d_rep1 SRR13866853_GSM5137583_LRNA_2i_2d_rep1 SRR13866854_GSM5137584_LRNA_2i_2d_rep2 SRR13866855_GSM5137585_LRNA_2i_7d_rep1
do
	read1=$TMPDIR/fastq/${sample}_Mus_musculus_RNA-Seq_1.fastq.gz
	read2=$TMPDIR/fastq/${sample}_Mus_musculus_RNA-Seq_2.fastq.gz
	STAR --runThreadN 6 --genomeDir $TMPDIR/index \
	--readFilesIn $read1 $read2 --readFilesCommand zcat \
	--outFileNamePrefix $sample \
	--outSAMtype BAM SortedByCoordinate \  #output sorted by coordinate
	--outSAMattributes Standard 
done

# Copy output files from TMPDIR back to HOMEDIR, and clean up temporary folder
cp -R $TMPDIR/star_final $HOMEDIR/
cp $HOMEDIR/error_mapping.txt $HOMEDIR/log/
cp $HOMEDIR/output_mapping.txt $HOMEDIR/log/

rm -r $TMPDIR/indices
rm $HOMEDIR/error_mapping.txt
rm $HOMEDIR/output_mapping.txt
```
