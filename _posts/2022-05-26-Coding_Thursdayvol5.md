Reads alignment.

Pair-end reads are aligned to mm9/mm10 genome references (GENCODE) using STAR v …

1. Start interactive session.
srun --pty /bin/bash

2. If you do not have assembly data yet, you should download it. Navigate to the correct genome, i.e. use primary assembly rather than top_level, both soft-masked and not masked repeats should work 

wget http://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz

wget http://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.chr.gtf.gz


2. Aligning reads using STAR will require two steps, i.e. creating a genome index and mapping reads to the genome. 
