---
title: Accessing and downloading raw sequencing data.
date: 2022-04-21 
categories: ["Coding Thursday"]
background: /assets/theme/images/sra.png
---

1. GEO does not store raw data.
As confusing as it seems, GEO does not store raw sequence data. Various formats of processed data are uploaded to GEO, including count tables etc. 

For example, upon unpacking GSE168378 repository data downloaded last week:

```bash
tar -xvf GSE168378_RAW.tar
```

we find that it contains bigwig files:

```bash
ls
filelist.txt                              GSM5137579_FRNA_mTORi_1d_rep1.mm9.rev.bw  GSM5137586_LRNA_SL_rep1.mm9.rev.bw
GSM5137573_FRNA_2i_2d_rep1.mm9.fwd.bw     GSM5137580_FRNA_mTORi_1d_rep2.mm9.fwd.bw  GSM5137587_LRNA_SL_rep2.mm9.fwd.bw
GSM5137573_FRNA_2i_2d_rep1.mm9.rev.bw     GSM5137580_FRNA_mTORi_1d_rep2.mm9.rev.bw  GSM5137587_LRNA_SL_rep2.mm9.rev.bw
GSM5137574_FRNA_2i_2d_rep2.mm9.fwd.bw     GSM5137581_FRNA_mTORi_2d_rep1.mm9.fwd.bw  GSM5137588_LRNA_SL_rep3.mm9.fwd.bw
GSM5137574_FRNA_2i_2d_rep2.mm9.rev.bw     GSM5137581_FRNA_mTORi_2d_rep1.mm9.rev.bw  GSM5137588_LRNA_SL_rep3.mm9.rev.bw
GSM5137575_FRNA_2i_7d_rep1.mm9.fwd.bw     GSM5137582_FRNA_mTORi_2d_rep2.mm9.fwd.bw  GSM5137589_LRNA_mTORi_1d_rep1.mm9.fwd.bw
GSM5137575_FRNA_2i_7d_rep1.mm9.rev.bw     GSM5137582_FRNA_mTORi_2d_rep2.mm9.rev.bw  GSM5137589_LRNA_mTORi_1d_rep1.mm9.rev.bw
GSM5137576_FRNA_SL_rep1.mm9.fwd.bw        GSM5137583_LRNA_2i_2d_rep1.mm9.fwd.bw     GSM5137590_LRNA_mTORi_1d_rep2.mm9.fwd.bw
GSM5137576_FRNA_SL_rep1.mm9.rev.bw        GSM5137583_LRNA_2i_2d_rep1.mm9.rev.bw     GSM5137590_LRNA_mTORi_1d_rep2.mm9.rev.bw
GSM5137577_FRNA_SL_rep2.mm9.fwd.bw        GSM5137584_LRNA_2i_2d_rep2.mm9.fwd.bw     GSM5137591_LRNA_mTORi_2d_rep1.mm9.fwd.bw
GSM5137577_FRNA_SL_rep2.mm9.rev.bw        GSM5137584_LRNA_2i_2d_rep2.mm9.rev.bw     GSM5137591_LRNA_mTORi_2d_rep1.mm9.rev.bw
GSM5137578_FRNA_SL_rep3.mm9.fwd.bw        GSM5137585_LRNA_2i_7d_rep1.mm9.fwd.bw     GSM5137592_LRNA_mTORi_2d_rep2.mm9.fwd.bw
GSM5137578_FRNA_SL_rep3.mm9.rev.bw        GSM5137585_LRNA_2i_7d_rep1.mm9.rev.bw     GSM5137592_LRNA_mTORi_2d_rep2.mm9.rev.bw
GSM5137579_FRNA_mTORi_1d_rep1.mm9.fwd.bw  GSM5137586_LRNA_SL_rep1.mm9.fwd.bw
```

The raw sequence data is deposited in the SRA database, which is linked to GEO.  

2. The relation between GEO and SRA. 

GSExxxxxx (GEO study) with many GSMxxxxxx samples is directly linked to SRPxxxxxx (SRA Study) with many SRA experiments.
The raw data can be found in the SRRxxxxxx accession (SRA run). 

3. Downloading raw data for the given GEO accession.
There are various ways to download sequencing data from SRA. ENA (on EBI) contains a mirror copy of SRA-stored data [ENA](www.ebi.ac.uk/ena). To downlaod raw data associated with GSE168378, I will use a SRA explorer tool, which creates scripts that enable downloading SRR accession from EBI. 
- navigate to [SRA-explorer](https://sra-explorer.info/) website
- on GEO website, find SRA identifier for your GSE accession, SRPxxxxxx (here, for GSE168378 it is SRP309506) 
- stick accession code to SRA explorer search box
- tick samples you wish to download
- press: `add to collection` button, and then `checkout`
- SRA explorer gives number of options to download the data, here we will use bash script:

```bash
#!/usr/bin/env bash
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/043/SRR13866843/SRR13866843_1.fastq.gz -o SRR13866843_GSM5137573_FRNA_2i_2d_rep1_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/043/SRR13866843/SRR13866843_2.fastq.gz -o SRR13866843_GSM5137573_FRNA_2i_2d_rep1_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/044/SRR13866844/SRR13866844_1.fastq.gz -o SRR13866844_GSM5137574_FRNA_2i_2d_rep2_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/044/SRR13866844/SRR13866844_2.fastq.gz -o SRR13866844_GSM5137574_FRNA_2i_2d_rep2_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/048/SRR13866848/SRR13866848_1.fastq.gz -o SRR13866848_GSM5137578_FRNA_SL_rep3_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/048/SRR13866848/SRR13866848_2.fastq.gz -o SRR13866848_GSM5137578_FRNA_SL_rep3_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/046/SRR13866846/SRR13866846_1.fastq.gz -o SRR13866846_GSM5137576_FRNA_SL_rep1_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/046/SRR13866846/SRR13866846_2.fastq.gz -o SRR13866846_GSM5137576_FRNA_SL_rep1_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/047/SRR13866847/SRR13866847_1.fastq.gz -o SRR13866847_GSM5137577_FRNA_SL_rep2_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/047/SRR13866847/SRR13866847_2.fastq.gz -o SRR13866847_GSM5137577_FRNA_SL_rep2_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/045/SRR13866845/SRR13866845_1.fastq.gz -o SRR13866845_GSM5137575_FRNA_2i_7d_rep1_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/045/SRR13866845/SRR13866845_2.fastq.gz -o SRR13866845_GSM5137575_FRNA_2i_7d_rep1_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/049/SRR13866849/SRR13866849_1.fastq.gz -o SRR13866849_GSM5137579_FRNA_mTORi_1d_rep1_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/049/SRR13866849/SRR13866849_2.fastq.gz -o SRR13866849_GSM5137579_FRNA_mTORi_1d_rep1_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/051/SRR13866851/SRR13866851_1.fastq.gz -o SRR13866851_GSM5137581_FRNA_mTORi_2d_rep1_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/051/SRR13866851/SRR13866851_2.fastq.gz -o SRR13866851_GSM5137581_FRNA_mTORi_2d_rep1_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/050/SRR13866850/SRR13866850_1.fastq.gz -o SRR13866850_GSM5137580_FRNA_mTORi_1d_rep2_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/050/SRR13866850/SRR13866850_2.fastq.gz -o SRR13866850_GSM5137580_FRNA_mTORi_1d_rep2_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/053/SRR13866853/SRR13866853_1.fastq.gz -o SRR13866853_GSM5137583_LRNA_2i_2d_rep1_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/053/SRR13866853/SRR13866853_2.fastq.gz -o SRR13866853_GSM5137583_LRNA_2i_2d_rep1_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/054/SRR13866854/SRR13866854_1.fastq.gz -o SRR13866854_GSM5137584_LRNA_2i_2d_rep2_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/054/SRR13866854/SRR13866854_2.fastq.gz -o SRR13866854_GSM5137584_LRNA_2i_2d_rep2_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/052/SRR13866852/SRR13866852_1.fastq.gz -o SRR13866852_GSM5137582_FRNA_mTORi_2d_rep2_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/052/SRR13866852/SRR13866852_2.fastq.gz -o SRR13866852_GSM5137582_FRNA_mTORi_2d_rep2_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/057/SRR13866857/SRR13866857_1.fastq.gz -o SRR13866857_GSM5137587_LRNA_SL_rep2_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/057/SRR13866857/SRR13866857_2.fastq.gz -o SRR13866857_GSM5137587_LRNA_SL_rep2_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/060/SRR13866860/SRR13866860_1.fastq.gz -o SRR13866860_GSM5137590_LRNA_mTORi_1d_rep2_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/060/SRR13866860/SRR13866860_2.fastq.gz -o SRR13866860_GSM5137590_LRNA_mTORi_1d_rep2_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/059/SRR13866859/SRR13866859_1.fastq.gz -o SRR13866859_GSM5137589_LRNA_mTORi_1d_rep1_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/059/SRR13866859/SRR13866859_2.fastq.gz -o SRR13866859_GSM5137589_LRNA_mTORi_1d_rep1_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/055/SRR13866855/SRR13866855_1.fastq.gz -o SRR13866855_GSM5137585_LRNA_2i_7d_rep1_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/055/SRR13866855/SRR13866855_2.fastq.gz -o SRR13866855_GSM5137585_LRNA_2i_7d_rep1_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/058/SRR13866858/SRR13866858_1.fastq.gz -o SRR13866858_GSM5137588_LRNA_SL_rep3_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/058/SRR13866858/SRR13866858_2.fastq.gz -o SRR13866858_GSM5137588_LRNA_SL_rep3_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/056/SRR13866856/SRR13866856_1.fastq.gz -o SRR13866856_GSM5137586_LRNA_SL_rep1_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/056/SRR13866856/SRR13866856_2.fastq.gz -o SRR13866856_GSM5137586_LRNA_SL_rep1_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/061/SRR13866861/SRR13866861_1.fastq.gz -o SRR13866861_GSM5137591_LRNA_mTORi_2d_rep1_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/061/SRR13866861/SRR13866861_2.fastq.gz -o SRR13866861_GSM5137591_LRNA_mTORi_2d_rep1_Mus_musculus_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/062/SRR13866862/SRR13866862_1.fastq.gz -o SRR13866862_GSM5137592_LRNA_mTORi_2d_rep2_Mus_musculus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/062/SRR13866862/SRR13866862_2.fastq.gz -o SRR13866862_GSM5137592_LRNA_mTORi_2d_rep2_Mus_musculus_RNA-Seq_2.fastq.gz
```

DONE!
