---
title: Accessing and downloading raw sequencing data.
date: 2022-04-21 
categories: ["Coding Thursday"]
background: /assets/theme/images/geo.png
---

1. As confusing as it seems, GEO does not store raw sequence data. Various formats of processed data are uploaded to GEO, including count tables etc. 

For example, upon unpacking GSE168378 repository data downloaded last week:

```tar -xvf GSE168378_RAW.tar```

we find that it contains bigwig files:

```ls```

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

The raw sequence data is deposited in the GEO- linked SRA database.  

2. The relation between GEO and SRA
GSE (GEO study) with many GSM samples : SRP (SRA Study) with many SRA experiments.
The raw data can be found in the SRR accession (SRA run). 
3. Download data Use the tool SRA explorer that allows creating scripts that enable downloading SRR accession from ebb
- stick accession code (SRA identifier for your GSE accession, SRPxxxxxx)
- tick samples to download
- add to collection
- checkout
- number of options: bash script
