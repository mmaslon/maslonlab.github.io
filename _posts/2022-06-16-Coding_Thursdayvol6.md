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


### References

* Shao, R. et al. Distinct transcription kinetics of pluripotent cell states. Molecular Systems Biology 18:e10407, doi.org/10.15252/msb.202110407 (2022).

* Schwalb, B. et al. TT-seq maps the human transient transcriptome. Science 352(6290):1225-8, doi: 10.1126/science.aad9841 (2016).

* Dobin, A. et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15-21, doi:10.1093/bioinformatics/bts635 (2013).

* Li, H. et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078-2079, doi:10.1093/bioinformatics/btp352 (2009).

* http://broadinstitute.github.io/picard/.
