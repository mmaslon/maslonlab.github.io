The data I have been working with are the result of TT-seq experiment performed by Shao et al., who set out to unravel how nascent transcription responds to cell state transitions. 
In this method, RNA is pulse labelled with 4sU for few minutes, isolated from cells and fragmented. Fragmented RNA is biotinylated and purified using streptavidin columns. Purified RNA as well as a fraction of not purified RNA (to represent total RNA); LRNA and FRNA, respectively, are subject to RNA sequencing. Data generated in such experiment allows to estimate the rates of RNA synthesis and degradation. I refer you to the original TT-seq paper for further info (Schwalb et al., 2016)
In the next few posts, I will attempt to repeat Shao’s analysis. 


Shao et al. used ERCC spike-in RNAs as a reference for total and labeled RNA samples normalisation. Spike-ins are used to generate a scale-factor to account for differences in library size between multiple samples. 

Therefore, in the first step of analysis I will align paired Illumina sequence reads against a reference genome as well as against spike-in sequences using STAR. 
I will mark duplicates and index the genome BAM. This is step is not essential, but can provide additional info on sample quality.


### References

* Dobin, A. et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15-21, doi:10.1093/bioinformatics/bts635 (2013).

* Li, H. et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078-2079, doi:10.1093/bioinformatics/btp352 (2009).

* Quinlan AR. BEDTools: The Swiss-Army Tool for Genome Feature Analysis. CurrProtoc Bioinformatics. 2014 Sep 8;47:11.12.1-34.

* http://broadinstitute.github.io/picard/.
