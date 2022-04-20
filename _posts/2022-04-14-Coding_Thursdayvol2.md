---
title: Downloading data from GEO to HPC.
date: 2022-04-14 
categories: ["Coding Thursday"]
background: /assets/theme/images/geo.png
---

Most, if not all (published) NGS data is deposited to GEO repository and available to download and process. Today the objective is to get the data of interest from GEO repository onto our HPC cluster accounts. 

1. Log-in to HPC as usual using `ssh` command in the Terminal.
2. Data should be downloaded in the interactive mode using: 
```srun --pty /bin/bash``` 
Once in the interactive mode, create an appropriate directory for your data.
3. Locate the data of interest in GEO repository. In this example, I want to download the following data: `GSE168378`. 
To find the link to transfer the data, we should access the FTP site - you can locate it under “Tools” header on GEO homepage. In the `series/` folder of FTP site find and enter the appropriate directory (here `GSE168nnn/`) and then `GSE168378/`
The  data is available in the `suppl/` directory. Right click/double tap (if you are a Mac user) to get the link for this directory.
4. Use the `wget` command with the following options: `--recursive` (as we are copying a directory) and `--no-parent` and `nd` (to avoid copying any parent folders) 

```wget --recursive —-no-parent -nd ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE168nnn/GSE168378/suppl/```

DONE!
