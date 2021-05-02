# PhyloCap

R package For processing high-throughput sequencing data from targeted sequence capture.


# Installation of prerequisites 

PhyloCap uses several R packages and other outside programs for certain functions. 

1. R packages

From CRAN
- devtools
- ape
- stringr
- data.table
- seqinr
- foreach
- doparallel
- rdrop2

From BioConductor
- rsamtools
- genomicranges
- biostrings

2. Outside programs

- fastp: adaptor trimming and paired-end read merging
- bwa: read mapping
- spades: assembly
- BLAST: matching assembled contigs to targets, other utilities
- mafft: creating alignments
- trimal: trimming alignments
- IqTree: gene tree and concatenation trees
- GATK4: variant calling functions
- SamTools: variant calling and read mapping tools


First, you will want to clone this repository to your computer to obtain the setup files. Or alternatively go to the green "Code" button in top right of this repository and select "download ZIP".

```bash
git clone https://github.com/chutter/PhyloCap.git
```

Second, change your working directory in the terminal to the downloaded repository. 

```bash
cd /Download/Location/PhyloCap
```

The R packages and outside programs can be installed manually or more easily through the anaconda environment file provided (version numbers are provided in environment file if manual installation is desired). To install with the environment file, the easiest and quickest way is to first install the Anaconda package manager. Anaconda can be downloaded and installed for different operating systems from https://anaconda.org. Miniconda is recommended. Once installed, you can create a new environment for PhyloCap by: 

```bash
conda env create -f phylocap.yml -n phylocap
```

OR if a specific location for the environment directory is needed:

```bash
conda env create -f phylocap.yml -p /PLACE/YOUR/DIRECTORY/HERE/phylocap
```

And finally, you may delete the cloned GitHub directory after installing the prerequisites through the conda env file that manually installs the anaconda environment. There are some useful examples (also in the tutorial here), which could be saved.   


# Installation of R package

The main functions of PhyloCap are contained in an R package that has been tested on R version 4.0.2 and use the listed programs above along with custom scripts. To install PhyloCap from GitHub, you can use the R package devtools. Here are step-by-step instructions for installation:

1) Install devtools by typing in your R console: 

```R
install.packages("devtools", dependencies = TRUE)
```

2) Install PhyloCap by typing in your R console: 

```R
devtools::install_github("chutter/PhyloCap", update = "never", dependencies = FALSE)
```

3) If Devtools asks you to install the package dependencies, select "No", because the dependencies are installed above. If devtools asks you to update packages, you may choose to do so. I would recommend not to install packages from source if devtools asks you and to include the update parameter set to "never" in install_github. 

4) Devtools should finish and say the package loaded properly. Load the package in your R script with:

```R
library(PhyloCap)
```

And installation should be done! All the functions for PhyloCap should be ready to go! It is recommended to keep the install line above in your R script as the package is frequently updated for bugs and other features. 

< coming soon a function to test if they can found >


# PhyloCap pipeline tutorials 

[Installation: detailed installation instructions and trouble-shooting ](https://github.com/chutter/PhyloCap/wiki/Installation)

[Tutorial 1: Data preprocessing (cleaning raw reads, decontamination, assembly)](https://github.com/chutter/PhyloCap/wiki/Tutorial-1)

[Tutorial 2: Data alignment (target marker matching, alignment, alignment filtering)](https://github.com/chutter/PhyloCap/wiki/Tutorial-2)

Tutorial 3: Tree estimation (gene trees, concatenation, species trees) COMING SOON

[Troubleshooting: random issues that arise](https://github.com/chutter/PhyloCap/wiki/Troubleshooting)



