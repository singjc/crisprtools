# crisprtools

crisprtools is a package to allow for the analysis of CRISPR Screens.
It includes 3 algorithms to analyze CRISPR screens developed by other individuals (See [Authors](## Authors)):

1. drugZ - was developed in the Traver Hart lab, by Gang Wang 

2. BAGEL - was developed in the Traver Hart lab, by Traver Hart 

3. JACKS - was developed in the Leopold Parts lab, by Fellicity Allen 

It also includes an algorithm to filter out sgRNA's based on their consistency between replicates and conditions. It uses the Bhattacharyya Coefficient to calculate the overlap between two distributions of replicate samples. The bhatt.coeff.R script was written by Thomas Guillerme.

Will add more useful information soon...

## Getting Started

### Prerequisites

First, you need to have R and RStudio (latest version suggested) installed.

```
1. Install R from: https://www.r-project.org
2. Install RStudio from: https://www.rstudio.com/products/rstudio/download/
```


To run drugZ (python 3) and BAGEL (python 2) you will need to have python installed
```
To install python, you can install python using Anaconda.

...1. Download Anaconda python 3.x from: https://www.anaconda.com/download/#macos
...2. Once installed, you can create separate envrionments from the main installation. 
...3. In terminal type in (creates a separate python 3.x from main): conda create -n py3x python=3.x
...4. In terminal type in (creates a python 2.7 environment): conda create -n py27 python=2.7
```

### Installing crisprtools

```
You first will need to install devtools if you do not already have it
...1. install.packages("devtools")
...2. load library(devtools)

To install crisprtools package
...3. install_github("singjc/crisprtools")
```

## Current Functions in crisprtools
```
1. Current_CRISPR_algos
...1. drugZ.py
...2. BAGEL.py
...3. rjacks.R
2. rdrugZ.R (re-wrote drugZ into Rscript)
3. Raw_Counts_PreProcessing (calls bhatt.coeff.R)
4. algo_comparison_scatter_plot_single 
```

## Running the tests

Will add soon....

## Authors

* **Traver Hart** - *BAGEL* - [SourceCode](https://github.com/hart-lab/drugz) and [Publication](https://www.biorxiv.org/content/early/2017/12/12/232736)
* **Gang Wang** - *drugZ* - [SourceCode](https://github.com/hart-lab/bagel) and [Publication](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1015-8)
* **Fellicity Allen** - *JACKS* - [SourceCode](https://github.com/felicityallen/JACKS) and [Publication](https://www.biorxiv.org/content/early/2018/04/12/285114)

## Acknowledgments

* **Thomas Guillerme** - *bhatt.coef.R* - [SourceCode](https://github.com/TGuillerme/Total_Evidence_Method-Missing_data/blob/master/Functions/TreeComparisons/bhatt.coef.R)