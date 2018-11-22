# crisprtools

crisprtools is a package to allow for the analysis of CRISPR Screens.
It includes 3 algorithms to analyze CRISPR screens developed by other individuals (See Authors):

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

	1. Download Anaconda python 3.x from: https://www.anaconda.com/download/#macos
	2. Once installed, you can create separate envrionments from the main installation. 
	3. In terminal type in (creates a separate python 3.x from main): conda create -n py3x python=3.x
	4. In terminal type in (creates a python 2.7 environment): conda create -n py27 python=2.7
```

### Installing crisprtools

```
You first will need to install devtools if you do not already have it
	1. install.packages("devtools")
	2. load library(devtools)

To install crisprtools package
	3. install_github("singjc/crisprtools")
```

## Current Functions in crisprtools
```
1. Current_CRISPR_algos
	1. drugZ.py
	2. BAGEL.py
	3. rjacks.R
2. rdrugZ.R (re-wrote drugZ into Rscript)
3. Raw_Counts_PreProcessing (calls bhatt.coeff.R)
4. algo_comparison_scatter_plot_single 
```

## Examples

### Running Current_CRISPR_algos()

You can run all three algorithms in one call given you provide all the input arguments.
In this case, algo = c('drugZ', 'BAGEL', 'JACKS')

#### Running python drugZ algorithm

There is also an implemntation of the drugZ algorithm in R, callable function rdrugZ(), requires the same arguments as python version.

```
# If you are using RStudio and have library(rstudioapi)
readcount_input_file = selectFile() # <- will return the path and file name of selected file
Or
readcount_input_file = '/usr/bob/CRISPR_Screen/data/Screen1/read-counts/screen1_counts.txt'

# Define where you want to store and name with file extension drugZ results
drugZ_output_file = '/usr/bob/CRISPR_Screen/data/Screen1/drugZ_Results/screen1_drugZ.txt'
drugZ_foldchange_file = '/usr/bob/CRISPR_Screen/data/Screen1/drugZ_Results/screen1_drugZ-foldchange.txt'

# Define string array of what columns contain control and what columns contain treatment
drugZ_control = 'NT_repA,NT_repB'
drugZ_treatmnet = 'DRUG_repA,DRUG_repB'

# Note: arguments must be named, due to arguments being positional, you need to specify argument value correspondence
Current_CRISPR_algos(algo_call='drugZ', readcount_input_file=readcount_input_file, drugZ_output_file=drugZ_output_file, drugZ_foldchange_file=drugZ_foldchange_file, drugZ_control,drugZ_treatmnet=drugZ_control,drugZ_treatmnet)
```

#### Running python BAGEL algorithm
```
# If you are using RStudio and have library(rstudioapi)
readcount_input_file = selectFile() # <- will return the path and file name of selected file
Or
readcount_input_file = '/usr/bob/CRISPR_Screen/data/Screen1/read-counts/screen1_counts.txt'

# Define full path and file name (WITHOUT .Extension) for bagel fodlchange
bagel_foldchange_out_file = '/usr/bob/CRISPR_Screen/data/Screen1/BAGEL_Results/screen1_BAGEL-foldchange' # <- NOTE: no .txt at the end of screen1_BAGEL-foldchange

# Specify the column that contains the T0 sample, since BAGEL is implemented in python, which starts its indexing at 0 and not 1, you need to take this into account. In addition the first column, mostlikely the sgRNA column, is treated as an index rownames, so indexing starts at the 'second' column GENE (iindex = 0)
T0_control_col = 13 # <- this would in reality be 15 if you count the columns starting from 1 and counting the very first column

# Define where you want to store and name with file extension BAGEL results
bf_out_file = paste(dirname(readcount_input_file),'/TEST_bf.txt',sep='')
# String array of which column to include in BAGEL analysis
bf_columns = '1,2'

Current_CRISPR_algos(algo_call='BAGEL', readcount_input_file, bagel_foldchange_out_file=bagel_foldchange_out_file, T0_control_col=T0_control_col, bf_out_file=bf_out_file, bf_columns=bf_columns)
```

#### Running R JACKS algorithm
```
# If you are using RStudio and have library(rstudioapi)
readcount_input_file = selectFile() # <- will return the path and file name of selected file
Or
readcount_input_file = '/usr/bob/CRISPR_Screen/data/Screen1/read-counts/screen1_counts.txt'

# This should be the path and filename to the replicate mapping file
repmapfile <- '/usr/bob/CRISPR_Screen/data/Screen1/read-counts/repmap.txt'

# an example repmap file should look something as follows in a .txt tab-deliminated file.
#
# Replicate 	Sample
# NT_repA		CTRL
# NT_repB		CTRL
# Drug1_repA	Drug1
# Drug1_repB	Drug1
# Drug2_repA	Drug2
# Drug2_repB	Drug2


# Identify which comlumn in the repmap file denotes the Replicate column information
jacks_replicate_col="Replicate"

# Identify which comlumn in the repmap file denotes the Sample column information
jacks_sample_col="Sample"

# Indentify which columns in the read-count file contains the sgRNA column and the Gene column
jacks_grna_col="sgRNA"
jacks_gene_col="Gene"

# Identify which is your reference control sample in the repmap file
jacks_reference_sample="CTRL"

# These are Default values
jacks_count_prior=32
jacks_normalization='median'
jacks_window=800

# If you are using RStudio and have library(rstudioapi)
jacks_saveDir = selectDirectory()
Or
jacks_saveDor = '/usr/bob/CRISPR_Screen/data/Screen1/JACKS_Results/'

# These are default values
jacks_target_genes = NULL
jacks_n_iter = 5
jacks_reference_library = NA

Current_CRISPR_algos(algo_call='JACKS',
                     readcount_input_file=readcount_input_file,
                     repmapfile=repmapfile,
                     jacks_replicate_col=jacks_replicate_col,
                     jacks_sample_col=jacks_sample_col,
                     jacks_grna_col=jacks_grna_col,
                     jacks_gene_col=jacks_gene_col,
                     jacks_count_prior=jacks_count_prior,
                     jacks_normalization=jacks_normalization,
                     jacks_reference_sample=jacks_reference_sample,
                     jacks_saveDir=jacks_saveDir,
                     jacks_target_genes=jacks_target_genes,
                     jacks_n_iter=jacks_n_iter,
                     jacks_reference_library=jacks_reference_library)
```


## Authors

* **Traver Hart** - *BAGEL* - [SourceCode](https://github.com/hart-lab/drugz) and [Publication](https://www.biorxiv.org/content/early/2017/12/12/232736)
* **Gang Wang** - *drugZ* - [SourceCode](https://github.com/hart-lab/bagel) and [Publication](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1015-8)
* **Fellicity Allen** - *JACKS* - [SourceCode](https://github.com/felicityallen/JACKS) and [Publication](https://www.biorxiv.org/content/early/2018/04/12/285114)

## Acknowledgments

* **Thomas Guillerme** - *bhatt.coef.R* - [SourceCode](https://github.com/TGuillerme/Total_Evidence_Method-Missing_data/blob/master/Functions/TreeComparisons/bhatt.coef.R)