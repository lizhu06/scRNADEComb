# scRNADEcomb

scRNADEComb combines several commonly used methods for DE analysis of scRNA-seq data, including:

* MAST
* MAST_det: correct for cell detecting rate in MAST
* Monocle
* ROTS
* DEseq
* SCDE
* ZINB, ZINB_F, ZINB_R: our newly developped methods for dataset including multiple subjects (manuscript in preparation)


### Installing

Download the code using git-clone

```
git clone https://github.com/lizhu06/scRNADEcomb.git
```

### Function
The main function to call is scRNADEComb.R. Parameters to specify:

* count: gene count matrix, gene by sample
* x: vector, length of sample size, case control index
* sample_id (optional, only needed for ZINB_F and ZINB_R), vector, length of sample size, indicate sample id.
* method: vector of methods to implement. Choose among "MAST", "MAST_det", "ROTS", "Monocle",
		"ZINB", "ZINB_F", "ZINB_R", "SCDE", "DESeq2".
* time_name: The name of an R object saved to store the computing time of each methods
* log_s_hat (optional): offset, if not specified, use log transformed total UMI.
* dropout_rate_below (optional): dropout rate cutoff, genes with dropout rate above this cutoff will be filtered out.
* sd_above (optional): standard deviation (SD) cutoff, genes with SD below this cutoff will be filtered out.
* ncores=1: number of cores for parallel.
* seed=123: random seed
* saveDir: the directory to save the final results
* TMBfolder (only needed for ZINB, ZINB_R, ZINB_F): the location of SCMM.TMB.v3.R, SCMM.cpp and SCMMRan.cpp

## Running the tests

An example script (example.R) is provided. Change the directory accordingly.

