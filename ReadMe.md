## Code for reproducing simulations and analyses in "Fast approximation of small p-values in permutation tests by partitioning the permutation space" by Brian Segal, Thomas Braun, and Hui Jiang

The accompanying R package is at [https://github.com/bdsegal/fastPerm](https://github.com/bdsegal/fastPerm).

### Contents:

1. `pval_trend`: Code for Figure 1
2. `simulation_diffMean`: Code for the difference in means simulations
    1. Run code in the `sym` and `nonSym` sub-directories
    2. Run `plots_diff.R`
3. `simulation_ratioMean`: Code for the ratio of means simulations
    1. First, run code in the `sym` and `nonSym` sub-directories
    2. Second, run `plots_ratio.R`
4. `application_cancer`: Code for the analysis with cancer genomic data:
    1. Download data from [TCGA](https://tcga-data.nci.nih.gov/tcga/) and place the unzipped folders in a sub-directory called `data`. The data we downloaded were labeled as `unc.edu_LUAD.IlluminaHiSeq_RNASeqV2.Level_3.1.12.0` and `unc.edu_LUSC.IlluminaHiSeq_RNASeqV2.Level_3.1.8.0`.
    2. Run `cancerAnalysis_parallel.R`
    3. Run `cancerAnalysis_parallel_post.R`
    .4 Run `cancerAnalysis_parallel_replicateTop15.R`
5. `algorithm_schematic`: Code for making a small visual for explaining our resampling algorithm; not included in the paper

Note: The `run` files are batch scripts for submitting jobs via [SLURM](http://slurm.schedmd.com/).