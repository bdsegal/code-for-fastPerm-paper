# Code for reproducing simulations and analyses in Segal et al. (2017)

The accompanying R package is at [https://github.com/bdsegal/fastPerm](https://github.com/bdsegal/fastPerm).

## Contents:

1. `pval_trend`: Code for Figure 1
2. `simulation_diffMean_C1`: Code for the difference in means simulations with normal data and equal variances (Web Appendix C.1) and unequal variances (Web Appendix E)
    1. Run code in the `sym` and `nonSym` sub-directories
    2. Run `plots_diff.R` and similarly named scripts
    3. Run code in `other_methods` sub-directory for MCC and saddlepoint approximation
3. `simulation_ratioMean_C2`: Code for the ratio of means simulations with exponential data (Web Appendix C.2)
    1. Run code in the `sym` and `nonSym` sub-directories
    2. Run `plots_ratio.R` and similarly named scripts
4. `simulation_diffMean_C3`: Code for the difference in means simulations with gamma data (Web Appendix C.3)
    1. Run code in the `sym` and `nonSym` sub-directories
    2. Run `plots_gammaDiff_smallN.R` and similarly named script
5. `simulation_ratioMeanGamma_C4`: Code for the ratio of means simulations with gamma data (Web Appendix C.4)
    1. Run code in the `sym` and `nonSym` sub-directories
    2. Run `plots_ratio_smallN.R` and similarly named script
6. `application_cancer`: Code for the analysis with cancer genomic data:
    1. Download data from [TCGA](https://cancergenome.nih.gov/) and place the unzipped folders in a sub-directory called `data`. The data we downloaded were labeled as `unc.edu_LUAD.IlluminaHiSeq_RNASeqV2.Level_3.1.12.0` and `unc.edu_LUSC.IlluminaHiSeq_RNASeqV2.Level_3.1.8.0`.
    2. Run `cancerAnalysis_parallel.R`
    3. Run `cancerAnalysis_parallel_post.R`
    4. Run `cancerAnalysis_parallel_replicateTop15.R`
7. `sample_size': Code for obtaining sufficient sample sizes (Web Appendix F) 
8. `diff_gamma_saddle`: Code for checking saddle point approximation for difference in gamma random variables. Requires the [gammaDist](https://github.com/bdsegal/gammaDist) R package.
9. `algorithm_schematic`: Code for making a small visual for explaining our resampling algorithm; not included in the paper

Note: The `run` files are batch scripts for submitting jobs via [SLURM](http://slurm.schedmd.com/).

## References

Segal, B. D., Braun, T., Elliott, M. R. and Jiang, H. (2017). Fast approximation of small p-values in permutation tests by partitioning the permutations. Biometrics doi:[10.1111/biom.12731](http://onlinelibrary.wiley.com/doi/10.1111/biom.12731/full)