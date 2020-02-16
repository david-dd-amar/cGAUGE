# cGAUGE
Causal Graphical Analysis Using GEnetics

## Input

cGAUGE takes as input a set of summary statistics of variants and traits. These are of three types: 
(1) P-value, effect size, and effect size SE matrices from the standard GWAS results (variants as rows and traits are columns)
(2) An object with the p-values of all conditional independence tests for each variant G vs. a trait X. We represent this object using a named list of lists in which element [[tr1]][[tr2]] is a matrix with the conditional independence results (p-values) for trait 1 conditioned on trait 2 (rows are variants).
(3) A binary matrix that specifies the results for all trait-trait CI tests. 

We generate 1-3 above using custom scripts that use R and PLINK. These are currently not available here, but will be added in the next version. Nevertheless, we provide our processed data (i.e., no individual level data) as RData files that can be used for further analysis, please see XX for details.

## Installation and packages

cGAUGE was developed and tested in R 3.4.0
The scripts require the following packages: MendelianRandomization (v 0.4.1 or later), limma (v 3.38.0 or later), and MRPRESSO (v 1.0). Installation time is negligible once these are loaded. However, loading the RData files (processed input data) may take a while (up to 20min) as some are large. 

## Output

