# cGAUGE: Causal Graphical Analysis Using GEnetics

cGAUGE is a set of tools that utilize conditional independence (CI) tests for improving causal inference among traits using genetic variables. 

There are three types of analyses that can be performed:
1. <em>ImpIV</em>: Filter out improper genetic instruments before Mendelian Randomization (MR) analysis. 
1. <em>UniqueIV</em>: Obtain a set of proper genetic instruments for MR. However, this is limited to trait pairs that have limited (conditional) dependence.
1. <em>ExSep</em>: Search for evidence for the exitence of a causal pathway between two traits.

The analyses above require three preprocessing steps that can be done using external tools like PLINK or R:
1. <em>GWAS</em>: Obtain genome-wide clumped or pruned association results for each trait.
1. <em>Trait skeleton</em>: Compute a "skeleton" among the traits, witout using the genetic data. Skeleton graphs represents traits that cannot be rendered independent by conditioning on other traits.
1. <em>CI tests</em>: For each significant variant obtained in the GWAS step at a significance level p<sub>1</sub> for a trait tr<sub>1</sub>, test if it is no longer significant with p>p<sub>2</sub> when conditioning on another trait tr<sub>2</sub>.

We provide several tools for computing the preprocessing steps above, but any custom scripts or tools can be used. Moreover, summary statistics from these steps can be used as well as input for cGAUGE. If such results are available you can skip the next section. 

## Preprocessing tools

1. For GWAS and CI tests on large-scale genetic data we use PLINK. Assume we have the following files: (1) pheno.phe: a phenotype file with a row per subject and three columns: family id (FID), individual id (IID), and the phenotype scores, which can be binary or continuous, (2) covars.phe: covariates file with a row per subject, the first two columns are FID and IID, and additional columns with covariates that we need to adjust for: sex, Array, and genetic principal components, and (3) a bed file with the individual-level genetic data. We then use the following command for logistic regression:
```
plink2 --bfile bed --logistic hide-covar firth-fallback \ 
  --pheno pheno.phe --covar covars.phe \
  --covar-name sex,age,Array,PC1,PC2,PC3,PC4,PC5 \
  --chr 1-22 --maf 0.01 
  --out ${output_path_name}
```
For linear regression we use `--linear hide-covar`.
1. For analysis of small datasets without plink we provide a few useful functions in [R/auxil_functions.R](R/auxil_functions.R) including: `run_lm(x,y,z,df)` which computes the effect size, standard error, and p-value for x~y|z - i.e., the linear effect of y on x when conditioned on z when all variables are available in the data frame df. A more complex wrapper is called `run_ci_test_one_is_numeric(x,y,z,df)` that assumes that either x or y are numeric (or both) and internally decides how to use correlation analysis or linear regression to compute the p-value for x,y|z. Finally, `run_ci_logistic_test(x,y,z,df)` can be used to get the logistic p-value when x is a binary variable.

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

