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

Available preprocessing tools:
1. 

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

