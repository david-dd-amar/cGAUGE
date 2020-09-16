# cGAUGE: Causal Graphical Analysis Using GEnetics

cGAUGE is a set of tools that utilize conditional independence (CI) tests for improving causal inference among traits using genetic variables. 

There are three types of analyses that can be performed:
1. <em>ImpIV</em>: Filter out improper genetic instruments before Mendelian Randomization (MR) analysis between two traits. 
1. <em>UniqueIV</em>: Obtain a set of proper genetic instruments for MR. However, this is limited to trait pairs that have limited (conditional) dependence.
1. <em>ExSep</em>: Search for evidence for the exitence of a causal pathway between two traits.

The analyses above require three preprocessing steps that can be done using external tools like PLINK or R. For notation let **T** be the set of traits and **G** be the set of genetic variants (after LD clumping or pruning).
1. <em>GWAS</em>: Obtain genome-wide association results for each genetic variable in **G** vs. every trait in **T**.
1. <em>Trait skeleton</em>: Compute a "skeleton" among the traits witout using the genetic data. Skeletons are undirected graphs in which edges represent traits that cannot be rendered independent by conditioning on other traits.
1. <em>CI tests</em>: For each significant variant obtained in the GWAS step at a significance level p<sub>1</sub> for a trait tr<sub>1</sub>, test if it is no longer significant with p>p<sub>2</sub> when conditioning on another trait tr<sub>2</sub>. Keep the results for every triplet (g in **G**, tr<sub>1</sub>, tr<sub>2</sub>).

We provide several tools for computing the preprocessing steps above, but any custom scripts or tools can be used. Moreover, summary statistics from these steps can be used as well as input for cGAUGE. If such results are available you can skip the next section. 

## Preprocessing tools

#### GWAS and CI tests

For GWAS and CI tests on large-scale genetic data we use PLINK. Assume we have the following files: (1) pheno.phe: a phenotype file with a row per subject and three columns: family id (FID), individual id (IID), and the phenotype scores, which can be binary or continuous, (2) covars.phe: covariates file with a row per subject, the first two columns are FID and IID, and additional columns with covariates that we need to adjust for: sex, age, Array, and genetic principal components, and (3) a bed file with the individual-level genetic data. We then use the following command for logistic regression:
```
plink2 --bfile bed --logistic hide-covar firth-fallback \ 
  --pheno pheno.phe --covar covars.phe \
  --covar-name sex,age,Array,PC1,PC2,PC3,PC4,PC5 \
  --chr 1-22 --maf 0.01 
  --out ${output_path_name}
```
For linear regression we use `--linear hide-covar`. When running CI tests of analysis 3 above, make sure to add the name of tr<sub>2</sub>, that is use: `--covar-name sex,age,Array,PC1,PC2,PC3,PC4,PC5,tr_2_name`

For analysis of without plink (e.g., for smaller datasets) we provide a few useful functions in [R/auxil_functions.R](R/auxil_functions.R) including: `run_lm(x,y,z,df)` which computes the effect size, standard error, and p-value for x~y|z - i.e., the linear effect of y on x when conditioned on z when all variables are available in the data frame df. A more complex wrapper is called `run_ci_test_one_is_numeric(x,y,z,df)` that assumes that either x or y are numeric (or both) and internally decides how to use correlation analysis or linear regression to compute the p-value for x,y|z. Finally, `run_ci_logistic_test(x,y,z,df)` can be used to get the logistic p-value when x is a binary variable.

#### Skeletons

Skeletons can be computed using the [pcalg](https://cran.r-project.org/web/packages/pcalg/index.html) R package. Specifically, the `skeleton` function can be used to get the skeleton of all variables in a data frame. If you use this function make sure to use `m.max=2` to avoid an exponential running time by testing all possible sets (that are conditioned upon). Moreover, for speedups you can use the functions above as input e.g., `indepTest=run_lm` to use linear regression instead of discrete tests. Note that you can also increase `numCores` to run the tests in parallel. However, this is limited to the number of cores in the machine. To better utilize resources in an HPC, we provide a useful R script [hpc_stanford/analyze_trait_pair_for_skeleton.R](hpc_stanford/analyze_trait_pair_for_skeleton.R) that receives as input: (1) an RData file with a data frame that has all sample-level data, (2) the name of tr<sub>1</sub>, (3) the name of tr<sub>2</sub>, and (4) the maximal set cardinallity (equivalent to `m.max` above). Thus, this script can be used to get the CI result for a single trait pair, and thus can be run in parallel for many or all pairs. 

## cGAUGE: Input

For individual level data over a set of traits and a set of genetic variants you can follow the preprocessing steps above to obtain all summary statistics that cGAUGE requires. We provide below these results for 96 traits and their genetics results from the UK-Biobank:

1. An object with the p-values of all conditional independence tests for each variant g in **G** vs. a trait x in **T**. We represent this object using a named list of lists in which element [[tr1]][[tr2]] is a matrix with the conditional independence results (p-values) for trait 1 conditioned on trait 2 (rows are variants). The results for the UK-Biobank data are available [here](https://drive.google.com/file/d/1XNZSYlDnepnPdLgG5qBrtTHrlo2Yq7IG/view?usp=sharing). Here is an example code for using the provided results:
```
# useful objects in this:
# code2gwas_res - a list with the GWAS results, limited to p < 1e-04
# code2clumped_list - a list with the LD clump results per trait
# cl_unified_list - the union of the LD clump sets (can be further pruned)
# snp_P_matrix: a matrix of p-values (row per variant, column per trait)
# sum_stat_matrix: a matrix of effect sizes (row per variant, column per trait)
# sum_stat_se_matrix: a matrix of effect size standard error (row per variant, column per trait)
# OPTIONAL
# code2pruned_list - a list with the LD prune results per trait
# pr_unified_list - the union of the LD prune sets (can be further pruned)
load("single_gwas_res.RData")
> snp_P_matrix[1:3,1:3]
                statins Alanine_aminotransferase   Albumin
rs761193    3.63141e-05                 0.895412 0.5172800
rs116112655 2.16605e-05                 0.241373 0.7728420
rs115045185 4.19488e-05                 0.455153 0.0189549
> sum_stat_matrix[1:3,1:3]
              statins Alanine_aminotransferase     Albumin
rs761193    0.9188275               -0.0132903  0.01290840
rs116112655 1.0116227                0.1315650 -0.00641066
rs115045185 1.0493625               -0.0928334 -0.05761930
```

2. A matrix that contains the skeleton analysis results [here](https://drive.google.com/file/d/1CGav4eGQLi-G1zCdqyrSXbGL8b_aseGM/view?usp=sharing). This link provides a square **|T|** X **|T|** matrix with the maximal p-value for each pair (tr<sub>1</sub>, tr<sub>2</sub>). That is, the maximal p-value obtained for the association of the pair (tr<sub>1</sub>, tr<sub>2</sub>) when trying to condition on another trait from **T**. Here is an example R code for loading this matrix and obtaining a skeleton graph:
```
> load("Gs_skeleton.RData")
> pmax_network[1:3,1:3]
                              statins Alanine_aminotransferase      Albumin
statins                            NA             1.074641e-33 2.132759e-05
Alanine_aminotransferase 1.074641e-33                       NA 6.701961e-01
Albumin                  2.132759e-05             6.701961e-01           NA
# Use p=1e-07 to get an undirected graph
> skel = pmax_network < 1e-07
# Number of edges in the binary network
> table(skel[lower.tri(skel)])
FALSE  TRUE 
 4082   670 
```
3. Numeric matrices with a row for each genetic variant and a column for each trait. For running both the cGAUGE filters and MR analysis you will need three matrices: (1) P-values, (2) effect sizes, and (3) effect size standard error. These are available for the UK-Biobank data in a single RData file [here](https://drive.google.com/file/d/1XNZSYlDnepnPdLgG5qBrtTHrlo2Yq7IG/view?usp=sharing).

```
load("genetic_CI_tests_results.RData")
tr1 = "Glucose"
tr2 = "Albumin"
tr1_given_tr2 = trait_pair_pvals[[tr1]][[tr2]]
# take the last column - contains the tests with adjustment for covariats
tr1_given_tr2_p = tr1_given_tr2[,ncol(tr1_given_tr2)]
# for comparison, take the p-values of tr1 without conditioning on tr2
tr1_ps = snp_P_matrix[,tr1]
shared_snps = intersect(names(tr1_ps),names(tr1_given_tr2_p))
plot(tr1_ps[shared_snps],tr1_given_tr2_p[shared_snps],pch=20,
     xlab = "Glucose p-values",ylab = "Glucose p-values, cond on Albumin")
```
The resulting figure shows that there are only mild effects on the GWAS p-values of tr1 (Glucose) when conditioning on tr2 (Albumin):


<img src="figures/gluc_cond_albumin.png" width="250">

## cGAUGE: Analysis and output




## Additional comments and contact info

cGAUGE was developed and tested in R 3.4.0

The scripts require the following packages: MendelianRandomization (v 0.4.1 or later), limma (v 3.38.0 or later), and MRPRESSO (v 1.0).

For suggestions please contact us at davidama AT stanford dot edu


