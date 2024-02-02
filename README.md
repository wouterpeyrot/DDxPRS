# DDx-PRS

The `DDxPRS` R function provides a tool for distuingishing different disorders based on polygenic prediction. The DDx-PRS (Differential Diagnosis-Polygenic Risk Score) method is described in detail in Peyrot et al. 2024 medRxiv. 


## Compute liability-scale case-control polygenic risk scores (PRS) 
Before applying DDx-PRS, liability-scale case-control PRS need to be computed for the disorders considered. This should be done with a Baysian method to achieve proper callibration of the case-control PRS. We used PRS-CS for this purpose (Ge at al. 2019 Nat Commun), but this can also be done with e.g. SBayesR (Lloyd-Jones et al. 2019 Nat Commun). For more details on computing liability-scale case-control PRS, see https://github.com/euffelmann/bpc. 

**Apply PRS-CS**

* PRS-CS can be downloaded here: https://github.com/getian107/PRScs
* Run PRS-CS-auto with default arguments, i.e. “--a=1 --b=0.5 --phi=None --beta_std=False --n_iter=1000 --thin=5 --n_burnin=500 ---seed=None”.
* In the input for PRS-CS, it is crucial to specify as sample-size the effective N (Neff) of the training GWAS results (i.e. the sum of Neff of the cohorts included in the meta-analyses of the training data; see Grotzinger at al. 2023 Biol Psychiatry for details). 

**Compute case-control PRS using Plink**

Using the SNP-effects thus obtained with PRS-CS, the PRS should be computed in two datasets: (1) in a population reference sample (e.g. 1000G) to assess (a) the liability-scale variance explained in every disorder by its case-control prs (see details below) and (b) to assess the correlations between the case-control prs; and (2) in the test sample in which you aim to apply `DDx-PRS`. The population reference sample, the test sample and the training GWAS results should all be from the same ancestry.
* Compute the PRS with the Plink command “--score header sum center” (it is important to use this exact command to attain proper callibration)
* When your test sample is very small, the allele frequencies used in"--score" can be computed in the reference sample and read with the Plink command “--read-freq” (see for details: https://www.cog-genomics.org/plink/1.9/score)
* The case-control PRS thus obtained are on the standardized observed scale with 50/50 case/control ascertainment, and need to be transformed to the liability-scale (see below).

**Transform observed-scale case-control PRS to the liability-scale**

For each disorder, the lifetime population prevelence (`K`) can be used to transform the observed-scale case-control PRS to the liability-scale in R as follows
```[r]
h2o_to_h2l<-function(K,P,h2o){
        ## Eq 23 of Lee et al. 2011 Am J Hum Genet
        t = -qnorm(K,0,1) ; z = dnorm(t) ;return(h2o*K*K*(1-K)*(1-K)/{z*z*P*(1-P)})
} 
prs_liab <- prs_obs5050*sqrt(h2o_to_h2l(K=K,P=0.5,h2o=1))
``` 
The case-control PRS needs to be transformed to liability-scale in both the population reference sample and in the test sample.


## Getting Started to apply `DDx-PRS`

Copy the `DDxPRS.R` file to your computer. Load DDxPRS and the following libraries.

```[r]
source(DDxPRS.R)
library(mvtnorm)
library(mvnfast)
``` 

If the R packages *mvtnorm* or *mvnfast* have not been installed in R, you can install them with the R command: `install.packages("...")`.



## Running `DDx-PRS`
Before applying DDxPRS, we recommend to study and run the example below. The disorder names should be provided in every input argument listed below. Disorders should be identically ordered for all input arguments.

The input arguments of `DDxPRS()` are:

* **prs_liab:** a dataframe with the liability-scale case-control PRS for each disorder (individuals in rows, prs in columns)

* **prs_r2l:** a vector with liability-scale variance explained in each disorder by its case-control prs, e.g. c(dis1=0.10,dis2=0.08,dis3=0.04). For each disorder, liability-scale variance explained can be approximated by assessing the variance of liability-scale case-control PRS in population referernce sample. See Uffelmann et al. 2024 MedRxiv (PMID: 38260678) for details.

* **snp_h2l:** a vector with liability-scale SNP-heritability for every disorder, e.g. c(dis1=0.24,dis2=0.19,dis3=0.09). Liability-scale SNP-heritabiities can be assessed with LD score regression, available at https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation

* **K:** a vector with lifetime population prevalence for each disorder, e.g. c(dis1=0.01,dis2=0.02,dis3=0.16).

* **crosstrait_cor.prs:** a dataframe of correlations between the case-control PRS. This should be estimated based on the PRS computed in the population reference sample. For illustation, see the example below.

* **crosstrait_rg:** a dataframe with the genetic correlations between the disorders considers. Genetic correlations can be assessed with cross-trait LD score regression, available at https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation

* **clinical.prior:** a vector with the clinical prior probabilities for the diagnostic categories specified in *liab.configuration* (see below), e.g. c(cat1=0.25,cat2=0.25,cat3=0.25,cat4=0.25). Note that the prior probabilities should add up to exactly 1.

* **liab.configuration:** a dataframe linking the configurations of liabilities to the diagnostic categories that you aim to predict. When considering n disorders, there exist 2^n possible configurations of liabilities (above or below the liability threshold for each disorder) (the number of rows of *liab.configuration*). The first n column-names should be the disorder names, and the next colums should have the name *diagnostic.category*, e.g. colnames: c("dis1","dis2","dis3","diagnostic.category"). For illustation, see the example below.





## !!! IGNORE EVERYTHING BELOW !!!





















* **A_name/B_name:** 2 or 3 letter disorder codes that will be used to plot the genetic distance between cases and controls in terms of F<sub>ST,causal</sub>.

* **sumstats_fileA1A0/sumstats_fileB1B0:** names of the gzip-files with case-control GWAS results. Column names should be: SNP, CHR, BP, EA, NEA, FRQ, OR, SE, P, Neff. EA denotes the effective allele, NEA denotes the non effective allele, FRQ denotes the frequency of the EA in controls, OR denotes the OR per EA, SE denotes the standard error of log(OR), Neff labels the effective sample size. When Neff is not known on a SNP-by-SNP basis, this column can be computed by 4/{(1/N_case)+(1/N_control)} estimated per cohort contributing to the meta-analyses and then added together across contributing cohorts. Note that the SE column is optional, although without SE some SNPs with very small p-values < 1e-310 may be discarded (as no z-scores can be estimated for them). Note that the case-control GWAS results of A1A0 and B1B0 will be merged based on SNP names, so make sure to align these adequately (by i.e. naming them as CHR-BP).

* **K_A1A0/K_B1B0:** the most likely lifetime disorder prevalences of disorder A and disorder B in the population, following the disorder definition used in the respective case-control GWAS. 

* **K_A1A0_high/K_B1B0_high:** upper bound of disorder prevalences. This parameter protects against potential false positive association at stress test SNPs (shared causal SNPs with the same allele frequency in cases of both disorders), and acknowledges that it is often hard to now the exact prevalence of the disorder the respective case-control GWAS. When a loose disorder definition has been applied this parameter should be set to a relatively larger value.

* **K_A1A0_low/K_B1B0_low:** lower bound of disorder prevalences. This parameter protects against potential false positive association at stress test SNPs (shared causal SNPs with the same allele frequency in cases of both disorders), and acknowledges that it is often hard to now the exact prevalence of the disorder the respective case-control GWAS. When a strict disorder definition has been applied this parameter should be set to a relatively lower value.

* **h2l_A1A0/h2l_B1B0:** SNP-based heritability on the liability scale, as estimated with e.g. stratified LD Score regression (Bulik-Sullivan et al. 2015A Nat Genet, PMID 25642630; Finucane et al. 2015 Nat Genet, PMID 26414678; Gazal et al. 2017 Nat Genet, PMID 28892061; Gazal et al. 2019 Nat Genet, PMID 31285579)*

* **rg_A1A0_B1B0:** genetic correlation between disorders A and B, as estimated with e.g. cross-trait LD score regression (Bulik-Sullivan et al. 2015B Nature Genetics; PMID: 26414676)

* **intercept_A1A0_B1B0:** intercept from cross-trait LD score regression (Bulik-Sullivan et al. 2015B Nature Genetics; PMID: 26414676) 

* ***m:*** approximation of number of independent effective loci. Our primary recommendation is to specify *m* based on published estimates of genome-wide polygenicity, such as the effective number of independently associated causal SNPs (O' Connor et al. 2019 Am J Hum Genet, PMID 31402091) or the total number of independently associated SNPs (Zhang et al. 2020 Nat Commun, PMID 32620889; Frei et al. 2019 Nat Commun, PMID 31160569; Zhang et al. 2018 Nat Genet, PMID 30104760; Zeng et al. 2018 Nat Genet, PMID 29662166). These values generally range from 1,000 for relatively sparse traits (e.g. autoimmune diseases) to 10,000 for highly polygenic traits (e.g. psychiatric disorders). When estimates of genome-wide polygenicity are not available, our recommendation is to specify *m*=1,000 for traits that are expected to have relatively sparse architectures (e.g. autoimmune diseases), *m*=10,000 for traits that are expected to have highly polygenic architectures (e.g. psychiatric disorders), and *m*=5,000 for traits with no clear expectation. When comparing disorders with different levels of polygenicity, our recommendation is to specify *m* based on the expected average across both disorders.

* **N_A1/N_B1:** total number of cases in the respective input case-control GWAS 

* **N_A0/N_B0:** total number of controls in the respective input case-control GWAS 

* **N_overlap_A0B0:** confirmed number of overlapping controls between the A1A0 and B1B0 GWAS samples. This number can increase power of CC-GWAS as it increases the modelled covariance of error-terms between the case-control GWAS results. When unknown, set to 0 to prevent inflated type I error.

* **subtype_data:** set to `FALSE` when comparing two different disorders (with different definitions of controls, e.g. when comparing psychiatric disorders). Set to `TRUE` when comparing subtypes of a disorder (with same definition of controls, e.g. when comparing subtypes of specific cancer). This setting adjusts the weights of the CC-GWAS<sub>Exact</sub> component to prevent inflated type I error rate at stress test SNP (with causal case-control effects but the same allele frequency among cases).

* **sumstats_fileA1B1:** set to `NA` when applying CC-GWAS based on case-control GWAS results only. Specify the name of the gzip-file with GWAS results from the direct case-case comparison, when applying CC-GWAS+. Format the file as described in *sumstats_fileA1A0/sumstats_fileB1B0* above.

* **N_A1_inA1B1/N_B1_inA1B1:** set to `NA` when applying CC-GWAS based on case-control GWAS results only. When applying CC-GWAS+, specify the number of cases in the direct case-case GWAS here.

* **intercept_A1A0_A1B1/intercept_B1B0_A1B1:** set to `NA` when applying CC-GWAS based on case-control GWAS results only. When applying CC-GWAS+, provide here the intercept from cross-trait LD score regression of A1A0 vs A1B1 respectively B1B0 vs A1B1 (i.e. the intercept from the "Genetic Covariance" section in the output from LD score regression; Bulik-Sullivan et al. 2015B Nature Genetics; PMID: 26414676)

* **save.all:** set to `FALSE` to save only the trimmed results that can directly be used for follow-up analyses (such as e.g. clumping). Set the `TRUE` to save all results, including information of all steps to protect against type I error; note that these results may include SNPs that should be removed to prevent type I error.

## Output files

The `CCGWAS()` function provides three output files. The `outcome_file.log` file provides a logfile of the analyses. This file also reports the CC-GWAS<sub>OLS</sub> weights and the CC-GWAS<sub>Exact</sub> weights. The `outcome_file.pdf` file provides a plot of the genetic distances between cases and controls of both disorders in terms of F<sub>ST,causal</sub>. The `outcome_file.results.gz` file reports results of the case-case association analyses. SNPs with significant case-case association are labelled as 1 in the **CCGWAS_signif** column. The other columns are

* **SNP, CHR, BP, EA, NEA:** as in the input case-control GWAS files.

* **OLS_beta, OLS_se, OLS_pval:** case-case association based on the CC-GWAS<sub>OLS</sub> component. The required level of significance of the CC-GWAS<sub>OLS</sub> component is 5x10<sup>-8</sup> (controlling type I error at null-null SNPs with no impact in either case-control comparison). 

* **Exact_beta, Exact_se, Exact_pval:** case-case association based on the CC-GWAS<sub>Exact</sub> component (based on the most likely lifetime disorder prevalences; see above). The required level of significance of the CC-GWAS<sub>Exact</sub> component is 10<sup>-4</sup> (controlling type I error at stress test SNPs, i.e. shared causal SNPs with the same allele frequency in cases of both disorders). 

* **CCGWAS_signif:** labels SNPs with significant case-case association as `1` (i.e. passing the required levels of significance for the OLS_pval, Exact_pval, Exact_ll_pval, Exact_lh_pval, Exact_hl_pval and Exact_hh_pval, without suggestive evidence for differential tagging of a nearby stress test SNP), and other SNPs as `0`. 

When setting `save.all=TRUE` (see above), CC-GWAS additionally outputs a file `outcome_file.results.ALL.gz`. In addition to the output columns above, this file also contains the following columns:

* **A1A0_beta, A1A0_se, A1A0_pval, B1B0_beta, B1B0_se, B1B0_pval:** case-control effects expressed on the standardized observed scale with a 50/50 case-control ascertainment. 

* **potential_tagging_stresstest:** reports results of the filtering step to exclude potential false positive associations due to differential tagging of a causal stress test SNP (shared causal SNPs with the same allele frequency in cases of both disorders). `NA` when at least one of the CC-GWAS<sub>OLS</sub> component and the CC-GWAS<sub>Exact</sub> component does not reach its required level of significance. `0` when both the CC-GWAS<sub>OLS</sub> component and the CC-GWAS<sub>Exact</sub> component reach their required level of significance, without evidence for differential tagging of a nearby stress test SNP. `1` when both the CC-GWAS<sub>OLS</sub> component and the CC-GWAS<sub>Exact</sub> component reach their required level of significance, but with suggestive evidence for differential tagging of a nearby stress test SNP (these SNPs are excluded from the significant results). 

* **Exact_ll_pval, Exact_lh_pval, Exact_hl_pval, Exact_hh_pval:** reports the p-values of perturbations of the CC-GWAS<sub>Exact</sub> component based on the specified ranges of the lifetime disorder prevalences (see above). *ll* corresponds to weights based on K_A1A0_low and K_B1B0_low (see above); *lh* is based on K_A1A0_low and K_B1B0_high; *hl* is based on K_A1A0_high and K_B1B0_low; *hh* is based on K_A1A0_high and K_B1B0_high. All of these p-values are required to pass the level of significance of 10<sup>-4</sup> to control type I error at stress test SNPs in the context of uncertainty about the population prevalences.

The file `outcome_file.results.ALL.gz` may also contain SNPs with `OLS_pval<5e-8` that should be removed to protect against type I error (which are per default removed from `outcome_file.results.gz`). When wanting to work with `outcome_file.results.ALL.gz`, these SNPs can be removed in R with:

```[r]
library(data.table)
d <- as.data.frame(fread("outcome_file.results.ALL.gz",header=TRUE))
d <- d[ {d$OLS_pval<5e-8 & d$CCGWAS_signif==0}==FALSE ,] 
``` 

## Using `CC-GWAS` results for follow-up analyses

We advise to use the results from the CC-GWAS<sub>OLS</sub> component (OLS_beta, OLS_se, OLS_pval) for clumping and for polygenic risk score analyses. We advise to use the results from the CC-GWAS<sub>Exact</sub> component (Exact_beta, Exact_se, Exact_pval) for genetic correlation analyses.

## Running the example in the *test* folder 

Download the `test.casecontrol.gwas.BIP.10snps.txt.gz`, and `test.casecontrol.gwas.SCZ.10snps.txt.gz` files from the *test* folder and place in your working directory. Run the `CCGWAS()` function with:

```[r]
library(data.table)
library(R.utils)
library(CCGWAS)
CCGWAS( outcome_file = "test.out" , A_name = "SCZ" , B_name = "BIP" , 
        sumstats_fileA1A0 = "./test.casecontrol.gwas.SCZ.10snps.txt.gz" ,
        sumstats_fileB1B0 = "./test.casecontrol.gwas.BIP.10snps.txt.gz" ,
        K_A1A0 = 0.004 , K_A1A0_high = 0.01 , K_A1A0_low = 0.004 ,  
        K_B1B0 = 0.01 , K_B1B0_high = 0.02 , K_B1B0_low = 0.005 , 
        h2l_A1A0 = 0.2 , h2l_B1B0 = 0.20 , rg_A1A0_B1B0 = 0.70 , intercept_A1A0_B1B0 = 0.2425 , m = 1e4 ,  
        N_A1 = 40675 , N_B1 = 20352 , N_A0 = 64643 , N_B0 = 31358 , N_overlap_A0B0 = 24265 )
        
``` 

This provides the results for 10 SNPs from the schizophrenia (SCZ) vs bipolar disorder (BIP) case-case comparison, described in detail in Peyrot & Price. 2021 Nature Genetics.
