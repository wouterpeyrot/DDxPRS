# DDx-PRS

The `DDxPRS` R function provides a tool for distuingishing different disorders based on polygenic prediction. The DDx-PRS (Differential Diagnosis-Polygenic Risk Score) method is described in detail in [Peyrot et al. 2024 medRxiv](https://www.medrxiv.org/content/10.1101/2024.02.02.24302228v1). Before applying `DDx-PRS`, the liability-scale case-control PRS need to be computed for each disorder in a population references sample (e.g. 1000G) and in your test sample.


## Compute liability-scale case-control polygenic risk scores (PRS) 
Before applying DDx-PRS, liability-scale case-control PRS need to be computed for each disorder considered. This should be done with a Baysian method to achieve proper callibration of the case-control PRS. We used PRS-CS for this purpose (Ge at al. 2019 Nat Commun), but this can also be done with e.g. SBayesR (Lloyd-Jones et al. 2019 Nat Commun). For more details on computing liability-scale case-control PRS, see https://github.com/euffelmann/bpc. 

**Apply PRS-CS**

* PRS-CS can be downloaded here: https://github.com/getian107/PRScs
* Run PRS-CS-auto with default arguments, i.e. “--a=1 --b=0.5 --phi=None --beta_std=False --n_iter=1000 --thin=5 --n_burnin=500 ---seed=None”.
* In the input for PRS-CS, it is crucial to specify as sample-size the effective N (Neff) of the training GWAS results (i.e. the sum of Neff of the cohorts included in the meta-analyses of the training data; see Grotzinger at al. 2023 Biol Psychiatry for details). 

**Compute case-control PRS using Plink**

Using the SNP-effects thus obtained with PRS-CS, the PRS should be computed in two datasets: (1) in a population reference sample (e.g. 1000G) to assess (a) the liability-scale variance explained in every disorder by its case-control prs (see details below) and (b) to assess the correlations between the case-control prs; and (2) in the test sample in which you aim to apply `DDx-PRS`. The population reference sample, the test sample and the training GWAS results should all be from the same ancestry.
* Compute the PRS with the Plink command “--score header sum center” (it is important to use this exact command to attain proper callibration)
* When your test sample is very small, the allele frequencies used in "--score" can be computed in the reference sample and read with the Plink command “--read-freq” (see for details: https://www.cog-genomics.org/plink/1.9/score)
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


## Getting started to apply `DDx-PRS`

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

* **liab.configuration:** a dataframe linking the configurations of liabilities to the diagnostic categories that you aim to predict. When considering n disorders, there exist 2^n possible configurations of liabilities (above or below the liability threshold for each disorder) (the number of rows of *liab.configuration*). The first n column-names should be the disorder names, and the next colum should be named *diagnostic.category*, e.g. colnames: c("dis1","dis2","dis3","diagnostic.category"). For illustation, see the example below.

## Output of `DDx-PRS`
The output of `DDx-PRS` consists of a list with three elements:

* **post_prob:** a dataframe with the following columns:
  * The posterior probabilities for each diagnostic category specified in *liab.configuration*
  * The posterior probabilities for each configuration of liabilities in *liab.configuration* 

* **liab.configuration:** (this output can safely be ignored) an updated version of *liab.configuration*, with the following columns added:
  * analyt_pop_proportion: the expected prevalence of each configuration of liabilities in the population, approximated based on the lifetime population prevalences for each disorder and the genetic correlations between disorders.
  * test_priorprob: the prior probability of each configuration of liabilities, combining the information of the *clinical.prior* for eacht diagnostic category and *analyt_pop_proportion*

* **mvnorm_list:** (this output can safely be ignored) a list with 2^n + 1 items when considering n disorders. The first item of the list contains an approximation of the means (*mu*) and variances/covariances (*sigma*) of the liability-scale case-control PRS and liabilities in the full population. The following items of the list contain the means and covariances for each configuration of liabilities specified in *liab.configuration*.

## Running the example in the *test* folder 

Download the `test.population.reference.sample.txt` and `test.test.sample.txt` files from the *test* folder and place in your working directory. Run the `DDxPRS()` function with:

```[r]
rm(list=ls())
source("DDxPRS.R")
library(data.table)
library(mvtnorm)
library(mvnfast) 

h2o_to_h2l<-function(K,P,h2o){
	## Eq 23 of Lee et al 2011 Am J Hum Genet
	t = -qnorm(K,0,1) ; z = dnorm(t) ; return(h2o*K*K*(1-K)*(1-K)/{z*z*P*(1-P)})
} 

disorder_names <- c("dis1","dis2","dis3")
K <- c(0.01,0.02,0.16)                ## population prevelance ; ordered in line with disorder_names
snp_h2l <- c(0.24,0.19,0.09)          ## SNP-h2l, assessed with e.g. sLDSC ; ordered in line with disorder_names
crosstrait_rg <- as.data.frame(rbind( ## rg, assessed with e.g. cross-trait LDSC; ordered in line with disorder_names
   c( 1.00 , 0.7  , 0.35 )
  ,c( 0.7  , 1.00 , 0.45 )
  ,c( 0.35 , 0.45 , 1.00 )
)) 
names(K) <- names(snp_h2l) <- names(snp_h2l) <- rownames(crosstrait_rg) <- colnames(crosstrait_rg) <- disorder_names

## load PRS data & transpose from standardized 50/50-scale to liability-scale
ref.sample <- as.data.frame(fread("test/test.population.reference.sample.txt"))
test.sample <- as.data.frame(fread("test/test.test.sample.txt"))
for(dis in disorder_names){
	ref.sample [,paste("prs_",dis,"_liab",sep="")] <- ref.sample [,paste("prs_",dis,"_5050",sep="")]*sqrt(h2o_to_h2l(K=K[dis],P=0.5,h2o=1))
	test.sample[,paste("prs_",dis,"_liab",sep="")] <- test.sample[,paste("prs_",dis,"_5050",sep="")]*sqrt(h2o_to_h2l(K=K[dis],P=0.5,h2o=1))
}

## compute liability-variance explained by prs based on reference sample
prs_r2l <- diag(var(ref.sample[,c("prs_dis1_liab","prs_dis2_liab","prs_dis3_liab")])) ## ordered in line with disorder_names
names(prs_r2l) <- disorder_names

## estimate correlation between prs based on reference sample
crosstrait_cor.prs <- as.data.frame(cor(ref.sample[,c("prs_dis1_liab","prs_dis2_liab","prs_dis3_liab")])) ## ordered in line with disorder_names
rownames(crosstrait_cor.prs) <- colnames(crosstrait_cor.prs) <- disorder_names

## prepare prs_liab with prs on liability scale
prs_liab <- test.sample[,c("prs_dis1_liab","prs_dis2_liab","prs_dis3_liab")] ## ordered in line with disorder_names

## specify liab.configuration & clinical.prior
liab.configuration <- as.data.frame(rbind( ## disorder-names ordered in line with disorder_names
   c( dis1=1 , dis2=1 , dis3=1 , diagnostic.category = "cat1")
  ,c( dis1=1 , dis2=1 , dis3=0 , diagnostic.category = "cat1")
  ,c( dis1=1 , dis2=0 , dis3=1 , diagnostic.category = "cat1")
  ,c( dis1=1 , dis2=0 , dis3=0 , diagnostic.category = "cat1")
  ,c( dis1=0 , dis2=1 , dis3=1 , diagnostic.category = "cat2")
  ,c( dis1=0 , dis2=1 , dis3=0 , diagnostic.category = "cat2")
  ,c( dis1=0 , dis2=0 , dis3=1 , diagnostic.category = "cat3")
  ,c( dis1=0 , dis2=0 , dis3=0 , diagnostic.category = "cat4")
))
clinical.prior <- c(cat1=0.25,cat2=0.25,cat3=0.25,cat4=0.25) ## lables in diagnostic.category should correspond to names in liab.configuration

## run DDxPRS
output <- DDxPRS( prs_liab=prs_liab 
                  ,prs_r2l=prs_r2l 
                  ,snp_h2l=snp_h2l 
                  ,K=K  
                  ,crosstrait_cor.prs=crosstrait_cor.prs 
                  ,crosstrait_rg=crosstrait_rg
                  ,clinical.prior 
                  ,liab.configuration=liab.configuration )

##############
## Examples to assess performance of output
##############

library(locfit) 
library(pROC)
ICI <- function(Y,P){
  ## Austin PC. The Integrated Calibration Index (ICI) and related metrics for quantifying the calibration of logistic regression models. Statistics in Medicine 2019
  na.index <- is.na(Y) | is.na(P) ; Y <- Y[na.index==FALSE] ; P <- P[na.index==FALSE]
  loess.calibrate <- locfit(Y~P) 
  P.calibrate <- predict (loess.calibrate, newdata = P)
  ICI <- mean (abs(P.calibrate-P),na.rm=TRUE)   
  return(ICI)
}

Y <- as.numeric(test.sample$DDx=="cat1") ; P <- output$post_prob$prob_cat1
ICI(Y=Y,P=P)
roc(response=Y,predictor=P,quiet=TRUE)$auc
        
``` 

Note that this example consists of simulated data for the analyses described in Peyrot et al. 2024 medRxiv, i.e. predicting the diagnostic categories schziophrenia (cat1), bipolar disorders (cat2), major depressive disorder (cat3) and controls (cat4) based on the case-control PRS of schizophrenia (dis1), bipolar disorder (dis2) and major depressive disorders (dis3).
