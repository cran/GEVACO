## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- results="hide", warning=FALSE, message=FALSE----------------------------
library(GEVACO) # load the library

## -----------------------------------------------------------------------------
data(cov_example) # read in the example covariate data

head(cov_example) # view the format of the covariate file


## -----------------------------------------------------------------------------

data(geno_example) # read in the example genomic file

                   # we included the first 20 SNPs to meet our criteria, each column 
head(geno_example) # is a different SNP


## -----------------------------------------------------------------------------

results <- GxEscreen(dat=cov_example, geno=geno_example, nsim = 1e5, K=7) # run the screen
              # view the first few results: these are the p.values for each SNP in the final
head(results) # genomic file


