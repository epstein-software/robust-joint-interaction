
## We wrote this R plugin to be consistent with PLINK guidelines
## For more general details on constructing R plugins, please see the PLINK website
##
## Within the R plugin function, we define
##
## PHENO:    vector of phenotypes contained in the PLINK *.ped file
## GENO:     matrix of SNP genotypes
## CLUSTER:  vector of cluster membership codes (this optional argument is ignored by this particular R plugin)
## COVAR:    matrix of covariates listed in the PLINK *.cov file
##
##==================================================================================================
## 
## Assumptions: 
##
## -For joint testing of SNP and SNP-environment interaction, we assume the first covariate listed in the *.cov file 
## is the environmental factor to be tested in the interaction term
##
## -We adjust only for the main effects of the remaining covariates listed in the *.cov file	
##
##==================================================================================================

## The R library geepack must be installed on local machine before using the R plugin:
library(geepack)

Rplink <- function(PHENO, GENO, CLUSTER, COVAR)
{
    f1 <- function(s){    # s denotes the genotype vector at the tested SNP

        ## Define the interaction variable between the test SNP and the environmental exposure
        ## As mentioned in Assumptions, we assume interaction with the first covariate listed in the
        ## covariate file:
        ge_int <- s*COVAR[,1] 

        ## The R geepack library requires each observation to have a subject id variable associated with it:
        subid <- seq(1, length(PHENO)) 

        ## Important: Need to remove subjects with missing data before analysis
        ## The R geepack library will crash if there are any missing data (NAs) in the dataset to be analyzed:
        olddat <- cbind(PHENO, s, COVAR, ge_int, subid) 

        ## Find rows containing missing data:
        badid <- unique(na.action(na.omit(olddat)))

        ## Remove rows containing missing data from dataset:
        if(length(badid) >= 1){
            PHENO  <- PHENO[-c(badid)]
            s      <- s[-c(badid)]
            COVAR  <- COVAR[-c(badid),]
            ge_int <- ge_int[-c(badid)]
            subid  <- subid[-c(badid)]
        }

        COVAR_mat <- as.matrix(COVAR) # formatting change needed for geepack

        ## Fit the regression model
        ## For continuous outcomes, specify linear regression model using family=gaussian()
        ## For binary outcomes, specify logistic regression model using family=binomial()
        ##
        ## While not specified, the model does include an intercept:
        gee.fit <- geese(PHENO ~ s + ge_int + COVAR_mat, id=subid, family=gaussian(), corstr='independence')

        ## Construct the joint test of SNP and SNP-environment interaction based on the fitted model
        ## Must specify the A matrix described in supplementary information of the Almli et al. (2014; JAMA-Psychiatry paper)
        ##
        ## The A matrix has 2 rows and P columns, where P is the number of regression parameters (excluding intercept) 
        ## in the 'gee.fit' model above:
        temp_matrix <- cbind(s, ge_int, COVAR)
        ncolA <- ncol(temp_matrix) + 1         # Need to add the 1 here to allow for the intercept
        Amat <- matrix(0, nrow=2, ncol=ncolA)

        ## In the 'gee.fit' model above:
        ##   the first parameter was the intercept,
        ##   the second parameter was the main effect of the SNP genotype,
        ##   the third parameter was the interaction effect between SNP genotype and environment
        Amat[1,2] <- 1   # Signifies the first parameter we want to test is the effect size of the SNP genotype
        Amat[2,3] <- 1   # Signifies the second parameter we want to test is the effect size of the SNP-environment interaction

        Vbeta_robust <- gee.fit$vbeta     # Estimate of the robust variance-covariance matrix of regression coefficients
        Vbeta_model <- gee.fit$vbeta.naiv # Estimate of the model-based variance-covariance matrix of regression coefficients

        beta_est <- gee.fit$beta  # Estimate of the regression coefficients from the regression model


        ## Using this information, we can construct the joint test as shown in the supplemental material 
        ## of the JAMA-Psychiatry paper:
        Ab <- Amat %*% beta_est

        AVA_robust <- Amat %*% Vbeta_robust %*% t(Amat)
        AVA_model  <- Amat %*% Vbeta_model %*% t(Amat)

        # the function solve takes the inverse of a matrix:
        test_robust <- t(Ab) %*% solve(AVA_robust) %*% Ab
        test_model  <- t(Ab) %*% solve(AVA_model) %*% Ab

        p_robust <- pchisq(test_robust, 2, lower.tail=F)  # Calculate the p-value of the robust test
        p_model  <- pchisq(test_model, 2, lower.tail=F)   # Calculate the p-value of the model-based test

        r <- c(p_model,p_robust)

        ## Have the function return the p-values of the two tests:
        c(length(r), r)

    }  # f1()

    apply(GENO,2,f1)

} # Rplink()
