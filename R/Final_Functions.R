#' Joint test for genetic and gene-environment interaction effects for a single SNP
#'
#' Function to test for the joint genetic and gene-environment
#' interaction effects for a single variant using a LRT model
#'
#' @references Crainiceanu, C. M., & Ruppert, D. (2004). Likelihood ratio tests in linear mixed models with one variance component.
#' Journal of the Royal Statistical Society Series B-Statistical Methodology, 66, 165-185. doi: 10.1111/j.1467-9868.2004.00438.x
#'
#' @param dat a data frame with covariate information. Column 1 should be phenotype,
#' column 2 should be the environmental factor of interest, columns 3 and later should be
#' additional covariates
#' @param snp_geno a vector containing genotypic information of SNP of interest to be tested
#' @param nsim the number of replicates in obtaining the p-value (standard 1e5)
#' @param K the number of knots used to control the flexibility in modeling GxE interaction
#'
#' @importFrom stats lm quantile
#' @importFrom nlme lme pdIdent
#' @importFrom RLRsim exactLRT
#'
#' @return empirical p-value obtained as the proportion of T0 that are greater than
#' the observed test statistic T
#' @export

GxEtest <- function(dat, snp_geno, nsim = 1e5, K){
    X1 = dat

    #Rename columns to make them consistent: V1, V2, ..., Vn
    colnames(X1) = paste0("V",1:ncol(X1))

    n <- nrow(dat)
    group <- rep(1,n)

    X2 <- data.frame(snp_geno, X1$V2*snp_geno)


    X <- cbind(X1,X2)
    #Rename columns to make them consistent: V1, V2, ..., Vn
    colnames(X) = paste0("V",1:ncol(X))

    # moving default knots function to here
    # calculates knots in terms of vector of environmental factor
    default_knots_fx <- function(x,num.knots)
    {
        if (missing(num.knots)) num.knots <- max(5,min(floor(length(unique(x))/4),35))
        return(quantile(unique(x),seq(0,1,length=(num.knots+2))[-c(1,(num.knots+2))]))
    }

    knots <- default_knots_fx(X1$V2,K)
    Z2 <- outer(X1$V2, knots,"-")
    Z2 <- Z2*(Z2>0)

    Z2 <- Z2*snp_geno #random effects for genetic interaction effects

    Znames <- paste("Z",1:K,sep="")

    dimnames(Z2)[[2]] = Znames

    X = data.frame(cbind(X, group, Z2))

    #Use <paste(paste0("V",3:ncol(X)), collapse = '+')> to represent the list of covariates
    fixed_formula = paste("V1~",paste(names(X)[startsWith(names(X), "V")][-1], collapse="+"))
    random_formula = paste("~-1+",paste(names(X)[startsWith(names(X), "Z")], collapse="+"))

    fit_full = eval(parse(text = paste0("nlme::lme(",fixed_formula, ", random=list(group=nlme::pdIdent(", random_formula, ")), data=X)")))


    fit_reduced <- lm(V1~.,data=X1)

    p_val <- suppressMessages(RLRsim::exactLRT(fit_full,fit_reduced,seed = 2019,nsim = nsim)$p.value) #suppress messages from the function to make it clean
    return(p_val)
}

#' Gene-Environment Interaction: Genome-wide Screen
#'
#' Function to test for the joint genetic and gene-environment
#' interaction effects for a set of variant using a LRT model
#'
#'
#' @param dat a data frame with covariate information. Col 1 should be phenotype,
#' col 2 should be environmental factor, col 3 and later should be additional covariates
#' @param geno a genotype matrix with 0-1-2 coding
#' @param nsim the number of replicates in obtaining the p-value (standard 1e5)
#' @param K the number of knots used to control the flexibility in modeling GxE interaction
#' @return a vector containing the p-value from the LRT associated with each SNP
#' @examples
#' GxEscreen(cov_example, geno_example, nsim=1e5, K=7)
#' @export


GxEscreen <- function(dat, geno, nsim=1e5, K=7){
    p_vals <- c()

    for(i in 1:ncol(geno))
        p_vals <- c(p_vals, GxEtest(dat=dat, snp_geno=geno[, i], nsim = nsim, K=K))

    return(p_vals)
}


