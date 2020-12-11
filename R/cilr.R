library(fitdistrplus)
library(mixtools)
library(compositions)
library(rsample)
library(purrr)
library(glue)
library(rlist)
library(Rfast)

#' @title Function to calculate cilr scores under different outputs 
#' @param X Count data in matrix form as an n x p matrix  
#' @param A Set identification in matrix form as an p x m matrix 
#' @param resample Boolean indicating whether resampling methods are used to estimate the true score 
#' @param output The preferred output which are 'sig', 'pvalue', 'cdf', 'zscore'
#' @param distr Character indicating which distribution is used in resampling option 
#' @param nperm Number of permutations for resampling option 
#' @param init The initialization vector for estimating distribution
#' @param adj  Boolean, indicating whether correlation adjustment will be performed 
#' @param thresh Threshold if sig is returned to return significance 
cilr <- function(X, A, resample, output = c("cdf","zscore", "pval", "sig"), 
                 distr = c("mnorm", "norm"), nperm=5, init=NULL, adj=TRUE, thresh=0.05, 
                 preprocess=TRUE, pcount=1, transform = "prop", ...){
    output <- match.arg(output)
    print(glue("Output is {opt}", opt = output))
    distr <- match.arg(distr)
	# check type 
    if (is.matrix(X) == F){
        message("Coercing X to matrix")
        X <- as.matrix(X)
    }
    if (is.vector(A) == F){
        message("Coercing A to a numeric vector")
        A <- as.matrix(A)
    }
    # if (nperm > 5 && pryr_object_size(X) >= 100){
    #     cost <- pryr::object_size(X) * nperm 
    #     warning("Beware! Resampling will cost {c} memory!", c = cost)
    # }
	# processing
	if(preprocess == T){
	    X <- process(X, pcount = pcount, transform = transform)
	}

	# second, generate cilr scores  
	p <- ncol(X) # number of features 
	n <- nrow(X) # number of samples 
	R <- matrix(0, ncol = ncol(A), nrow = nrow(X))
	# generate bootstrapped samples  
	if (resample == T){
		if (missing(distr)){
			warning("No distribution chosen, defaulting to the normal distribution")
		} 
        if (missing(init) | is.null(init)){
            message("Default initialization")
        }
		# generate permuted values and bootstrap values that are unpermuted
		perm <- map(1:nperm, ~{
			X[,sample(1:p, replace = FALSE)]
		})
		perm <- do.call(rbind, perm)
		if (adj == TRUE){
		 	unperm <- rsample::bootstraps(X, times = nperm) 
		 	train_list <- purrr::map(unperm$splits, training)
		 	unperm <- do.call(rbind, train_list)
		 	rm(train_list)
		 	gc()
		}
	}

	# loop through each set to generate scores and estimate null distribution 
	if (resample == T){
	  message("Resampling to get null distribution")
	} 
	if (adj == T){
	  message("Adjusting for correlation...")
	} else {
	  message("Not adjusting for correlation...")
	}
	for (i in seq(ncol(A))){
		score <- get_score(X, A[,i])
		if (resample == T){
			sc_perm <- get_score(perm, A[,i])
			perm_dist <- estimate_distr(sc_perm, distr = distr, init = init, ...)
			rm(sc_perm)
			gc()
			if (adj == TRUE){ # if correlation adjustment is true
				sc_unperm <- get_score(unperm, A[,i] %>% as.vector())
                #saveRDS(sc_unperm, "cache/testing.rds")
                #print(
                #    estimate_distr(as.vector(sc_unperm), distr = "norm", init = init, ...)
                #)
				unperm_dist <- estimate_distr(sc_unperm, distr = distr, init = init, ...)
				if (distr == "norm"){
					final_distr <- list(mean = perm_dist$mean, sd = unperm_dist$sd)
				} else if (distr == "mnorm"){
					final_distr <- get_adj_mnorm(perm = perm_dist, unperm = unperm_dist)
				}
			} else {
				final_distr <- perm_dist
			}
			score <- scale_scores(score, method = output, param = final_distr, 
                                thresh = thresh)
		}
		R[,i] <- score
	}
	return(R) 
}

#' This function is going to take X and a vector of index for columns to get score for 
#' a certain set
#' @param X Matrix of processed data  
#' @param idx A vector size equals number of columns of X to indicate which taxa is part of which set 
get_score <- function(X, idx){
	# check type 
	if (is.matrix(X) == F){
		message("Coercing X to matrix")
		X <- as.matrix(X)
	}
	if (is.vector(idx) == F){
		message("Coercing idx to a numeric vector")
		idx <- as.vector(idx)
	}
	# get total columns and define size ans scale funciton 
	p <- ncol(X)
	size <- sum(idx == 1)
	scale <- sqrt(size * (p - size)/p)
	# calculate geometric mean 
	num <- geometricmeanRow(x = as.matrix(X[,idx == 1]))
	denom <- geometricmeanRow(x = as.matrix(X[,idx == 0]))
	# return the ilr like statistic 
	return(scale*(log(num/denom)))
}

#' This function is going to estimate from a list of supported distribution
#' @param data This is the data used to estimate a distribution 
#' @param distr The distribution to be fitted  
#' @param init The initialization (list)
#' @param ... Other paremteres passed to fitdistrplus or normalmixEM
estimate_distr <- function(data, distr = c("mnorm", "norm"), init, ...){
	# matching argument 
	distr <- match.arg(distr)
	dist <- tryCatch({
		if (missing(init) | is.null(init)){
			if (distr == "norm"){
				init <- NULL
			} else if (distr == "mnorm"){
				init <- list(lambda = NULL, mu = NULL, sigma = NULL)
			}
		}
		message(glue("Fitting the {d} distribution", d = distr))
		if (distr %in% c("norm")){
			fit <- fitdistrplus::fitdist(data, distr = distr, method = "mle", start = init)
			list(mean = fit$estimate[['mean']], sd = fit$estimate[['sd']])
		} else if (distr == "mnorm"){
			params <- rlist::list.append(init, x = data)
			params <- c(params, list(...)) # grabbing arguments to put in normalmixEM
			fit <- do.call(normalmixEM, params)
			list(mu = fit$mu, sigma = fit$sigma, lambda = fit$lambda)
		}
	}, error = function(cond){
		message("There were errors with the fitting procedure, adjusting the fitting parameters might help")
		message(cond)		
		return(NULL)
	})
	# if you cannot estimate the distribution of dist 
	if (is.null(dist)){
		return(NULL)
	} else {
		return(dist)
	}
}

# this function rescales the scores into significance (0,1), p-value, zscore or cdf values 
scale_scores <- function(scores, method = c("cdf","zscore", "pval", "sig"), param, thresh=0.05){
    # detect if parameter length is concordant with distribution type 
  	if (length(param) > 2){
  		f <- "pmnorm"
  	} else {
  		f <- "pnorm"
  	}
  	if(method %in% c("cdf", "sig","pval")){
  		param <- rlist::list.append(q = as.vector(scores), param)
  		scale <- do.call(f, param)
  		if (method == "pval"){
              scale <- 1 - scale
  		} else if (method == "sig"){
              scale <- 1 - scale
              scale <- ifelse(scale <= thresh, 1, 0)
  		}
  	} else if (method == "zscore"){
  		if (f == "pmnorm"){
  			mean <- get_mean(mu = param$mu, lambda = param$lambda)
  			sd <- get_sd(sigma = param$sigma, mu = param$mu, mean = mean, lambda = param$lambda)
  		} else if (f == "pnorm"){
  			mean <- param$mean
  			sd <- param$sd
  		}
  		scale <- (scores - mean) * 1/sd
  	}
  	return(scale)
}

# function to get the adjusted mixed normal using BFGS optimization 
get_adj_mnorm <- function(perm, unperm){
	perm_mean <- get_mean(mu = perm$mu, lambda = perm$lambda)
	unperm_mean <- get_mean(mu = unperm$mu, lambda = unperm$lambda)
	unperm_sd <- get_sd(sigma = unperm$sigma, mu = unperm$mu, 
						mean = unperm_mean, lambda = unperm$lambda)
	perm_sd <- get_sd(sigma = perm$sigma, mu = perm$mu, mean = perm_mean,
						lambda = perm$lambda)
	perm_skew <- Rfast::skew(rnormmix(n=1e4, lambda = perm$lambda, 
							sigma = perm$sigma, mu = perm$mu))
	# define objective function 
	obj_function <- function(sigma, mu, lambda, mean, sd, skew){
		s1 <- sigma[1]
		s2 <- sigma[2]
		m1 <- mu[1]
		m2 <- mu[2]
		l1 <- lambda[1]
		l2 <- lambda[2]
		estimate_sd <- sqrt((s1 + m1 - mean)*l1 + (s2 + m2 - mean)*l2)
		#estimate_skew <- Rfast::skew(rnormmix(n = 1e5, lambda = c(l1,l2), 
		#									sigma = c(s1,s2), mu = c(m1,m2)))
		#obj <- sqrt(sum((c(estimate_sd, estimate_skew) - c(sd, skew))^2))
		obj <- sqrt(sum((estimate_sd - sd)^2))
		return(obj)
	}

	opt <- optim(par = c(0.1,0.1), obj_function, mu = perm$mu, lambda = perm$lambda, 
					mean = perm_mean, sd = unperm_sd, skew = perm_skew)
	estim_sd <- get_sd(sigma = opt$par, lambda = perm$lambda, mu = perm$mu, mean = perm_mean)
	#estim_skew <- Rfast::skew(rnormmix(n=1e5, lambda = perm$lambda, mu = perm$mu, sigma = opt$par))
	print(glue("Total sd is {x} and estimated sd is {y}", x = unperm_sd, y = estim_sd))
	#print(glue("Total skew is {x} and estimated skew is {y}", x = perm_skew, y = estim_skew))
	param <- list(mu = perm$mu, sigma = opt$par, lambda = perm$lambda)
	return(param)
}

# cdf function of the mixture normal distribution 
#' @param q values presumably distributed via a mixture normal 
#' @param mu A vector of mu values equals to the number of components
#' @param sigma A vector of sigma values equals to the number of components
#' @param lambda A vector of lambda values equals to the number of components 
#' @param log A boolean indicating whether returning probabilities are in log format 
pmnorm <- function(q, mu, sigma, lambda, log = FALSE){
	q <- as.vector(q)
	n_components <- length(sigma)
	print(paste(n_components, "components!"))
	comp <- vector(mode = "list", length = n_components)
	for (i in 1:n_components){
		comp[[i]] <- lambda[i] * pnorm(q, mu[i], sigma[i], log.p = log)
	}
	return(Reduce("+", comp))
}


dmnorm <- function(x, mu, sigma, lambda, log=FALSE){
  x <- as.vector(x)
  n_components <- length(sigma)
  print(paste(n_components, "components!"))
  comp <- vector(mode = "list", length = n_components)
  for (i in 1:n_components){
    comp[[i]] <- lambda[i] * dnorm(x, mu[i], sigma[i], log = log)
  }
  return(Reduce("+", comp))
}


get_sd <- function(sigma, mu, mean, lambda){
  sqrt(sum((sigma + mu - mean)*lambda))
}

get_mean <- function(mu, lambda){
  sum(lambda * mu)
}