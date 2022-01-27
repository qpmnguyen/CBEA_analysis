# Simulation design
library(tidyverse)
library(glue)
library(phyloseq)
library(MASS)
library(VGAM)
library(furrr)
library(mia)
library(BiocSet)
library(Matrix)

#' Shorthand to create parameters
create_parameters <- function(params) {
    par <- cross_df(params)
    par <- par %>%
        mutate(id = seq(nrow(par))) %>%
        group_by(id) %>%
        nest()
    par <- par %>% transmute(param = data)
    return(par)
}


#' Simulating according to zero inflated negative binomial distribution
#' @param n_samp Number of samples
#' @param spar Additive sparsity
#' @param b_rho Baseline inter-taxa correlation
#' @param s_rho Correlation within the set
#' @param rho_ratio Ratio of correlation between the set and baseline
#' @param n_tax Number of taxa
#' @param n_inflate Number of taxa with inflated counts
#' @param n_sets Number of sets to simulate
#' @param prop_set_inflate Number of proportions of sets that are inflated
#' @param prop_inflate The number of differentially abundant taxa within a set that is inflated
#' @param samp_prop The proportion of samples have an inflated taxa
#' @param eff_size Effect size that is a multiplier to the base mu of a sample
#' @param method Method of simulation. Can be "compensation", "non-compensation" or "normal"
#' @param parameters Path to rds file containing estimated values
#' @param vary_params Whether we stochastically draw parameter values for each feature
#' TODO: More than one correlation structure
#' TODO: Adjust rho_ratio calculation when b_rho = 0 (no background correlation)
#' TODO: Add shuffling to make sure that everything is randomized
#' TODO: Add a way to sample from empirical distribution of zinb values
zinb_simulation <- function(n_samp, spar, s_rho, eff_size,
                            b_rho = 0, n_tax = 300, n_inflate = 50, n_sets = 1,
                            prop_set_inflate = 1, prop_inflate = 1, samp_prop = 0.5,
                            method = "compensation", vary_params = TRUE, parameters = NULL) {
    # generate the the diagnonal matrix
    sigma <- diag(n_tax)
    sigma[sigma == 0] <- b_rho
    print(dim(sigma))
    set_sigma <- sigma[1:n_inflate, 1:n_inflate]
    set_sigma[set_sigma != 1] <- s_rho
    sigma[1:n_inflate, 1:n_inflate] <- set_sigma

    # First create mvrnorm variables with correlation set by sigma
    margins <- pnorm(MASS::mvrnorm(n = n_samp, mu = rep(0, n_tax), Sigma = sigma))
    # Second, set marginals
    # default marginals for negative binomial is from size of 0.595 and mu of 6.646 from HMP data
    true_size <- round(n_inflate * prop_inflate, 0)

    # Create parameters based on situation
    if (vary_params == T) {
        if (is.null(parameters)) {
            message("Randomly sample means from 1 to 10 and sizes from 1 to 5")
            means <- runif(n_tax, 3, 5)
            sizes <- runif(n_tax, 1, 3)
        } else {
            estimated <- readRDS(file = parameters)
            means <- sample(estimated$mean, size = n_tax, replace = T)
            sizes <- sample(estimated$size, size = n_tax, replace = T)
        }
    } else {
        message("Setting mean to be constant at 3.06 and size at 1.67 estimated from HMP data...")
        means <- rep(3.05, n_tax)
        sizes <- rep(1.67, n_tax)
    }


    # first n_samp * samp_prop samples will always be inflated
    inf_size <- round(n_samp * samp_prop, 0)

    # Number of inf_taxa
    inf_tax <- round(n_inflate * n_sets * prop_set_inflate, 0)

    # inflating samples
    suppressMessages(
        inf_samples <- map_dfc(seq(n_tax), .f = function(.x) {
            if (method == "compensation") {
                prop <- 1 / (eff_size + 1)
                # randomly select first prop of those selected to be upregulated
                idx <- sample(seq(inf_tax), size = round(prop * inf_tax, 0), replace = F)
                remainder <- seq(inf_tax)[-idx]
                a <- sum(means[idx])
                b <- sum(means[remainder])
                # randomly select second prop of those selected to be downregulated
                if (.x %in% idx) {
                    result <- qnbinom(
                        p = margins[seq(inf_size), .x], size = sizes[.x],
                        mu = means[.x] * eff_size
                    )
                } else if (.x %in% remainder) {
                    result <- qnbinom(
                        p = margins[seq(inf_size), .x], size = sizes[.x],
                        mu = means[.x] * ((a / b) * (1 - eff_size) + 1)
                    )
                } else {
                    result <- qnbinom(
                        p = margins[seq(inf_size), .x], size = sizes[.x],
                        mu = means[.x]
                    )
                }
            } else if (method == "no_compensation") {
                # randomly select first prop of those selected to be upregulated
                idx <- sample(seq(inf_tax), size = inf_tax / 2, replace = F)
                remainder <- seq(inf_tax)[-idx]
                # randomly select second prop of those selected to be downregulated
                if (.x %in% idx) {
                    result <- qnbinom(
                        p = margins[seq(inf_size), .x], size = sizes[.x],
                        mu = means[.x] * eff_size
                    )
                } else if (.x %in% remainder) {
                    result <- qnbinom(
                        p = margins[seq(inf_size), .x], size = sizes[.x],
                        mu = means[.x] / eff_size
                    )
                } else {
                    result <- qnbinom(
                        p = margins[seq(inf_size), .x], size = sizes[.x],
                        mu = means[.x]
                    )
                }
            } else if (method == "normal") {
                if (.x %in% seq(inf_tax)) {
                    result <- qnbinom(
                        p = margins[seq(inf_size), .x], size = sizes[.x],
                        mu = means[.x] * eff_size
                    )
                } else {
                    result <- qnbinom(
                        p = margins[seq(inf_size), .x], size = sizes[.x],
                        mu = means[.x]
                    )
                }
            }
            return(result)
        })
    )
    # not inflated samples
    suppressMessages(
        notinf_samples <- map_dfc(seq(n_tax), .f = function(.x) {
            result <- qnbinom(
                p = margins[-seq(inf_size), .x], size = sizes[.x],
                mu = means[.x]
            )
            return(result)
        })
    )
    message("Completed loop!")
    if (inf_size == n_samp) {
        abundance <- inf_samples
    } else {
        abundance <- rbind(inf_samples, notinf_samples)
    }
    abundance <- as.matrix(abundance)
    if (!is.null(spar)) {
        zeroes <- rbinom(length(abundance), size = 1, prob = 1 - spar)
        abundance <- abundance * zeroes
    }
    label <- c(rep(1, inf_size), rep(0, n_samp - inf_size))

    colnames(abundance) <- glue("Tax{i}", i = seq(n_tax))
    rownames(abundance) <- glue("Samp{i}", i = seq(n_samp))
    if (n_sets > 1) {
        A <- diag(n_sets)
        vec <- as.matrix(rep(1, n_inflate))
        A <- kronecker(A, vec)
    } else {
        # TODO if n_sets * n_inflate != n_tax, then there are issues.
        message("Only one set!")
        A <- rep(0, n_tax)
        A[1:n_inflate] <- 1
        A <- as.matrix(A)
    }
    colnames(A) <- glue("Set{i}ss", i = 1:n_sets)
    rownames(A) <- colnames(abundance)
    sets_inf <- rep(0, n_sets)
    sets_inf[seq(round(n_sets * prop_set_inflate, 0))] <- 1
    # convert A into BiocSet and abundance into physeq with label and sets_inf
    set_list <- map(seq_len(ncol(A)), ~{
        rownames(A)[which(A[,.x] == 1)]
    })
    names(set_list) <- colnames(A)
    set <- BiocSet::BiocSet(set_list)
    
    # convert abund into physeq
    obj <- TreeSummarizedExperiment(list(Counts = t(as.data.frame(abundance))))
    

    output <- list(obj = obj, set = set, label = label, sets_inf = sets_inf)
    return(output)
}

#'  @param modeltype What type of simulations
#'  @param snr Signal to noise ratio (or effect size)
#'  @param sat Saturation proportion (proportion of sets associated with outcome)
#'  @param ... parameters to pass to the zinb simulation function
sim_prediction <- function(type = c("regr", "classif"), snr = 2, sat = 0.1, ...) {
    type <- match.arg(type)
    # handle defaults
    def <- list(
        n_samp = 300, spar = 0.2, s_rho = 0.5, eff_size = 1,
        b_rho = 0, n_tax = 2000, n_inflate = 50, n_sets = 40,
        prop_set_inflate = 0.5, prop_inflate = 1, samp_prop = 0.5,
        method = "normal", vary_params = FALSE, parameters = NULL
    )
    if (missing(...)) {
        args <- def
    } else {
        sup <- list(...)
        args <- merge_lists(defaults = def, supplied = sup)
    }
    if (type == "classif") {
        sd_beta <- sqrt(2)
    } else if (type == "regr") {
        sd_beta <- 1
    }
    # generate baseline data based on arguments
    baseline <- do.call(zinb_simulation, args)
    # index of samples that are related to the outcome based on model saturation
    set_names <- colnames(baseline$A)
    # baseline$X <- baseline$X %>% as.matrix() %>% acomp() %>% unclass() %>% as.data.frame()
    if (sat > 0) {
        # First, let's sample from the sets the size of saturation
        n_core_sets <- round(length(set_names) * sat, 0)
        if (n_core_sets %% 2 != 0) {
            n_core_sets <- n_core_sets + 1
        }
        important_sets <- sample(seq(length(set_names)), size = n_core_sets)
        # Then let's generate a list of the sets that are determined to be important, and then
        # retrieve the taxa that are important
        n_pos <- n_core_sets / 2
        index_list <- lapply(important_sets, function(x) {
            as.vector(which(baseline$A[, x] == 1))
        })
        y_mean <- rep(0, nrow(baseline$X))
        beta_0 <- 6 / sqrt(10)
        while (sd(y_mean) == 0) { # remake y until the standard deviation is larger than 0
            # beta are beta values per taxa
            beta <- rep(0, ncol(baseline$X))
            # beta_sets are beta valutes per set
            beta_sets <- rep(0, n_core_sets)
            beta_sets[1:n_pos] <- runif(n_pos, 1.5, 2)
            beta_sets[-c(1:n_pos)] <- runif(n_pos, -2, -1.5)
            for (i in seq(length(index_list))) {
                # mat <- .sparseDiagonal(length(index_list[[i]]), shape = "s")
                # mat[mat == 0] <- beta_corr
                beta[index_list[[i]]] <- beta_sets[i]
            }
            y_mean <- beta_0 + as.matrix(baseline$X) %*% beta
        }
        y <- y_mean + rnorm(nrow(baseline$X), mean = 0, sd = sd(y_mean) * (1 / snr))
        if (type == "classif") {
            prob <- plogis(y) # logistic function
            n_samp <- length(prob)
            lower_bound <- round(0.4 * n_samp, 0)
            upper_bound <- round(0.6 * n_samp, 0)
            b_y <- rep(0, n_samp)
            iter <- 0
            while (iter < 2e4 & (sum(b_y == 1) < lower_bound | sum(b_y == 1) > upper_bound)) {
                iter <- iter + 1
                b_y <- as.factor(rbinom(n_samp, size = 1, prob = prob))
            }
            if (iter >= 2e4) {
                message("Did not converge on even classes")
            }
            y <- b_y
        }
    } else {
        y <- rnorm(nrow(baseline$X))
        beta <- rep(0, ncol(baseline$X))
    }
    print("End generation")
    output <- list(outcome = y, predictors = list(X = baseline$X, A = baseline$A), beta = beta)
}
