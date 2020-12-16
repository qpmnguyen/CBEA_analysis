# Title     : Unit testing for diff_ab functions
# Created by: Quang Nguyen
# Created on: 12/16/2020


library(testthat)
library(phyloseq)
library(DESeq2)
library(corncob)
source("../R/cilr.R")
source("../R/diff_ab_functions.R")
source("../R/utils.R")
source("../R/simulations.R")
context("Testing functionality of differential abundance functions")
dat <- zinb_simulation(n_samp = 200, s_rho = 0, spar = 0.2, eff_size = 2, n_inflate = 50, n_sets = 6)
phy <- sim2phylo(dat)

phy <- tax_glom(phy, taxrank = "GENUS")
phy <- transform_sample_counts(phy, function(x) ifelse(x == 0, 1, x))

# testing for deseq
deseq <- phyloseq_to_deseq2(phy, ~as.factor(group))
mod <- DESeq2::DESeq(deseq, test = "LRT", fitType = "local", reduced = ~1)
res <- DESeq2::results(mod)
sig_deseq <- res$pvalue
names(sig_deseq) <- as.vector(tax_table(phy)[rownames(res), "GENUS"])

# testing for corncob
mod <- differentialTest(phy, formula = ~ group, phi.formula = ~ group, formula_null = ~1,
                phi.formula_null = ~ group, test = "LRT", boot = FALSE,
                fdr_cutoff = 0.05)
sig_corncob <- mod$p
names(sig_corncob) <- as.vector(tax_table(phy)[names(sig_corncob), "GENUS"])


test_that("Testing main function", {
    res_corncob <- diff_ab(phy, method = "corncob", output = "pvalue")
    res_cilr <- diff_ab(phy, method = "cilr_wilcox", output = "pvalue")
    res_deseq2 <- diff_ab(phy, method = "deseq2", output = "pvalue")

    # test for correction
    expect_identical(res_corncob, sig_corncob)
    expect_identical(res_cilr, sig_cilr)
    expect_identical(res_deseq2, sig_deseq)

    # testing for output formatting
    expect_gt(length(unique(res_corncob)), 2)
    expect_equal(length(unique(diff_ab(phy, method = "corncob", output = "sig"))), 2)

})

test_that("model_interface testing",{
    # getting results from corncob_test and deseq2_test
    corncob_test <- model_interface(phy, method = "corncob", "GENUS")
    deseq2_test <- model_interface(phy, method = "deseq2", agg_level = "GENUS")
    cilr_test <- model_interface(phy, method = "cilr_wilcox", agg_level = "GENUS", pcount = 0, transform = NULL)

    # Testing for correct formatting
    expect_vector(corncob_test)
    expect_vector(deseq2_test)
    expect_vector(cilr_test)
    expect_named(corncob_test)
    expect_named(deseq2_test)
    expect_named(cilr_test)


    # Testing for throwing error
    expect_error(model_interface(phy, method = "cilr", agg_level = "GENUS"))

    # testing for corncob
    expect_identical(corncob_test, sig_corncob)
    # testing for deseq2
    expect_identical(deseq2_test, sig_deseq)

})