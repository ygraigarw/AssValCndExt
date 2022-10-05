# Associated values from conditional extreme value models : AssValCntExt

This repository contains MATLAB code and supplementary figures supporting the article

Estimation of associated values (AssVal) from conditional extreme (CndExt) value models by Towe, Randell, Kensler, Feld and Jonathan.

A preprint of the paper is available at https://ygraigarw.github.io/TowEA22_AssValCndExt.pdf

# Summary of paper

The design and reanalysis of offshore and coastal structures usually requires the estimation of return values for dominant metocean variables (such as significant wave height) and associated values for other variables (such as peak spectral period or wind speed) from a finite sample of data; these are typically estimated using extreme value analysis. Yet the parameters of extreme value models can only be estimated with error from finite data. Different choices available to summarise uncertain information about the characteristics of the tail of a multivariate distribution in a small number of summary statistics (such as return values and associated values) complicates their estimation, especially for small sample sizes: choices regarding the ordering of mathematical operations lead to estimators of return values and associated values with different finite sample bias and variance characteristics. The current work extends a previous study (Jonathan et al. 2021) into the performance of estimators for marginal return values in the presence of sampling uncertainty, to estimators of associated values based on the bivariate conditional extremes model (Heffernan and Tawn 2004) and competitors. Using a large designed simulation experiment, we explore the performance of combinations of 12 different estimators and three bivariate model candidates. The rich set of results from the simulation experiment are reported and explained in detail. Briefly: (a) calculation of associated values is only always feasible from small samples using 2 of the 12 estimators, which should be preferred; (b) estimators exploiting the median rather than the mean to summarise a distribution are more robust, and should also be preferred, especially for small sample sizes; (c) extreme value models incorporating appropriate descriptions of marginal and dependence provide better estimation of associated values for larger sample size; and (d) summarising the joint tail of metocean variables (in terms of return values and associated values) should be avoided where possible, in favour of probabilistic risk analysis of structural failure incorporating full uncertainty propagation.

# Supplementary figures

The directory SupplementaryFigures contains analogues of Figures 9 and 10 in the paper, for every possible combination of estimator of associated values (defined in Equations 6 and 7 of the paper) and sample size. Each figure shows the empirical distribution of fractional bias as a function of four design covariates, for the conditional extremes and simple linear regression (RgrEst) estimation schemes. 

In the box-whisker structure, thin horizontal lines in the box represent sample 25\%, 50\% and 75\% quantiles, and whiskers 2.5\% and 97.5\% quantiles. The thick horizontal line corresponds to the sample mean. When the whisker locations for a given plot panel lie outside the interval [-3,3], the vertical extent of the panel is restricted to [-3,3] for each of viewing. 

The paper provides a discussion.

# Software

The directory Code contains MATLAB functions allowing the full simulation study to be reproduced. Look at A_DriverPrlSmlPpr.m ("DRIVER for the PaRaLlel SiMuLation in the PaPeR") to get going. 

The full simulation takes of the order of 10 days with 40 parallel workers.

The driver file is currently set up to run a quick "check" analysis, so that a new user of the software can test that functions execute and data are saved appropriately.
