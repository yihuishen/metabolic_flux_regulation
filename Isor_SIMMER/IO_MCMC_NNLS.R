## perform MCMC-NNLS fitting of multiomics data to all possible reaction kinetics

library(ggplot2)
library(simmer) # https://github.com/shackett/simmer

options(stringsAsFactors = F)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

load('../../../result/rxnForms_filled.RData')

fitMetOnly <- T # T, only fit metabolites; F, include both enzyme and metabolites in fitting

{markov_pars <- list()
  markov_pars$sample_freq <- 20 #what fraction of markov samples are reported (this thinning of samples decreases sample autocorrelation)
  markov_pars$n_samples <- 1000 #how many total markov samples are desired
  markov_pars$burn_in <- 200 #how many initial samples should be skipped
} # markov parameters

if (!exists('errRxn')){
  if (!file.exists('../../../result/errRxn.txt')) {
    errRxn <- data.frame(rMech = character(0), reason = character(0))
  } else {
    errRxn <- read.delim('../../../result/errRxn.txt', sep = '\t')
  }
} # some reaction mechanisms may fail in the MCMC_NNLS

for (rxn in 1:length(rxnForms)) {

  if((rxn %% 10) == 1) {
    message(paste('Analyzing rxnForms',rxn,'-',rxn+9))
  }
  rm(list = ls()[!ls() %in% c('rxnForms','markov_pars','errRxn','rxn','fitMetOnly')])
  rxnSummary <- rxnForms[[rxn]]


  # for isozyme not detected, log2(enzyme) will be set to -Inf
  # since function 'lik_calc_fittedSD' considers all isozymes in nnls fitting
  fil <- apply(rxnSummary$enzymeComplexes,1,FUN = function(x){all(is.nan(x))})
  rxnSummary$enzymeComplexes[fil,] <- -Inf

  if (!length(rxnSummary$fit) == 0) {
    message(paste('reaction-mod',rxnSummary$listEntry,'result already exists'))
    next
  }

  if (rxnSummary$rxnID %in% gsub('-rm','',errRxn)) {
    message(paste('reaction-mod',rxnSummary$listEntry,' skipped according to errRxn list'))
    next
  }

  kinetically_differing_isoenzymes <- F

  ## ----evaluate_one_rMech--------------------------------------------------

  rxnEquations <- format_raw_equations(rxnSummary, kinetically_differing_isoenzymes)
  summarize_kinetic_parameters(rxnSummary, rxnEquations, kinetically_differing_isoenzymes)
  omic_data <- format_omic_data(kineticPars, all_species, rxnSummary, kinetically_differing_isoenzymes)
  if(fitMetOnly) {
    omic_data$enzyme_abund[,]=1
  }
  err1 <- try(finalize_reaction_equations(rxnEquations,
                                          all_species,
                                          kinetically_differing_isoenzymes),
              silent = T)

  kineticParPrior <- build_kinetic_parameter_priors(rxnSummary, kineticPars, omic_data)
  message(paste('running MCMC-NNLS for reaction-mod',rxnSummary$listEntry))
  err2 <- try(fit_reaction_equation_MCMC_NNLS(markov_pars,
                                              kineticPars,
                                              kineticParPrior,
                                              rxnEquations,
                                              all_species,
                                              omic_data,
                                              kinetically_differing_isoenzymes),
              silent = T)
  if (inherits(err1,'try-error') | inherits(err2,'try-error') ) {
    errRxn <- rbind(errRxn,data.frame(rMech = rxnSummary$listEntry, reason = 'MCMC-NNLS failed'))
    message('failed')
    next
  }

  #####
  ## fitted fluxes
  markov_fit <- matrix(NA, ncol = nrow(rxnSummary$flux), nrow = nrow(markov_par_vals))
  nnls_par_vals <- matrix(NA,
                          ncol = ncol(omic_data$enzyme_abund),
                          nrow = nrow(markov_par_vals))
  colnames(nnls_par_vals) <- colnames(omic_data$enzyme_abund)
  r = rep(NA, nrow(markov_par_vals))

  for (i in 1:nrow(markov_par_vals)){
    cur_par <- markov_par_vals[i,]
    n_c <- nrow(omic_data$met_abund)
    par_stack <- rep(1, n_c) %*% t(cur_par)
    colnames(par_stack) <- kineticPars$formulaName

    par_stack <- 2^par_stack
    occupancy_vals <- data.frame(omic_data$met_abund, par_stack)

    if(!(kinetically_differing_isoenzymes)){
      predOcc <- eval(rxnEquations[["l_occupancyExpression"]], occupancy_vals)
      # predict occupancy as a function of metabolites and kinetic constants based upon the occupancy equation
      enzyme_activity <- (predOcc %*% t(rep(1, sum(all_species$SpeciesType == "Enzyme"))))*2^omic_data$enzyme_abund
      # occupany of enzymes * relative abundance of enzymes
    }else{
      enzyme_activity <- NULL
      for(isoenzyme in names(rxnSummary$rxnForm)){
        predOcc <- eval(rxnEquations[["l_occupancyExpression"]][[isoenzyme]], occupancy_vals)
        enzyme_activity <- cbind(enzyme_activity, predOcc %*% t(rep(1, sum(occEqtn_complex_match$occEqtn == isoenzyme)))*2^omic_data$enzyme_abund[,colnames(omic_data$enzyme_abund) %in% occEqtn_complex_match$complex[occEqtn_complex_match$occEqtn == isoenzyme]])
      }
    }

    flx <-  (omic_data$flux$FVAmax + omic_data$flux$FVAmin)/2
    flux_fit <- nnls::nnls(enzyme_activity,flx)
    markov_fit[i,] <- flux_fit$fitted
    nnls_par_vals[i,] <- flux_fit$x
    r[i] <- cor(flux_fit$fitted,flx)
  }

  markov_fit <- as.data.frame(markov_fit)
  colnames(markov_fit) <- rownames(rxnSummary$flux)

  ## check the distribution of fitted flux with boundaries from data
  # set <- 'P4'
  # figure <- ggplot(markov_fit, aes_string(set)) +
  #   geom_area(stat='bin') +
  #   geom_vline(xintercept = rxnSummary$flux[set,'FVAmin'])+
  #   geom_vline(xintercept = rxnSummary$flux[set,'FVAmax'])
  # figure

  pos <- which(lik_track == max(lik_track))
  if (length(pos) == markov_pars$n_samples) {
    message('cannot find maximum likelihood from MCMC-NNLS')
    errRxn <- rbind(errRxn,data.frame(rMech = rxnSummary$listEntry,
                                      reason = 'cannot find maximum likelihood from MCMC-NNLS'))
    next
  }
  fit_5 <- numeric(length= ncol(markov_fit))
  fit_95 <- numeric(length= ncol(markov_fit))

  for (i in 1:ncol(markov_fit)) {
    tmp <- sort(markov_fit[,i])
    fit_5[i] <- tmp[round(nrow(markov_fit)*0.05)]
    fit_95[i] <- tmp[round(nrow(markov_fit)*0.95)]
  }

  rxnSummary$kineticPars <- kineticPars
  rxnSummary$rxnEquations <- rxnEquations
  rxnSummary$omic_data <- omic_data

  rxnSummary$fit <- list()
  rxnSummary$fit$prior <- kineticParPrior
  rxnSummary$fit$r <- r
  rxnSummary$fit$fit_5 <- fit_5
  rxnSummary$fit$fit_95 <- fit_95
  rxnSummary$fit$markov_par_vals <- markov_par_vals
  rxnSummary$fit$markov_fit <- markov_fit
  rxnSummary$fit$lik_track <- lik_track
  rxnSummary$fit$nnls_par_vals <- nnls_par_vals

  rxnForms[[rxn]] <- rxnSummary

  message('successful!')

  if(((rxn %% 10) == 0)) {
    message(paste('Saving rxnForms',ifelse((rxn %% 10) == 0, rxn-9, rxn - (rxn %% 10) + 1),'-',rxn))
    save(rxnForms, file = ifelse(fitMetOnly,'../../../result/rxnForms_fitted_enzyme_excluded.RData',
                                 '../../../result/rxnForms_fitted.RData'))
    write.table(unique(errRxn),file='../../../result/errRxn.txt',
                sep = '\t', quote = F, row.names = F, col.names = F)
  } # save every 10 analysis

}

message(paste('Saving rxnForms',ifelse((rxn %% 10) == 0, rxn-9, rxn - (rxn %% 10) + 1),'-',rxn))
save(rxnForms, file = ifelse(fitMetOnly,'../../../result/rxnForms_fitted_enzyme_excluded.RData',
                             '../../../result/rxnForms_fitted_included.RData'))
write.table(unique(errRxn),file='../../../result/errRxn.txt',
            sep = '\t', quote = F, row.names = F, col.names = F)

