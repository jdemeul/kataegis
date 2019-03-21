### Code for generating the background odds using simulated events

#STATICS
library(readr)
library(ggplot2)
PCFDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20180309_pcf_rerun/"
CCLUSTDIR <- "/srv/shared/vanloo/ICGC-consensus-clustering/consensus_subclonal_reconstruction_v1.3_20190221/"
CTIMINGDIR <- "/srv/shared/vanloo/ICGC-consensus-clustering/consensus_subclonal_reconstruction_v1.3_20190221_probgain/"
CCCFDIR <- "/srv/shared/vanloo/ICGC-consensus-clustering/consensus_subclonal_reconstruction_v1.3_20190221_mutccf/"
KATRESULTSFILE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_Kataegis_Results_all.txt"
SIMDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/simulations/"
NSIMS <- 10000

source(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/kataegis_functions.R", local = T)


# katcalls <- allpcfout_clean
## TEMP: only LOH events
allpcfout_clean <- read.delim(file = KATRESULTSFILE, as.is = T)
allpcfout_clean_inform <- allpcfout_clean[allpcfout_clean$informative > 0 & allpcfout_clean$is_preferred, ]
# allpcfout_clean_inform <- allpcfout_clean[allpcfout_clean$informative > 0 & allpcfout_clean$is_preferred & allpcfout_clean$nMin == 0, ]
samplesnhist <- allpcfout_clean_inform[!duplicated(allpcfout_clean_inform$sample), c("sample", "histology")]


### start function defs
# for every sample, do
simulate_events_forodds <- function(sample_id, pcfdir = PCFDIR, cclustdir = CCLUSTDIR, cccfdir = CCCFDIR, ctimingdir = CTIMINGDIR, katcalls, nreps = 10, simdir = SIMDIR) {
  print(sample_id)
  # katcalls in sample
  # must have made sure that no samples without any informative foci are excluded from input samplelist
  katcallsub <- katcalls[katcalls$sample == sample_id, ]
  katcallsizes <- katcallsub$informative

  # pcf data
  pcf_file <- list.files(path = pcfdir, pattern = sample_id, full.names = T)
  muts_all_pcf <- read.delim(file = pcf_file, header = T, sep = "\t", as.is = T,
                             colClasses = c("character", "integer", rep("character", 2), "integer"))
  colnames(muts_all_pcf) <- c("chromosome", "pos", "ref", "alt", "foci")
  muts_all_pcf$chromosome <- as.character(muts_all_pcf$chromosome)
  
  # pcf_file <- list.files(path = pcfdir, pattern = sample_id, full.names = T)
  # muts_all_pcf <- read.delim(file = pcf_file, header = T, sep = "\t", as.is = T,
  #                            colClasses = c("character", "integer", rep("character", 3),
  #                                           rep(c("integer", "integer", "numeric"), 2),
  #                                           "integer", "integer", rep("numeric", 4)))
  
  # clustering consensus data
  cassignments <- read.delim(file = gzfile(file.path(cclustdir, paste0(sample_id, "_cluster_assignments.txt.gz"))), header = T, sep = "\t", as.is = T)
  cccfs <- read.delim(file = gzfile(file.path(cccfdir, paste0(sample_id, "_mutation_ccf.txt.gz"))), header = T, sep = "\t", as.is = T)
  ctiming <- read.delim(file = gzfile(file.path(ctimingdir, paste0(sample_id, "_prob_gained.txt.gz"))), header = T, sep = "\t", as.is = T)
  colnames(ctiming) <- c("chromosome", "position", "mut_type", "timing", "chromosome2", "position2", "svid", "prob_clonal_early", "prob_clonal_late", "prob_subclonal")
  
  if (all(c(nrow(cassignments), nrow(cccfs), nrow(ctiming)) == nrow(cassignments)) ) {
    mutsdf <- cbind(cassignments[ , c("chromosome", "position", "mut_type", "cluster_1")],
                    ctiming[ , c("timing", "prob_clonal_early", "prob_clonal_late", "prob_subclonal")],
                    cccfs[ , c("major_cn", "minor_cn", "mult")])
  } else {
    print("simulate_events_forodds: should not get here")
    mutsdf <- merge(x = cassignments[ , c("chromosome", "position", "mut_type", "cluster_1")], y = ctiming[ ,c("chromosome", "position", "timing", "prob_clonal_early", "prob_clonal_late", "prob_subclonal")], by = c("chromosome", "position"))
    mutsdf <- merge(x = mutsdf, y = cccfs[ , c("chromosome", "position", "major_cn", "minor_cn", "mult")], by = c("chromosome", "position"))
    mutsdf <- mutsdf[!duplicated(paste0(mutsdf$chromosome, "_", mutsdf$position, "_", mutsdf$cluster_1)),]
  }
  
  csubcl <- read.delim(file = gzfile(file.path(cclustdir, paste0(sample_id, "_subclonal_structure.txt.gz"))), header = T, sep = "\t", as.is = T)
  
  # subsetting for background
  mutsdf <- mutsdf[!paste(mutsdf$chromosome, mutsdf$position, sep = "_") %in% paste(muts_all_pcf$chromosome, muts_all_pcf$pos, sep = "_") &
                     mutsdf$mut_type == "SNV" & !is.na(mutsdf$cluster_1) & !is.na(mutsdf$prob_clonal_early), ]
  
  # fixing/regularising pGain/pSingle/pSubclonal
  # values < 1e-15 (approx .Machine$double.eps) set to 1e-15 and renormalised such that pSingle + pGain + pSub = 1
  mutsdf$prob_clonal_early <- ifelse(mutsdf$major_cn > 1 & mutsdf$prob_clonal_early < 1e-15, 1e-15, mutsdf$prob_clonal_early)
  mutsdf$prob_clonal_late <- ifelse(mutsdf$prob_clonal_late < 1e-15, 1e-15, mutsdf$prob_clonal_late)
  # if >= 1 subclone present, also regularise these probs
  if (nrow(csubcl) > 1) {
    mutsdf$prob_subclonal <- ifelse(mutsdf$prob_subclonal < 1e-15, 1e-15, mutsdf$prob_subclonal)
  }
  # renorm
  mutsdf[ , c("prob_clonal_early", "prob_clonal_late", "prob_subclonal")] <- mutsdf[ , c("prob_clonal_early", "prob_clonal_late", "prob_subclonal")] / rowSums(mutsdf[ , c("prob_clonal_early", "prob_clonal_late", "prob_subclonal")])

  # clean up
  rm(cassignments, cccfs, ctiming, muts_all_pcf)
  rownames(mutsdf) <- NULL
  
  # get samples for all nreps and all events in one go
  ## get idxs of nogain/gain/gainloh variants
  nogainidxs <- which(mutsdf$major_cn < 2)
  gainidxs <- which(mutsdf$major_cn > 1 & mutsdf$minor_cn > 0)
  gainlohidxs <- which(mutsdf$major_cn > 1 & mutsdf$minor_cn == 0)
  
  # if no variants of certain type, set to NA to avoid error in sample()
  if (length(nogainidxs) == 0) nogainidxs <- NA
  if (length(gainidxs) == 0) gainidxs <- NA
  if (length(gainlohidxs) == 0) gainlohidxs <- NA

  # If no matching variants where required for sampling (should be rare), just sample from the lot
  # This could slightly increase bias while reducing variance
  samplesize <- sum(katcallsizes)*nreps
  if ((sum(katcallsub$nNoGain) > 0 & is.na(nogainidxs)) ||
      (sum(katcallsub$nGain) > 0 & is.na(gainidxs)) ||
      (sum(katcallsub$nGainLOH) > 0 & is.na(gainlohidxs))) {
    print(paste0(sample_id, ": using entire background"))
    sampleidxs <- sample(x = 1:nrow(mutsdf), size = samplesize, replace = T)
  } else {
    # construct logicals indicating whether matched SNVs should be sampled from a nogain/gain/gainloh region
    is_gain_snv <- rep(x = rep(x = rep(x = c(F,T), times = length(katcallsizes)), times = c(rbind(katcallsub$nNoGain, katcallsub$nGain + katcallsub$nGainLOH))), times = nreps)
    is_loh_snv <- rep(x = rep(x = rep(x = c(F,T), times = length(katcallsizes)), times = c(rbind(katcallsub$nNoGain + katcallsub$nGain, katcallsub$nGainLOH))), times = nreps)
    # sample accordingly
    sampleidxs <- ifelse(is_gain_snv, 
                         ifelse(is_loh_snv, sample(x = gainlohidxs, size = samplesize, replace = T), sample(x = gainidxs, size = samplesize, replace = T)),
                         sample(x = nogainidxs, size = samplesize, replace = T))
  }
  
  sampledmuts <- mutsdf[sampleidxs, c("cluster_1", "timing", "prob_clonal_early", "prob_clonal_late", "prob_subclonal", "major_cn", "minor_cn", "mult")]
  # compute required "summary stats"
  sampledevents <- by(data = sampledmuts, INDICES = rep(1:(length(katcallsizes)*nreps), rep(x = katcallsizes, times = nreps)), FUN = get_odds_input_sims)
  # turn in to nreps x dataframes
  sampleddfs <- lapply(X = 1:nreps, FUN = function(i, n, evlist) do.call(rbind, evlist[((i-1)*n+1):(i*n)]), n = length(katcallsizes), evlist = sampledevents)
  
  saveRDS(object = sampleddfs, file = file.path(simdir, paste0(sample_id, "_oddssim.rds")))
  return(NULL)
  # return(sampleddfs)
}


## get the probs from the df with muts corresponding to a single simulated event
get_odds_input_sims <- function(df) {
  # compute normalised likelihoods
  posterior <- c(p_clonal = exp(sum(log(df$cluster_1))), p_subclonal = exp(sum(log(1 - df$cluster_1))) )
  posterior <- posterior/sum(posterior)
  
  # compute "probabilities" of clonal early/late/NA
  # classes_if_clonal <- ifelse(df$timing != "subclonal", df$timing,
  #                             ifelse(df$major_cn < 2, "clonal [NA]",
  #                                    ifelse(df$mult >= 2, "clonal [early]",
  #                                           ifelse(df$minor_cn == 0, "clonal [late]", "clonal [NA]"))))
  # classes_if_clonal <- factor(x = classes_if_clonal, levels = c("clonal [NA]", "clonal [early]", "clonal [late]"), labels = c("w_clonal_NA", "w_clonal_early", "w_clonal_late"))
  # w_classes_if_clonal <- c(by(data = df$cluster_1, INDICES = classes_if_clonal, FUN = sum)) / sum(df$cluster_1)
  # w_classes_if_clonal[is.na(w_classes_if_clonal)] <- 0
  earlylateprobs <- get_probs_earlylate(df)
  
  return(c(posterior, earlylateprobs))
}


## get the odds per tumour type from a single simulated eventsdf
get_sim_odds <- function(numer, denom, is_powered, histofactor, weights) {
  numer <- numer[is_powered]
  denom <- denom[is_powered]
  # sample the simulated events ~ weight (to emulate bootstrap iteration)
  simodds <- c(by(data = data.frame(numer = numer, denom = denom, weights = weights), INDICES = histofactor, FUN = odds_helper))
  return(simodds)
}


## get the internally normalised odds per tumour type from a single simulated eventsdf and the observed eventsdf
get_odds_internalnorm <- function(numer, denom, is_powered, histofactor, weights, numerobs, denomobs) {
  numer <- numer[is_powered]
  denom <- denom[is_powered]
  numerobs <- numerobs[is_powered]
  denomobs <- denomobs[is_powered]
  # sample the simulated events ~ weight (to emulate bootstrap iteration)
  simodds <- c(by(data = data.frame(numer = numer, denom = denom, numerobs = numerobs, denomobs = denomobs, weights = weights), INDICES = histofactor, FUN = odds_helper_internalnorm))
  return(simodds)
}


## helper function to compute the odds using a resampled (weighted) version of the eventsdf for a single tumour type  
odds_helper <- function(df) {
  f <- sample(x = 1:nrow(df), size = nrow(df), replace = T, prob = df$weights)
  odds <- sum(df[f, "numer"], .5) / sum(df[f, "denom"], .5)
  return(odds)
}


## helper function to compute the odds using a resampled (weighted) version of the eventsdf and the observed eventsdf for a single tumour type
odds_helper_internalnorm <- function(df) {
  f <- sample(x = 1:nrow(df), size = nrow(df), replace = T, prob = df$weights)
  odds <- (sum(df[f, "numerobs"], .5) / sum(df[f, "numer"], .5)) / (sum(df[f, "denomobs"], .5) / sum(df[f, "denom"], .5))
  return(odds)
}

#### end functions


##### simulation of events ~ observations
# debug(simulate_events_forodds)
# SAMPLE <- "6cfce053-bfd6-4ca0-b74b-b2e4549e4f1f"
# simlist <- lapply(X = SAMPLE, FUN = simulate_events_forodds, katcalls = allpcfout_clean_inform, nreps = NSIMS)
# simlist <- lapply(X = samplesnhist$sample, FUN = simulate_events_forodds, katcalls = allpcfout_clean_inform, nreps = NSIMS)
# simlist <- mclapply(X = samplesnhist$sample[1], FUN = simulate_events_forodds, katcalls = allpcfout_clean_inform, nreps = NSIMS, mc.preschedule = T, mc.cores = 1)
mclapply(X = samplesnhist$sample, FUN = simulate_events_forodds, katcalls = allpcfout_clean_inform, nreps = NSIMS, simdir = SIMDIR, mc.preschedule = F, mc.cores = 18)

# simlist <- mclapply(X = samplesnhist$sample, FUN = simulate_events_forodds, katcalls = allpcfout_clean_inform, nreps = NSIMS, mc.preschedule = F, mc.cores = 16)
simlist <- lapply(X = samplesnhist$sample, FUN = function(sampleid, simdir) readRDS(file = file.path(simdir, paste0(sampleid, "_oddssim.rds"))), simdir = SIMDIR)
# reformat simlist from a list of samples to a list of simulations
simlist <- lapply(X = 1:NSIMS, FUN = function(i, sims) as.data.frame(do.call(rbind, lapply(sims, '[[', i))), sims = simlist)



###### Try internal normalisation of the odds ratio using the matched simulated events
# get required info from original dataframe of results
is_powered <- allpcfout_clean_inform$nrpcc >= 10
eventweights <- c((table(allpcfout_clean_inform[is_powered, "sample"])^-1)[allpcfout_clean_inform[is_powered, "sample"]])
histofactor <- allpcfout_clean_inform[is_powered, "histology"]

# compute the odds
simodds_cvs_intnorm <- do.call(rbind, lapply(X = simlist, FUN = function(x) get_odds_internalnorm(numer = x$p_clonal, denom = x$p_subclonal,
                                                                                                  numerobs = allpcfout_clean_inform$p_clonal, denomobs = allpcfout_clean_inform$p_subclonal,
                                                                                                  is_powered = is_powered, histofactor = histofactor, weights = eventweights)))
simodds_cvs_intnorm <- t(apply(X = simodds_cvs_intnorm, MARGIN = 2, FUN = function(x) quantile(x = x, probs = c(.025, .5, .975))))
colnames(simodds_cvs_intnorm) <- c("lower", "median", "upper")

# write out + visualise
write.table(x = simodds_cvs_intnorm, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_background_timing_odds_clonalVsubclonal_intnorm_cnstrat.txt", quote = F, sep = "\t")
p1 <- ggplot(data = as.data.frame(simodds_cvs_intnorm), mapping = aes(x = rownames(simodds_cvs_intnorm))) + geom_pointrange(mapping = aes(y = median, ymin = lower, ymax = upper)) + theme(axis.text.x = element_text(angle = 90)) + scale_y_log10(breaks = c(0.01, 0.1,1,10), labels = c(0.01,0.1,1,10))
p1


##### CLONAL EARLY v CLONAL LATE
# line below was quick patch
# allpcfout_clean_inform[, c("w_clonal_early", "w_clonal_late", "w_clonal_NA")] <- allpcfout_clean_inform[, c("w_clonal_early", "w_clonal_late", "w_clonal_NA")] / rowSums(allpcfout_clean_inform[ , c("w_clonal_early", "w_clonal_late", "w_clonal_NA")])
eventweights <- c((table(allpcfout_clean_inform$sample)^-1)[allpcfout_clean_inform$sample])
histofactor <- allpcfout_clean_inform$histology

# compute the odds
# simodds_evl_intnorm <- do.call(rbind, lapply(X = simlist, FUN = function(x) get_odds_internalnorm(numer = x$w_clonal_early*x$p_clonal, denom = x$w_clonal_late*x$p_clonal,
#                                                                                                   numerobs = allpcfout_clean_inform$w_clonal_early*allpcfout_clean_inform$p_clonal,
#                                                                                                   denomobs = allpcfout_clean_inform$w_clonal_late*allpcfout_clean_inform$p_clonal,
#                                                                                                   is_powered = T, histofactor = histofactor, weights = eventweights)))
simodds_evl_intnorm <- do.call(rbind, lapply(X = simlist, FUN = function(x) get_odds_internalnorm(numer = x$p_early, denom = x$p_late,
                                                                                                  numerobs = allpcfout_clean_inform$p_early,
                                                                                                  denomobs = allpcfout_clean_inform$p_late,
                                                                                                  is_powered = T, histofactor = histofactor, weights = eventweights)))
sum(!complete.cases(simodds_evl_intnorm))
# which(!complete.cases(simodds_evl_intnorm))
# which(!complete.cases(simlist[[which(!complete.cases(simodds_evl_intnorm))]]))
simodds_evl_intnorm <- t(apply(X = simodds_evl_intnorm, MARGIN = 2, FUN = function(x) quantile(x = x, probs = c(.025, .5, .975), na.rm = T)))
colnames(simodds_evl_intnorm) <- c("lower", "median", "upper")

# write out + visualise
write.table(x = simodds_evl_intnorm, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_background_timing_odds_earlyVlate_intnorm_cnstrat.txt", quote = F, sep = "\t")
p2 <- ggplot(data = as.data.frame(simodds_evl_intnorm), mapping = aes(x = rownames(simodds_evl_intnorm))) + geom_pointrange(mapping = aes(y = median, ymin = lower, ymax = upper)) + theme(axis.text.x = element_text(angle = 90)) + scale_y_log10(breaks = c(0.1,1,10,100), labels = c(0.1,1,10,100))
p2

# save(list = c("simlist", "allpcfout_clean_inform"), file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_background_timing_odds_simlist+katcalls.RData")
# rm(simlist)
load("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_background_timing_odds_simlist+katcalls.RData")




# ##### temp visualisation stuff
# cvs <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/20180303_background_timing_odds_clonalVsubclonal.txt", as.is = T)
# cvs_loh <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/20180303_background_timing_odds_clonalVsubclonal_LOHonly.txt", as.is = T)
# cvs_intnorm <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/20180303_background_timing_odds_clonalVsubclonal_intnorm.txt", as.is = T)
# cvs_intnorm_loh <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/20180303_background_timing_odds_clonalVsubclonal_intnorm_LOHonly.txt", as.is = T)
#
# evl <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/20180303_background_timing_odds_earlyVlate.txt", as.is = T)
# evl_loh <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/20180303_background_timing_odds_earlyVlate_LOHonly.txt", as.is = T)
# evl_intnorm <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/20180303_background_timing_odds_earlyVlate_intnorm.txt", as.is = T)
# evl_intnorm_loh <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/20180303_background_timing_odds_earlyVlate_intnorm_LOHonly.txt", as.is = T)
#
# timres <- list(cvs, cvs_loh, cvs_intnorm, cvs_intnorm_loh, evl, evl_loh, evl_intnorm, evl_intnorm_loh)
# timres <- lapply(timres, function(x) {
#   y <- cbind.data.frame(x, histology = rownames(x))
#                  rownames(y) <- NULL
#                  return(y)})
# timresdf <- do.call(rbind.data.frame, timres)
# timresdf$cvs <- rep(c(T, T, T, T, F, F, F, F), sapply(timres, nrow))
# timresdf$intnorm <- rep(c(F, F, T, T, F, F, T, T), sapply(timres, nrow))
# timresdf$lohonly <- rep(c(F, T, F, T, F, T, F, T), sapply(timres, nrow))
# timresdf$grp <- interaction(timresdf$cvs, timresdf$intnorm, timresdf$lohonly)
#
# library(ggplot2)
# p5 <- ggplot(data = timresdf[timresdf$cvs, ], mapping = aes(x = histology)) + geom_pointrange(mapping = aes(x = histology, y = mean, ymin = lower, ymax = upper, colour = interaction(intnorm, lohonly), alpha = .5), position = position_dodge(width = .75)) + theme(axis.text.x = element_text(angle = -90))
# p5



##################### Replacing the Venn diagrams
load("/srv/shared/vanloo/home/mfittall/ICGC/Chromoplexy/Final_2019_2_27/compiled_odds_sample_tables_chromoplexy_2019_02_27.RDa")
ctcalls <- read.delim(file = "/srv/shared/vanloo/home/mtarabichi/PCAWG/chromothripsis/tableCT.Step6.022019.txt", as.is = T)
katresults <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_Kataegis_calls_JD_allcolumns.txt", as.is = T)
katresults_repsamples <- katresults[katresults$is_preferred & !is.na(katresults$signature), ]

source(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/kataegis_functions.R", local = T)
# source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/pcawg.colour.palette.R")

CLEANHISTOLOGYFILE <- "/srv/shared/vanloo/ICGC_annotations/summary_table_combined_annotations_v2.txt"

histology_all <- read_histology(histologyfile = CLEANHISTOLOGYFILE)
# 
# venndiagdf <- histology_all[histology_all$is_preferred, c("samplename", "histology_abbreviation")]
# venndiagdf$histology_abbreviation <- factor(x = venndiagdf$histology_abbreviation, levels = sort(unique(venndiagdf$histology_abbreviation)))
# venndiagdf$has_kat <- venndiagdf$samplename %in% katresults_repsamples$sample
# venndiagdf$has_cplexy <- venndiagdf$samplename %in% Chromo_odds$all_samples$sample
# venndiagdf$has_ct <- venndiagdf$samplename %in% ctcalls[ctcalls$FinalCalls == "Chromothripsis", "samplename"]
# 
# venndiagdf$inter <- interaction(venndiagdf[, c("has_kat", "has_cplexy", "has_ct")])
# venndiagdf$inter <- factor(x = venndiagdf$inter, levels = rev(c("TRUE.FALSE.FALSE", "TRUE.FALSE.TRUE", "TRUE.TRUE.TRUE", "TRUE.TRUE.FALSE", "FALSE.TRUE.FALSE", "FALSE.TRUE.TRUE", "FALSE.FALSE.TRUE", "FALSE.FALSE.FALSE")))
# 
# samplecountsdf <- do.call(rbind, by(data = venndiagdf, INDICES = venndiagdf$histology_abbreviation, FUN = function(x) data.frame(npunct = sum(x$inter != "FALSE.FALSE.FALSE"), ntot = nrow(x))))
# samplecountsdf$histology_abbreviation <- rownames(samplecountsdf)


library(ggplot2)
# 
# p1 <- ggplot(data = venndiagdf) + geom_bar(mapping = aes(x = histology_abbreviation, fill = inter), position = position_fill())
# p1 <- p1 + scale_fill_manual(values = c("TRUE.FALSE.FALSE" = rgb(0.8941176, 0.1019608, 0.1098039),
#                                         "TRUE.FALSE.TRUE" = rgb(0.5549020,0.2980392,0.4156863),
#                                         "TRUE.TRUE.TRUE" = rgb(0.4705882,0.4274510,0.3738562),
#                                         "TRUE.TRUE.FALSE" = rgb(0.5980392,0.3941176,0.2000000),
#                                         "FALSE.TRUE.FALSE" = rgb(0.3019608,0.6862745,0.2901961),
#                                         "FALSE.TRUE.TRUE" = rgb(0.2588235,0.5901961,0.5058824),
#                                         "FALSE.FALSE.TRUE" = rgb(0.2156863,0.4941176,0.7215686),
#                                         "FALSE.FALSE.FALSE" = rgb(1,1,1,0)), guide = F)
# p1 <- p1 + theme_minimal() + theme(axis.text = element_blank(), panel.grid = element_blank(), axis.title = element_blank())
# p1 <- p1 + ylim(y = c(0, 1.8))
# p1 <- p1 + geom_segment(data = data.frame(xstart = 1:length(levels(venndiagdf$histology_abbreviation)) - .3, xend = .3 + 1:length(levels(venndiagdf$histology_abbreviation))), mapping = aes(x = xstart, xend = xend, y = 1.45, yend = 1.45))
# p1 <- p1 + geom_text(data = samplecountsdf, mapping = aes(x = histology_abbreviation, y = 1.65, label = npunct))
# p1 <- p1 + geom_text(data = samplecountsdf, mapping = aes(x = histology_abbreviation, y = 1.25, label = ntot))
# p1
# ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190218_timing_results_samplefractions.pdf", plot = p1, width = 12, height = 1)
# 


#### new plot to top odds

stackbardf <- histology_all[histology_all$is_preferred, c("samplename", "histology_abbreviation")]
stackbardf$histology_abbreviation <- factor(x = stackbardf$histology_abbreviation, levels = sort(unique(stackbardf$histology_abbreviation)))
stackbardf$has_kat <- stackbardf$samplename %in% katresults_repsamples$sample
stackbardf$has_cplexy <- stackbardf$samplename %in% Chromo_odds$all_samples$sample
stackbardf$has_ct <- stackbardf$samplename %in% ctcalls[ctcalls$FinalCalls == "Chromothripsis", "samplename"]


# fion here
create_plotdf <- function(singlettdf) {

row1 <- data.frame(histology_abbreviation = singlettdf$histology_abbreviation[1],
                   event = "kataegis",
                   xmin = 0, xmax = 0,
                   ymin = 0, ymax = sum(singlettdf$has_kat)/nrow(singlettdf))

row2 <- data.frame(histology_abbreviation = singlettdf$histology_abbreviation[1],
                   event = "chromoplexy",
                   xmin = 0, xmax = 0,
                   ymin = 0, ymax = sum(singlettdf$has_kat & singlettdf$has_cplexy)/nrow(singlettdf))
row3 <- data.frame(histology_abbreviation = singlettdf$histology_abbreviation[1],
                   event = "chromoplexy",
                   xmin = 0, xmax = 0,
                   ymin = row1$ymax, ymax = row1$ymax + sum(!singlettdf$has_kat & singlettdf$has_cplexy)/nrow(singlettdf))

row4 <- data.frame(histology_abbreviation = singlettdf$histology_abbreviation[1],
                   event = "chromothripsis",
                   xmin = 0, xmax = 0,
                   ymin = 0, ymax = sum(singlettdf$has_kat & singlettdf$has_cplexy & singlettdf$has_ct)/nrow(singlettdf))
row5 <- data.frame(histology_abbreviation = singlettdf$histology_abbreviation[1],
                   event = "chromothripsis",
                   xmin = 0, xmax = 0,
                   ymin = row2$ymax, ymax = row2$ymax + sum(singlettdf$has_kat & !singlettdf$has_cplexy & singlettdf$has_ct)/nrow(singlettdf))
row6 <- data.frame(histology_abbreviation = singlettdf$histology_abbreviation[1],
                   event = "chromothripsis",
                   xmin = 0, xmax = 0,
                   ymin = row1$ymax, ymax = row1$ymax + sum(!singlettdf$has_kat & singlettdf$has_cplexy & singlettdf$has_ct)/nrow(singlettdf))
row7 <- data.frame(histology_abbreviation = singlettdf$histology_abbreviation[1],
                   event = "chromothripsis",
                   xmin = 0, xmax = 0,
                   ymin = row3$ymax, ymax = row3$ymax + sum(!singlettdf$has_kat & !singlettdf$has_cplexy & singlettdf$has_ct)/nrow(singlettdf))

# row2 <- data.frame(histology_abbreviation = singlettdf$histology_abbreviation[1],
#                    event = "chromothripsis",
#                    xmin = 0, xmax = 0,
#                    ymin = 0, ymax = sum(singlettdf$has_kat & singlettdf$has_ct)/nrow(singlettdf))
# row3 <- data.frame(histology_abbreviation = singlettdf$histology_abbreviation[1],
#                    event = "chromothripsis",
#                    xmin = 0, xmax = 0,
#                    ymin = row1$ymax, ymax = row1$ymax + sum(!singlettdf$has_kat & singlettdf$has_ct)/nrow(singlettdf))
# 
# row4 <- data.frame(histology_abbreviation = singlettdf$histology_abbreviation[1],
#                    event = "chromoplexy",
#                    xmin = 0, xmax = 0,
#                    ymin = 0, ymax = sum(singlettdf$has_kat & singlettdf$has_ct & singlettdf$has_cplexy)/nrow(singlettdf))
# row5 <- data.frame(histology_abbreviation = singlettdf$histology_abbreviation[1],
#                    event = "chromoplexy",
#                    xmin = 0, xmax = 0,
#                    ymin = row2$ymax, ymax = row2$ymax + sum(singlettdf$has_kat & !singlettdf$has_ct & singlettdf$has_cplexy)/nrow(singlettdf))
# row6 <- data.frame(histology_abbreviation = singlettdf$histology_abbreviation[1],
#                    event = "chromoplexy",
#                    xmin = 0, xmax = 0,
#                    ymin = row1$ymax, ymax = row1$ymax + sum(!singlettdf$has_kat & singlettdf$has_ct & singlettdf$has_cplexy)/nrow(singlettdf))
# row7 <- data.frame(histology_abbreviation = singlettdf$histology_abbreviation[1],
#                    event = "chromoplexy",
#                    xmin = 0, xmax = 0,
#                    ymin = row3$ymax, ymax = row3$ymax + sum(!singlettdf$has_kat & !singlettdf$has_ct & singlettdf$has_cplexy)/nrow(singlettdf))
outdf <- rbind(row1, row2, row3, row4, row5, row6, row7)
return(outdf)
}

plotdf <- do.call(rbind, lapply(X = split(x = stackbardf, f = stackbardf$histology_abbreviation), FUN = create_plotdf))
plotdf$histology_abbreviation <- factor(x = plotdf$histology_abbreviation)
plotdf$event <- factor(x = plotdf$event, levels = c("kataegis", "chromoplexy", "chromothripsis"))
plotdf$xmin <- as.numeric(x = plotdf$histology_abbreviation) + as.numeric(x = plotdf$event)/4 -.125
plotdf$xmax <- as.numeric(x = plotdf$histology_abbreviation) + as.numeric(x = plotdf$event)/4 +.125

plotdf2 <- do.call(rbind, by(data = plotdf, plotdf$histology_abbreviation, FUN = function(x) data.frame(histology_abbreviation = x$histology_abbreviation[1],
                                                                                         xmin = min(x$xmin), xmax = max(x$xmax), ymin = 0, ymax = 1)))
plotdf <- plotdf[plotdf$ymax > plotdf$ymin, ]

samplecountsdf <- do.call(rbind, by(data = stackbardf, INDICES = stackbardf$histology_abbreviation, FUN = function(x) data.frame(npunct = sum(x$has_kat | x$has_ct | x$has_cplexy), ntot = nrow(x))))
samplecountsdf$xpos <- as.numeric(factor(x = rownames(samplecountsdf)))

p1 <- ggplot(data = plotdf) + geom_rect(data = plotdf2, mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey90")
p1 <- p1 + geom_rect(mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = event), show.legend = F) + theme_minimal()
# p1 <- p1 + geom_segment(data = plotdf2, mapping = aes(x = xmin + .125, xend = xmax - .125 , y = 1.45, yend = 1.45))
# p1 <- p1 + geom_text(data = samplecountsdf, mapping = aes(x = xpos+.5, y = 1.65, label = npunct))
p1 <- p1 + geom_text(data = samplecountsdf, mapping = aes(x = xpos+.5, y = 1.15, label = ntot))
p1 <- p1 + theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank()) + 
  scale_fill_manual(values = c(kataegis = '#e41a1c', chromothripsis = '#377eb8', chromoplexy = '#4daf4a')) + ylim(c(0,1.2))
p1

ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190218_timing_results_samplefractions_NEW.pdf", plot = p1, width = 12, height = 1)



### plotting final odds ratios CT, CP, KAT
katodds_cvs <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_background_timing_odds_clonalVsubclonal_intnorm_cnstrat.txt", as.is = T)
katodds_cvs$histology <- rownames(katodds_cvs)
katodds_cvs$type <- "kataegis"
katodds_cvs$comp <- "clonal/subclonal"

katodds_evl <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_background_timing_odds_earlyVlate_intnorm_cnstrat.txt", as.is = T)
katodds_evl$histology <- rownames(katodds_evl)
katodds_evl$type <- "kataegis"
katodds_evl$comp <- "early/late"

load("/srv/shared/vanloo/home/mfittall/ICGC/Chromoplexy/Final_2019_2_27/compiled_odds_sample_tables_chromoplexy_2019_02_19.RDa")
cpodds_cvs <- as.data.frame(Chromo_odds$Clonal_subclonal_odds)
cpodds_cvs$histology <- rownames(cpodds_cvs)
colnames(cpodds_cvs) <- c("median", "lower", "upper", "histology")
cpodds_cvs$type <- "chromoplexy"
cpodds_cvs$comp <- "clonal/subclonal"

cpodds_evl <- as.data.frame(Chromo_odds$EvL_odds)
cpodds_evl$histology <- rownames(cpodds_evl)
colnames(cpodds_evl) <- c("median", "lower", "upper", "histology")
cpodds_evl$type <- "chromoplexy"
cpodds_evl$comp <- "early/late"

load(file = "/srv/shared/vanloo/home/mtarabichi/PCAWG/chromothripsis/oddsCS.nonas.22022019.alloldremoved.Rda")
ctodds_cvs <- cbind(as.data.frame(do.call(rbind, oddsCS)), names(oddsCS), stringsAsFactors = F)
colnames(ctodds_cvs) <- c("lower", "median", "upper", "histology")
ctodds_cvs$type <- "chromothripsis"
ctodds_cvs$comp <- "clonal/subclonal"

load(file = "/srv/shared/vanloo/home/mtarabichi/PCAWG/chromothripsis/oddsEL.nonasv2.22022019.alloldremoved.Rda")
ctodds_evl <- cbind(as.data.frame(do.call(rbind, oddsEL)), names(oddsEL), stringsAsFactors = F)
colnames(ctodds_evl) <- c("lower", "median", "upper", "histology")
ctodds_evl$type <- "chromothripsis"
ctodds_evl$comp <- "early/late"


oddsplotdf <- rbind(katodds_cvs, katodds_evl, cpodds_cvs, cpodds_evl, ctodds_cvs, ctodds_evl)
oddsplotdf$is_signif <- oddsplotdf$upper < 1 | oddsplotdf$lower > 1
oddsplotdf <- oddsplotdf[!is.na(oddsplotdf$median), ]
rownames(oddsplotdf) <- NULL
oddsplotdf[oddsplotdf$histology == "Kidney-RCC.clearcell", "histology"] <- "Kidney-RCC-Clear"
oddsplotdf[oddsplotdf$histology == "Kidney-RCC.papillary", "histology"] <- "Kidney-RCC-Pap"
oddsplotdf[oddsplotdf$histology == "Skin-Melanoma.acral", "histology"] <- "Skin-Melanoma-Acral"
oddsplotdf[oddsplotdf$histology == "Skin-Melanoma.cutaneous", "histology"] <- "Skin-Melanoma-Cut"
oddsplotdf[oddsplotdf$histology == "Skin-Melanoma.mucosal", "histology"] <- "Skin-Melanoma-Mucosal"

oddsplotdf$histology <- factor(x = oddsplotdf$histology, levels = sort(unique(histology_all$histology_abbreviation)))
oddsplotdf$type <- factor(oddsplotdf$type, levels = c("kataegis", "chromoplexy", "chromothripsis"))
# oddsplotdf$comp <- factor(oddsplotdf$comp, levels = c("early/late", "clonal/subclonal"))

samplesizes <- as.data.frame(do.call(rbind, by(data = venndiagdf[, c("has_kat", "has_ct", "has_cplexy")], INDICES = venndiagdf$histology_abbreviation, FUN = colSums)))
katexclude <- rownames(samplesizes)[samplesizes$has_kat < 3]
ctexclude <- rownames(samplesizes)[samplesizes$has_ct < 3]
cpexclude <- rownames(samplesizes)[samplesizes$has_cplexy < 3]

oddsplotdf <- oddsplotdf[!((oddsplotdf$histology %in% katexclude & oddsplotdf$type == "kataegis") |
                           (oddsplotdf$histology %in% ctexclude & oddsplotdf$type == "chromothripsis") |
                           (oddsplotdf$histology %in% cpexclude & oddsplotdf$type == "chromoplexy")), ]


p1 <- ggplot(data = oddsplotdf)
p1 <- p1 + geom_hline(yintercept = 1, colour = "grey", alpha = .5)
# p1 <- p1 + geom_pointrange(mapping = aes(y = median, ymin = lower, ymax = upper, x = as.numeric(histology) + (as.numeric(type)-2)/4, colour = type, alpha = is_signif), size = .3)
p1 <- p1 + geom_pointrange(mapping = aes(y = median, ymin = lower, ymax = upper, x = histology, colour = type, alpha = is_signif), position = position_dodge(width = .75), size = .3)
p1 <- p1 + geom_vline(xintercept = 1.5:(length(levels(oddsplotdf$histology)) - .5), colour = "grey", alpha = .5)
p1 <- p1 + scale_alpha_manual(values = c('TRUE' = .85, 'FALSE' = .15))
p1 <- p1 + scale_color_manual(values = c(kataegis = rgb(0.8941176, 0.1019608, 0.1098039), chromoplexy = rgb(0.3019608,0.6862745,0.2901961), chromothripsis = rgb(0.2156863,0.4941176,0.7215686)))
p1 <- p1 + scale_y_log10(breaks = c(0.01,0.1,1,10,100), labels = c(0.01,0.1,1,10,100)) + scale_x_discrete(drop = F)
p1 <- p1 + theme_minimal() + theme(axis.text.x = element_text(angle = 90), panel.grid.major.x = element_blank(), legend.position = "none")
p1 <- p1 + facet_wrap(~comp, nrow = 2)
p1 <- p1 + labs(x = "", y = "relative odds ratio")
p1
ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_greaterthan3_timing_results.pdf", plot = p1, width = 12, height = 4.5)


