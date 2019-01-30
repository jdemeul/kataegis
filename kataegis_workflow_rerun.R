## statics
library(doParallel)

## load functions
source(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/kataegis_functions.R")


# input
# SAMPLE <- "f221cbb5-eefa-187f-e040-11ac0c481708"
# SAMPLE <- "097a7d36-905b-72be-e050-11ac0d482c9a"
PCFDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20180309_pcf_rerun/"
SVDIR <- "/srv/shared/vanloo/ICGC-structural-variants/"
SNVDIR <- "/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/"
CNBREAKSDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/consensus_breakpoints/consensus_bp.basic.20161220"
# RHOPSIFILE <- "/srv/data/vanloo/jdemeul/ICGC/kataegis/clustering_consensus/1_purity_ploidy/purity_ploidy.txt"
RHOPSIFILE <- "/srv/shared/vanloo/ICGC_consensus_copynumber/consensus.20170217.purity.ploidy.txt.gz"
# MCNDIR <- "/srv/data/vanloo/jdemeul/ICGC/kataegis/clustering_consensus/0_multiplicity/"
# CLUSTDIR <- "/srv/data/vanloo/jdemeul/ICGC/kataegis/clustering_consensus/2_subclones/"
# CCLUSTDIR <- "/srv/shared/vanloo/ICGC-consensus-clustering/20180319_consensus_subclonal_reconstruction_beta1.5_svfix"
CCLUSTDIR <- "/srv/shared/vanloo/ICGC-consensus-clustering/consensus_subclonal_reconstruction_v1.1_20181121/"
# CTIMINGDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/20180310_consensus_subclonal_reconstruction_beta1.5_svfix_formattingpilot"
# CCCFDIR <- "/srv/shared/vanloo/ICGC-consensus-clustering/WM_release20170325/mutcff/"
CCCFDIR <- "/srv/shared/vanloo/ICGC-consensus-clustering/consensus_subclonal_reconstruction_v1.1_20181121_mutccf/"
CTIMINGDIR <- "/srv/shared/vanloo/ICGC-consensus-clustering/consensus_subclonal_reconstruction_v1.1_20181121_probgain/"
PHASING <- "/srv/shared/vanloo/ICGC_snv_mnv/phasing_final_consensus_12oct/"
# HISTOLOGYFILE <- "/srv/data/vanloo/khaase/ICGC/annotation/pcawg_specimen_histology_August2016_v6.fixed.tsv"
# RELEASEDATAFILE <- "/srv/data/vanloo/khaase/ICGC/annotation/release_may2016.v1.4.tsv"
WGTRINUCCONTFILE <-  "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/wg_trinuc_freq.txt"
CLEANHISTOLOGYFILE <- "/srv/shared/vanloo/ICGC_annotations/summary_table_combined_annotations_v2.txt"
# CLEANHISTOLOGYFILE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/merged_histology_v6_release_v1.4.txt"
CIRCOSDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/circos"
# RESULTSBASE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/analysed_rerun_annotmuts"
RESULTSBASE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/final_rerun_annotmuts"
# SIGPERSAMPLEFILE <- "/srv/shared/vanloo/ICGC_signatures/PCAWG_signatures_in_samples.20170302.txt.gz"
SIGFILE <- "/srv/shared/vanloo/ICGC_signatures/20180322_release/sigProfiler_SBS_signatures.csv"
KATSIGFILE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20180322_kataegis_signature_patterns.csv"
SIGSINSAMPLESFILE <- "/srv/shared/vanloo/ICGC_signatures/20180322_release/PCAWG_sigProfiler_SBS_signatures_in_samples_waliqID.csv"
RELEVANT_SIGS_ONLY <- T

# kataegis_threshold <- 1e3


NTHREADS <- 18
KATSAMPLES <- sub(pattern = "_kataegis_annotated.txt", replacement = "", x = list.files(path = PCFDIR))


## Run once
# write.table(x=data.frame(as.list(get_wg_trinuc_normalisation_factors())), file=WGTRINUCCONTFILE, sep = "\t", quote = F)
# write.table(x=get_clean_histology(histologyfile = HISTOLOGYFILE, releasedatafile = RELEASEDATAFILE), file=CLEANHISTOLOGYFILE, sep = "\t", quote = F)


## load one-time results
histology_all <- read_histology(histologyfile = CLEANHISTOLOGYFILE)

trinuc_freq <- unlist(read.delim(file = WGTRINUCCONTFILE)[1,])
pcawg_sigs_all <- load_signatures(SIGFILE, mergesigs = c(7,10,17))
kataegis_sigs <- load_signatures(signatures_file = KATSIGFILE, pad = T)
active_sigs_per_sample <- get_active_sigs_per_sample(sigs_per_sample_file = SIGSINSAMPLESFILE, mergesigs = c(7,10,17))

clp = makeCluster(NTHREADS)
registerDoParallel(clp)

foreach(i=1:length(KATSAMPLES), .packages = c("BSgenome", "BSgenome.Hsapiens.1000genomes.hs37d5", "ggplot2"), .verbose = T) %dopar% {
# foreach(i=1:2, .verbose = T) %dopar% {
# for (i in 1:length(KATSAMPLES)) {
  # i <- 3
  # source(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/kataegis_functions.R")
  
  SAMPLE <- KATSAMPLES[i]
  # print(SAMPLE)
  
  # output
  OUTDIR <- file.path(RESULTSBASE, SAMPLE)
  dir.create(path = file.path(OUTDIR), showWarnings = F)
  
  histology <- histology_all[histology_all$samplename == SAMPLE, "histology_abbreviation"]
  if (length(histology) == 0)
    histology <- NA
  
  ## load data
  alldata <- tryCatch( 
    {
      read_all_data(sample_id = SAMPLE, pcfdir = PCFDIR, svdir = SVDIR, snv_mnvdir = SNVDIR, cnbreaksdir = CNBREAKSDIR, 
                    rhopsifile = RHOPSIFILE, cclustdir = CCLUSTDIR, cccfdir = CCCFDIR, ctimingdir = CTIMINGDIR, phasingdir = PHASING)
    }, error = function(err) {
      print(paste("ERROR: sample ", SAMPLE, err))
      return(NULL)
    }, finally = {}
  )
  
  if (is.null(alldata)) next
  
  
  if (RELEVANT_SIGS_ONLY) {
    active_sigs <- active_sigs_per_sample[[SAMPLE]]
    if (is.null(active_sigs)) {
      samehistosamples <- histology_all[histology_all$histology_abbreviation == histology, "samplename"]
      active_sigs <- unique(unlist(active_sigs_per_sample[samehistosamples]))
    }
    if (grepl(pattern = "Lymph", x = histology) & !"SBS9" %in% active_sigs) {
      active_sigs <- c(active_sigs, "SBS9")
    }
    # focal_sigs <- paste0("Signature.", c("2", "13", "17a", "17b", "19", "28"))
    # pcawg_sigs <- pcawg_sigs_all[union(active_sigs, focal_sigs)]
    pcawg_sigs <- c(pcawg_sigs_all[active_sigs], kataegis_sigs)
  } else {
    pcawg_sigs <- pcawg_sigs_all
  }
  
  ## reidentify/subdivide pcf segments
  # alldata[["cpcf"]]$foci <- get_foci(mutations = alldata[["cpcf"]], max_distance = 2*kataegis_threshold)
  
  ## plotting
  if (nrow(alldata[["snv"]]) > 0) {
    plot_rainfall(mutations = alldata[["snv"]], sample = SAMPLE, outdir = OUTDIR, svbreaks = alldata[["sv"]],
                  plot_position = T, bsgenome = genome)
  }
  
  write.table(x = alldata[["cpcf"]], file = file.path(OUTDIR, paste0(SAMPLE, "_all_muts.txt")), sep = "\t", col.names = T, row.names = F, quote = F)
  # alldata[["sum_cpcf_all"]] <- do.call(rbind, by(data = alldata[["cpcf"]], INDICES = alldata[["cpcf"]]$foci,
  #                                                         FUN = summarize_segment_full, sample = SAMPLE, histology = histology, katobj = alldata, #subclones = alldata[["subcl"]]$ccf, svs = alldata[["sv"]], cnbreaks = alldata[["cn"]],
  #                                                         trinuc_freq = trinuc_freq, outdir = OUTDIR, simplify = T, plotting = F))
  alldata[["sum_cpcf"]] <- do.call(rbind, by(data = alldata[["cpcf"]], INDICES = alldata[["cpcf"]]$foci,
                                                 FUN = summarize_segment_full, sample = SAMPLE, histology = histology, katobj = alldata, #subclones = alldata[["subcl"]]$ccf, svs = alldata[["sv"]], cnbreaks = alldata[["cn"]],
                                                 trinuc_freq = trinuc_freq, outdir = OUTDIR, plotting = F, simplify = T))
  # alldata[["sum_cpcf"]] <- do.call(rbind, sapply(X = alldata$sum_cpcf_all, FUN = '[[', 1, simplify = F))
  # alldata[["sum_cpcf_muts"]] <- do.call(rbind, sapply(X = alldata$sum_cpcf_all, FUN = '[[', 2, simplify = F))
  write.table(x = unique(alldata[["sum_cpcf"]][, -c(1:4)]), file = file.path(OUTDIR, paste0(SAMPLE, "_kataegis_cpcf.txt")), sep = "\t", col.names = T, row.names = F, quote = F)
  write.table(x = alldata[["sum_cpcf"]][, c(1:9)], file = file.path(OUTDIR, paste0(SAMPLE, "_all_muts_annot.txt")), sep = "\t", col.names = T, row.names = F, quote = F)
  
}

stopCluster(clp)
