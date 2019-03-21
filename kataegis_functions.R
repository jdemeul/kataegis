## libraries
library(ggplot2)
library(BSgenome)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
# library(penalized)
# library(igraph)
genome <- BSgenome.Hsapiens.1000genomes.hs37d5

## flags up runs of non-standard SNVs
#
# find_alternative_kataegis <- function(mutations,
#                                       min_no_mutations = 5,
#                                       max_distance = 1e5,
#                                       exclude = c("CT", "GA", "CG", "GC", "CA", "GT")) {
#   
#   is_same_chrom <- c(T, mutations[-nrow(mutations), "chromosome"] == mutations[-1, "chromosome"])
#   is_same_stretch <- c(T, abs(diff(mutations$pos)) <= max_distance)
#   focus_start <- c(1, which(!is_same_chrom | !is_same_stretch))
#   focus_end <- c(which(!is_same_chrom | !is_same_stretch) - 1, nrow(mutations))
#   focus_length <- diff(c(1, which(!is_same_chrom | !is_same_stretch), nrow(mutations) + 1))
#   focus_factor <- rep(1:length(focus_length), focus_length)
#   return(list(focus_start, focus_end, focus_length, focus_factor))
# }

## get the different foci
get_foci <- function(mutations,
                     max_distance = 2.5e3) {
  if (length(mutations) == 0) 
    return(vector(mode = "integer"))
  is_same_chrom <- c(T, mutations[-nrow(mutations), "chromosome"] == mutations[-1, "chromosome"])
  is_same_stretch <- c(T, abs(diff(mutations$pos)) <= max_distance)
  focus_length <- diff(c(1, which(!is_same_chrom | !is_same_stretch), nrow(mutations) + 1))
  foci <- rep(1:length(focus_length), focus_length)
  return(foci)
}

## get the different foci
get_foci_stranded <- function(mutations,
                     max_distance = 2.5e3) {
  if (length(mutations) == 0) 
    return(vector(mode = "integer"))
  row.names(mutations) <- NULL
  mut_sorted <- mutations[order(mutations$ref, mutations$chromosome, mutations$pos), ]
  is_same_chrom <- c(T, mut_sorted[-nrow(mut_sorted), "chromosome"] == mut_sorted[-1, "chromosome"])
  is_same_stretch <- c(T, abs(diff(mut_sorted$pos)) <= max_distance)
  is_same_strand <- c(T, mut_sorted[-nrow(mut_sorted), "ref"] == mut_sorted[-1, "ref"])
  focus_length <- diff(c(1, which(!is_same_chrom | !is_same_stretch | !is_same_strand), nrow(mut_sorted) + 1))
  foci <- rep(1:length(focus_length), focus_length)[order(as.integer(rownames(mut_sorted)))]
  return(foci)
}



split_focus_apobec <- function(focus, max_distance = 2e3) {
  focus_filt <- focus[focus$ref %in% c("C", "G"), ]
  focus_filt$foci <- paste0(focus_filt$foci, "_", get_foci_stranded(mutations = focus_filt, max_distance = max_distance))
  
  other_muts <- focus[focus$ref %in% c("A", "T"), ]
  if (nrow(other_muts) > 0) {
    other_muts$foci <- paste0(other_muts$foci, 
                              unlist(apply(X = other_muts, MARGIN = 1, 
                                           FUN = function(atmuts, cgmuts) cgmuts[which.min(abs(as.integer(atmuts[["pos"]]) - cgmuts$pos)), "refoci"], 
                                           cgmuts = focus_filt)))
  }
  foci_all <- rbind(focus_filt, other_muts)
  return(foci_all[order(foci_all$pos), ])
}



split_focus_apobec_nodist <- function(focus) {
  focus_filt <- focus[focus$ref %in% c("C", "G"), ]
  focus_filt$foci <- focus_filt$foci + ifelse(focus_filt$ref == "C", 1000, 2000)
  
  other_muts <- focus[focus$ref %in% c("A", "T"), ]
  if (nrow(other_muts) > 0) {
    other_muts$foci <- sapply(X = other_muts$pos, FUN = function(atmutpos, cgmuts) cgmuts[which.min(abs(cgmuts$pos - atmutpos)), "foci"], cgmuts = focus_filt)
    }
  foci_all <- rbind(focus_filt, other_muts)
  return(foci_all[order(foci_all$pos), ])
}




# 
# ## deprecated and replaced by filter_stand_nonstand
# find_alternative_patterns <- function(mutations,
#                                       foci,
#                                       min_no_mutations = 5,
#                                       standard = c("CT", "GA", "CG", "GC", "CA", "GT")) {
#   if (length(mutations) == 0) return(data.frame())
#   types <- by(data = paste0(mutations$ref, mutations$alt), INDICES = foci, FUN = table)
#   max_type_ocurrences <- sapply(types, max)
#   max_type <- sapply(types, function(x) names(which.max(x)))
#   of_interest <- names(which( (max_type_ocurrences >= min_no_mutations) & !(max_type %in% standard) ))
#   mutations_of_interest <- mutations[foci %in% of_interest, ]
#   return(mutations_of_interest)
# }
# 
# filter_stand_nonstand <- function(mutations,
#                                   foci,
#                                   min_no_mutations = 0,
#                                   standard = c("CT", "GA", "CG", "GC", "CA", "GT"),
#                                   return_standard = T) {
#   if (length(mutations) == 0) return(data.frame())
#   types <- by(data = paste0(mutations$ref, mutations$alt), INDICES = foci, FUN = table)
#   max_type_ocurrences <- sapply(types, max)
#   max_type <- sapply(types, function(x) names(which.max(x)))
#   if (return_standard) {
#     of_interest <- names(which( (max_type_ocurrences >= min_no_mutations) & (max_type %in% standard) ))
#   } else {
#     of_interest <- names(which( (max_type_ocurrences >= min_no_mutations) & !(max_type %in% standard) ))
#   }
#   mutations_of_interest <- mutations[foci %in% of_interest, ]
#   return(mutations_of_interest)
# }


# get the standard Sig2/13 APOBEC foci
# get_kataegis_apobec <- function(mutations,
#                                 min_no_mutations = 3,
#                                 standard = c("TC[ACGT]", "[ACGT]GA"),
#                                 max_distance = 1e4) {
#   #filter
#   has_apobec_sig <- grepl(pattern = paste0(standard[1], "|", standard[2]), x = mutations$trinuc, perl = T)
#   mutations_filt <- mutations[has_apobec_sig, ]
#   mutations_filt$foci <- get_foci(mutations = mutations_filt, max_distance = max_distance)
#   
#   informative_foci <- which(rle(mutations_filt$foci)$lengths >= min_no_mutations)
#   mutations_filt <- mutations_filt[mutations_filt$foci %in% informative_foci, ]
#   return(mutations_filt)
# }


get_kataegis_apobec_sens <- function(mutations,
                                min_no_mutations = 3,
                                standard = c("TC[ACGT]", "[ACGT]GA"),
                                max_distance = 2.5e3,
                                max_reinclusion_dist = 1e3,
                                min_frac_APOBEC_context = 0.5) {
  #filter
  mutations_filt$foci <- get_foci_stranded(mutations = mutations_filt, max_distance = max_distance)
  non_singlets <- rle(mutations_filt$foci)$values[rle(mutations_filt$foci)$lengths > 1]
  mutations_filt <- mutations_filt[mutations_filt$foci %in% non_singlets, ]
  if (nrow(mutations_filt) == 0) return(mutations_filt)
  
  other_muts <- mutations[!has_apobec_sig & mutations$chromosome %in% unique(mutations_filt$chromosome), ]
  if (nrow(other_muts) > 0) {
    mutations_filt <- readd_non_context(mutations_filt = mutations_filt, other_muts = other_muts, mutations_all = mutations,
                                       max_reinclusion_dist = max_reinclusion_dist, 
                                       min_frac_APOBEC_context = min_frac_APOBEC_context, standard = standard)
  }
  informative_foci <- rle(mutations_filt$foci)$values[rle(mutations_filt$foci)$lengths >= min_no_mutations]
  mutations_filt <- mutations_filt[mutations_filt$foci %in% informative_foci, ]
  return(mutations_filt)
}


get_kataegis <- function(mutations,
                                 min_no_mutations = 3,
                                 type,
                                 max_distance = 2.5e3,
                                 max_reinclusion_dist = 1e3,
                                 min_frac_context = 0.5) {
  if (type == "apobec")
    standard <- c("TC[ACGT]", "[ACGT]GA")
  else if (type == "s17")
    standard <- c("TT", "AA")
  else
    stop("not a known kataegis type")
  
  #filter
  has_sig <- grepl(pattern = paste0(standard[1], "|", standard[2]), x = mutations$trinuc, perl = T)
  mutations_filt <- mutations[has_sig, ]
  mutations_filt$foci <- get_foci_stranded(mutations = mutations_filt, max_distance = max_distance)
  # non_singlets <- rle(mutations_filt$foci)$values[rle(mutations_filt$foci)$lengths > 1]
  tab1 <- table(mutations$foci)
  non_singlets <- names(tab1[tab1 > 1])
  mutations_filt <- mutations_filt[mutations_filt$foci %in% non_singlets, ]
  if (nrow(mutations_filt) == 0) return(mutations_filt)
  
  other_muts <- mutations[!has_sig & mutations$chromosome %in% unique(mutations_filt$chromosome), ]
  if (nrow(other_muts) > 0) {
    mutations_filt <- readd_non_context(mutations_filt = mutations_filt, other_muts = other_muts, mutations_all = mutations,
                                        max_reinclusion_dist = max_reinclusion_dist, 
                                        min_frac_context = min_frac_context, standard = standard)
  }
  tab2 <- table(mutations$foci)
  informative_foci <- names(tab2[tab2 >= min_no_mutations])
  mutations_filt <- mutations_filt[mutations_filt$foci %in% informative_foci, ]
  return(mutations_filt)
}


get_kataegis_sb <- function(mutations,
                         min_no_mutations = 5,
                         max_distance = 2.5e3) {
  mutations$foci <- get_foci_stranded(mutations = mutations, max_distance = max_distance)
  tab <- table(mutations$foci)
  informative_foci <- names(tab[tab >= min_no_mutations])
  mutations <- mutations[mutations$foci %in% informative_foci, ]
  return(mutations)
}



# deprecated and replaced by more general version below
# readd_non_context <- function(mutations_filt, other_muts, mutations_all, max_reinclusion_dist = 1e3,
#                               min_frac_APOBEC_context = 0.5, standard = c("TC[ACGT]", "[ACGT]GA")) {
#   other_muts$foci <- apply(X = other_muts, MARGIN = 1, FUN = assign_to_focus_nearest_mut, in_context_foci = mutations_filt, max_reinclusion_dist = max_reinclusion_dist)
#   non_apobec_foci <- filter_excess_noncontext(mutations_filt = mutations_filt, other_muts = other_muts, min_frac_APOBEC_context = min_frac_APOBEC_context)
#   
#   muts_out <- rbind(mutations_filt, other_muts[!is.na(other_muts$foci) & other_muts$ref %in% c("C", "G"), ])
#   muts_out <- muts_out[!muts_out$foci %in% non_apobec_foci, ]
#   if (nrow(muts_out) > 0) {
#     muts_out <- do.call(rbind, by(data = muts_out, INDICES = muts_out$foci, FUN = filter_strand, standard = standard))
#     muts_out <- muts_out[order(muts_out$foci, muts_out$pos), ]
#   }
#   return(muts_out)
# }

readd_non_context <- function(mutations_filt, other_muts, mutations_all, max_reinclusion_dist = 1e3,
                              min_frac_context = 0.5, standard) {
  other_muts$foci <- apply(X = other_muts, MARGIN = 1, FUN = assign_to_focus_nearest_mut, in_context_foci = mutations_filt, max_reinclusion_dist = max_reinclusion_dist)
  noncontext_foci <- filter_excess_noncontext(mutations_filt = mutations_filt, other_muts = other_muts, min_frac_context = min_frac_context)
  
  if (any(grepl(pattern = "TT", x = standard)))
    same_strand_muts <- c("A", "T")
  else if (any(grepl(pattern = "TC", x = standard)))
    same_strand_muts <- c("C", "G")
  else
    same_strand_muts <- c("A", "T", "C", "G")
  muts_out <- rbind(mutations_filt, other_muts[!is.na(other_muts$foci) & other_muts$ref %in% same_strand_muts, ])
  muts_out <- muts_out[!muts_out$foci %in% noncontext_foci, ]
  if (nrow(muts_out) > 0) {
    muts_out <- do.call(rbind, by(data = muts_out, INDICES = muts_out$foci, FUN = filter_strand, standard = standard))
    muts_out <- muts_out[order(muts_out$foci, muts_out$pos), ]
  }
  return(muts_out)
}


filter_excess_noncontext <- function(mutations_filt, other_muts, min_frac_context = 0.5) {
  t1 <- table(other_muts[!is.na(other_muts$foci), "foci"])
  t2 <- table(mutations_filt[, "foci"])
  frac_in_context <- t2[names(t1)]/(t1+t2[names(t1)])
  noncontext_foci <- names(frac_in_context[frac_in_context <= min_frac_context])
  return(noncontext_foci)
}


filter_strand <- function(merged_focus, standard) {
  strand_counts <- c(plus = sum(grepl(pattern = standard[1], x = merged_focus$trinuc, perl = T)),
              minus = sum(grepl(pattern = standard[2], x = merged_focus$trinuc, perl = T)))
  if (which.max(strand_counts) == 1)
    bases <- c("C", "T")
  else 
    bases <- c("G", "A")
  filtered_focus <- merged_focus[merged_focus$ref %in% bases, ]
  return(filtered_focus)
}


assign_to_focus <- function(mutation, segments) {
  is_in_focus <- mutation[["chromosome"]] == segments$chr & mutation[["pos"]] >= segments$start & mutation[["pos"]] <= segments$end
  if (any(is_in_focus))
    return(segments[is_in_focus, "foci"][1])
  else
    return(NA)
}


assign_to_focus_nearest_mut <- function(mutation, in_context_foci, max_reinclusion_dist) {
  on_chrom_foci <- in_context_foci[in_context_foci$chromosome == mutation[["chromosome"]], ]
  if (nrow(on_chrom_foci) == 0)
    return(NA)
  nearest_mut <- on_chrom_foci[which.min(abs(as.integer(mutation[["pos"]]) - on_chrom_foci$pos)), ]
  if (abs(as.integer(mutation[["pos"]]) - nearest_mut$pos) > max_reinclusion_dist)
    return(NA)
  else 
    return(nearest_mut$foci)
}





# non-APOBEC kataegis sites = foci which do not contain mutations identified in the APOBEC foci
get_kataegis_other <- function(mutations,
                                   min_no_mutations = 5,
                                   typed_mutations) {
noninformative_foci <- which(rle(mutations$foci)$lengths < min_no_mutations)
typedmut_containing_foci <- unique(mutations[match(x = paste0(typed_mutations$chromosome, "_", typed_mutations$pos), 
                                                        table = paste0(mutations$chromosome, "_", mutations$pos)), "foci"])
nonstand <- mutations[!mutations$foci %in% union(typedmut_containing_foci, noninformative_foci), ]
return(nonstand)
}

# 
# filter_known_signal <- function(mutations,
#                                 signal = c("CT", "GA", "CG", "GC", "CA", "GT")) {
#   is_known <- paste0(mutations$ref, mutations$alt) %in% signal
#   return(mutations[is_known, ])
# }

# deprecated
# summarize a set of segments containing mutations
# summarize_segment <- function(focus, subclones, svs, cnbreaks, trinuc_freq = trinuc_freq, sample = SAMPLE, histology = histology) {
#   chrom <- unique(focus$chromosome)
#   # if (length(chrom) > 1) stop("focus stretches multiple chromosomes")
#   locus <- focus[c(1, nrow(focus)), "pos", drop = T]
#   types <- paste0(focus$ref, focus$alt)
#   typecounts <- c(sum(types %in% c("CA", "GT")), 
#                   sum(types %in% c("CG", "GC")),
#                   sum(types %in% c("CT", "GA")),
#                   sum(types %in% c("TA", "AT")), 
#                   sum(types %in% c("TC", "AG")),
#                   sum(types %in% c("TG", "AC")))
#   mutcount <- nrow(focus)
#   informative <- sum(!is.na(focus$subclonal.fraction))
#   in_context <- sum(grepl(pattern = "TC[ACGT]|[ACGT]GA", x = focus$trinuc, perl = T))
#   out_context <- mutcount - in_context
#   switchcount <- length(rle(focus$ref)$lengths)
#   if (all(is.na(focus$subclonal.fraction))) 
#     ccfinfo <- 1 
#   else 
#     ccfinfo <- mean(winsor(focus$subclonal.fraction), na.rm = T)
#   subclccf <- subclones[ which.min( abs(ccfinfo - subclones) ) ]
#   CN <- has_cnbreak(chr = chrom, startpos = locus[1], endpos = locus[2], cnbreaks = cnbreaks)
#   if (CN) 
#     cnsize <- NA
#   else
#     cnsize <- get_minimal_cnsegment_size(chrom = chrom, start = locus[1], end = locus[2], cnbreaks = cnbreaks)
#   SV <- get_dist_nearest_sv(chr = chrom, startpos = locus[1], endpos = locus[2], svbreaks = svs)
#   mutspectrum <- as.data.frame(t(get_focal_mutspectrum(mutations = focus, trinuc_freq = trinuc_freq)))
#   return(data.frame(sample = sample, histology = histology, chr = chrom, start = locus[1], end = locus[2],
#                     CT = typecounts[1], CG = typecounts[2], CA = typecounts[3], TG = typecounts[4], TC = typecounts[5], TA = typecounts[6],
#                     in_context = in_context, out_context = out_context, switches = switchcount,
#                     total = mutcount, informative = informative, katccf = ccfinfo, subclccf = subclccf,
#                     sv = SV, cn = CN, cnsize = cnsize, mutspectrum,
#                     stringsAsFactors = F))
# }


# # summarize a set of segments containing mutations
# summarize_segment <- function(focus, subclones, svs, cnbreaks, trinuc_freq = trinuc_freq, sample = SAMPLE, histology = histology) {
#   chrom <- unique(focus$chromosome)
#   # if (length(chrom) > 1) stop("focus stretches multiple chromosomes")
#   locus <- focus[c(1, nrow(focus)), "pos", drop = T]
#   types <- paste0(focus$ref, focus$alt)
#   # typecounts <- c(sum(types %in% c("CA", "GT")), 
#   #                 sum(types %in% c("CG", "GC")),
#   #                 sum(types %in% c("CT", "GA")),
#   #                 sum(types %in% c("TA", "AT")), 
#   #                 sum(types %in% c("TC", "AG")),
#   #                 sum(types %in% c("TG", "AC")))
#   mutcount <- nrow(focus)
#   informative <- sum(!is.na(focus$subclonal.fraction))
#   # in_context <- sum(grepl(pattern = "TC[ACGT]|[ACGT]GA", x = focus$trinuc, perl = T))
#   # out_context <- mutcount - in_context
#   switchcount <- length(rle(focus$ref)$lengths)
#   if (all(is.na(focus$subclonal.fraction))) 
#     ccfinfo <- 1
#   else 
#     ccfinfo <- mean(winsor(focus$subclonal.fraction), na.rm = T)
#   subclccf <- subclones[ which.min( abs(ccfinfo - subclones) ) ]
#   
#   nMaj <- Mode(focus$nMaj1)
#   nMin <- Mode(focus$nMin1)
#   no.chrs.bearing.focus <- Mode(focus$no.chrs.bearing.mut)
#   timing <- ifelse(subclccf >= .9 & subclccf <= 1.1,
#                    ifelse(nMaj >= 2,
#                           ifelse(no.chrs.bearing.focus >= 2, "clonal.early", 
#                                  ifelse(no.chrs.bearing.focus == 1 & nMin == 0, "clonal.late", "clonal.na")
#                           ),
#                           "clonal.na"),
#                    "subclonal")
#   
#   CN <- has_cnbreak(chr = chrom, startpos = locus[1], endpos = locus[2], cnbreaks = cnbreaks)
#   if (CN) 
#     cnsize <- NA
#   else
#     cnsize <- get_minimal_cnsegment_size(chrom = chrom, start = locus[1], end = locus[2], cnbreaks = cnbreaks)
#   SV <- get_dist_nearest_sv(chr = chrom, startpos = locus[1], endpos = locus[2], svbreaks = svs)
#   # mutspectrum <- as.data.frame(t(get_focal_mutspectrum(mutations = focus, trinuc_freq = trinuc_freq)))
#   return(data.frame(sample = sample, histology = histology, chr = chrom, start = locus[1], end = locus[2],
#                     # CT = typecounts[1], CG = typecounts[2], CA = typecounts[3], TG = typecounts[4], TC = typecounts[5], TA = typecounts[6],
#                     nMaj = nMaj, nMin = nMin, no.chrs.bearing.focus = no.chrs.bearing.focus, 
#                     # in_context = in_context, out_context = out_context,
#                     switches = switchcount, total = mutcount, informative = informative,
#                     katccf = ccfinfo, subclccf = subclccf, timing = timing,
#                     sv = SV, cn = CN, cnsize = cnsize,
#                     # mutspectrum,
#                     stringsAsFactors = F))
# }


# summarize a set of segments containing mutations
summarize_segment_full <- function(focus, katobj, trinuc_freq = trinuc_freq, sample = SAMPLE, histology = histology, signatures = pcawg_sigs, outdir = OUTDIR, plotting = F) {
  subclones <- katobj[["csubcl"]]$fraction_cancer_cells
  svs <- katobj[["sv"]]
  cnbreaks <- katobj[["cn"]]
  focus <- focus[order(focus$pos), ]
  
  # active_sigs <- get_focal_signature(focal_mutations = focus, signatures = signatures, trinuc_genome = trinuc_freq, n = 1)
  # active_sigs2 <- get_focal_signature(focal_mutations = focus, signatures = signatures, trinuc_genome = trinuc_freq, n = 2)
  # 
  sigs_out <- get_cosine_sig_sims(focal_mutations = focus, signatures = signatures, trinuc_genome = trinuc_freq)
  # sigs_order <- order(x = sigs_out$cos, decreasing = T)
  # sigs_order_multi <- order(x = sigs_out$dmulti, decreasing = T)
  sigs_rank1 <- integer(length = length(sigs_out$cos))
  sigs_rank2 <- integer(length = length(sigs_out$cos))
  sigs_rank1[order(sigs_out$cos, decreasing = T)] <- 1:length(sigs_out$cos)
  sigs_rank2[order(sigs_out$dmulti, decreasing = T)] <- 1:length(sigs_out$dmulti)
  sigs_rank <- setNames(object = sigs_rank1 + sigs_rank2, nm = names(sigs_out$cos))
  sigs_rank <- sort(x = sigs_rank, decreasing = F)[1:3]
  sigs_names <- names(sigs_rank)
  sigs_names[sigs_rank >= 8] <- NA
  
  sigs_cosines <- sort(sigs_out$cos, decreasing = T)
  sigs_multi <- sort(sigs_out$dmulti, decreasing = T)
    
  # coef_tstat <- active_sigs[["coef_tstat"]]
  
  stranded <- is_stranded(focal_mutations = focus, trinuc_genome = trinuc_freq)
  pstreak <- get_prob_streak(focal_mutations = focus, trinuc_genome = trinuc_freq)
  deaminated <- is_deaminated(focal_mutations = focus, trinuc_genome = trinuc_freq)
  aidprint <- has_AID_footprint(focal_mutations = focus, trinuc_genome = trinuc_freq)
  
  types <- table(factor(focus$ref, levels = c("A", "C", "G", "T")))
  switchcount <- length(rle(focus$ref)$lengths)
  
  kat_phased <- subset(x = alldata[["phasing"]], subset = Chr == unique(focus$chromosome) & Pos1 %in% focus$pos & Pos2 %in% focus$pos)
  if (nrow(kat_phased) == 0) {
    phasing_stats <- c(no_phased_muts = 0, no_subclonal_muts = 0, no_antiphased_muts = 0)
  } else { 
    phasing_stats <- get_phasing_graph_stats(kat_phased = kat_phased, focus = focus)
    if (plotting & sigs_names[1] %in% paste0("SBS", c(2, 13, 19, "17a", "17b", "APO", "CTT", "ALT")))
    # if (plotting & active_sigs[["signatures"]] %in% c("Signature.2", "Signature.13", "Signature.19", "Signature.17a", "Signature.17b"))
    # if (plotting == T & active_sigs[["signatures"]] %in% c("Signature.2", "Signature.13") & (phasing_stats["no_phased_muts"] >= 10 | phasing_stats["no_subclonal_muts"] > 2 | phasing_stats["no_antiphased_muts"] > 2))
        plot_phasing_graph(kat_phased = kat_phased, focus = focus, outdir = outdir, sample_id = sample)
  }
  
  if (all(na.omit(sigs_names) %in% paste0("SBS", c(2, 13, "APO"))) & types[["C"]] >= 4 & types[["G"]] >= 4 & types[["A"]] <= 2 & types[["T"]] <= 2) {
  # if (all(c("Signature.2", "Signature.13") %in% active_sigs2[["signatures"]]) & types[["C"]] >= 4 & types[["G"]] >= 4 & types[["A"]] <= 2 & types[["T"]] <= 2 & active_sigs2[["varexpl"]] > 25) {
  # if (all(c("Signature.2", "Signature.13") %in% active_sigs2[["signatures"]]) & types[["C"]] >= 4 & types[["G"]] >= 4 & switchcount <= 8 & active_sigs2[["varexpl"]] >= 25) {
      
    splitfocus <- split_focus_apobec_nodist(focus)
    splitfoci_sumarized <- do.call(rbind, by(data = splitfocus, INDICES = splitfocus$foci,
                      FUN = summarize_segment_full, sample = sample, histology = histology, katobj = katobj, #subclones = subclones, svs = svs, cnbreaks = cnbreaks,
                      trinuc_freq = trinuc_freq, signatures = pcawg_sigs, outdir = OUTDIR, plotting = F, simplify = T))

    ## sanity check
    if (all(splitfoci_sumarized$sig1 %in% paste0("SBS", c(2, 13, "APO"))))
      return(splitfoci_sumarized)
  }
  # typecounts <- c(sum(types %in% c("CA", "GT")), 
  #                 sum(types %in% c("CG", "GC")),
  #                 sum(types %in% c("CT", "GA")),
  #                 sum(types %in% c("TA", "AT")), 
  #                 sum(types %in% c("TC", "AG")),
  #                 sum(types %in% c("TG", "AC")))
  
  
  chrom <- unique(focus$chromosome)
  locus <- focus[c(1, nrow(focus)), "pos", drop = T]
  types <- paste0(focus$ref, focus$alt)
  mutcount <- nrow(focus)
  informative <- sum(!is.na(focus$ccf))
  # informative <- sum(!is.na(focus$subclonal.fraction))
  if (all(is.na(focus$ccf))) 
    ccfinfo <- 1
  else 
    ccfinfo <- mean(winsor(focus$ccf), na.rm = T)
  subclccf <- subclones[ which.min( abs(ccfinfo - subclones) ) ]
  # subclccf <- katobj$csubcl[katobj$csubcl$cluster == sub(pattern = "cluster_", replacement = "", x = names(which.max(apply(X = focus[ , grep(colnames(focus), pattern = "cluster_"), drop = F], MARGIN = 2, FUN = function(x) sum(log(x), na.rm = T))))), "ccf"]

  # prior_clonal <- katobj[["csubcl"]][1, "n_snvs"] / sum(katobj[["csubcl"]][, "n_snvs"])
  posterior <- c(p_clonal = exp(sum(log(focus$cluster_1), na.rm = T)), p_subclonal = exp(sum(log(1-focus$cluster_1), na.rm = T)) )
  # posterior <- c(p_clonal = exp(sum(log(focus$cluster_1), log(prior_clonal))), p_subclonal = exp(sum(log(1 - focus$cluster_1), log(1 - prior_clonal))) )
  posterior <- posterior/sum(posterior)
    
  nMaj <- Mode(focus$major_cn)
  nMin <- Mode(focus$minor_cn)
  
  no.chrs.bearing.focus <- Mode(focus$mult)
  # timing <- ifelse(subclccf >= .9,
  #                  ifelse(nMaj >= 2,
  #                         ifelse(no.chrs.bearing.focus >= 2, "clonal.early", 
  #                                ifelse(no.chrs.bearing.focus == 1 & nMin == 0, "clonal.late", "clonal.na")
  #                         ),
  #                         "clonal.na"),
  #                  "subclonal")
  
  if (informative != 0) {
    # timing <- names(which.max(table(factor(focus$timing, levels = c("clonal [NA]", "clonal [early]", "clonal [late]", "subclonal")))))
    # timing_if_clonal <- names(which.max(table(factor(focus$timing, levels = c("clonal [NA]", "clonal [early]", "clonal [late]")))))
    
    
    #### edits 20180110 for probabilistic assignment of clonal early/late/NA classes
    # classes_if_clonal <- ifelse(focus$timing == "subclonal", ifelse(focus$major_cn >= 2 & focus$minor_cn == 0, ifelse(focus$mult >= 2, "clonal [early]", "clonal [late]"), "clonal [NA]"), focus$timing)
    # classes_if_clonal <- ifelse(focus$timing != "subclonal", focus$timing,
    #                             ifelse(focus$major_cn < 2, "clonal [NA]",
    #                                    ifelse(focus$mult >= 2, "clonal [early]",
    #                                           ifelse(focus$minor_cn == 0, "clonal [late]", "clonal [NA]"))))
    # classes_if_clonal <- factor(x = classes_if_clonal, levels = c("clonal [NA]", "clonal [early]", "clonal [late]"))
    # w_classes_if_clonal <- c(by(data = focus$cluster_1, INDICES = classes_if_clonal, FUN = sum, na.rm = T)) / sum(focus$cluster_1, na.rm = T)
    # w_classes_if_clonal[is.na(w_classes_if_clonal)] <- 0

    #### edits 20180314 for probabilistic assignment of clonal early/late/NA classes using MutationTimer data
    informative_focus <- focus[!(is.na(focus$prob_clonal_early) | is.na(focus$prob_clonal_late) | is.na(focus$prob_subclonal)), ]
    siminfo <- table(factor(x = ifelse(informative_focus$major_cn < 2, "nogain", ifelse(informative_focus$minor_cn == 0, "gainloh", "gain")), levels = c("nogain", "gain", "gainloh")))
    probs_earlylate <- get_probs_earlylate(informative_focus)
    
    # timing <- names(which.max(table(factor(focus$timing, levels = c("clonal [NA]", "clonal [early]", "clonal [late]", "subclonal")))))
    # timing_focus <- c(table(factor(focus$timing, levels = c("clonal [NA]", "clonal [early]", "clonal [late]", "subclonal"))))
    # timing_chisq_p <- chisq.test(x = timing_focus, p = katobj[["timing_summary"]], rescale.p = T, simulate.p.value = T)$p.value
    
    # timing_binom_p <- pbinom(q = sum(!focus$timing %in% timing), size = nrow(focus),
                             # prob = 1 - katobj[["timing_summary"]][timing] / sum(katobj[["timing_summary"]]), lower.tail = T)
    
    # subcl_chisq_p <- chisq.test(x = colSums(focus[ , grep(colnames(focus), pattern = "cluster_")], na.rm = T),
                                # p = katobj[["csubcl"]]$n_snvs, rescale.p = T, simulate.p.value = T, B = 10000)$p.value
    # if (length(grep(colnames(focus), pattern = "cluster_")) > 1) {
      # subcl_chisq_p <- suppressWarnings(chisq.test(x = colSums(focus[ , grep(colnames(focus), pattern = "cluster_"), drop = F], na.rm = T), 
                                                   # p = katobj[["csubcl"]]$n_snvs, rescale.p = T)$p.value)
      # subcl_binom_p <- pbinom(q = round(sum(focus$cluster_1)), size = nrow(focus),
      # prob = prior_clonal, lower.tail = T)
    # } else {
      # subcl_chisq_p <- NA
    # }
  } else {
    # timing <- NA
    # timing_if_clonal <- NA
    # w_classes_if_clonal <- rep(x = NA, 3)
    siminfo <- list(nogain = 0, gain = 0, gainloh = 0)
    probs_earlylate <- c(p_early = NA, p_late = NA, p_sub = NA, p_na = NA)
    # timing_focus <- NA
    # timing_chisq_p <- NA
    # timing_binom_p <- NA
    # subcl_chisq_p <- NA
    # subcl_chisq_p_nosim <- NA
    # subcl_binom_p <- NA
  }
  
  if (nrow(katobj[["snv"]]) > 0) {
    focal_prob <- get_focal_prob(all_mutations = katobj[["snv"]], focal_mutations = focus)
  } else {
    focal_prob <- NA 
  }

  
  CN <- has_cnbreak(chr = chrom, startpos = locus[1], endpos = locus[2], cnbreaks = cnbreaks)
  if (CN) 
    cnsize <- NA
  else
    cnsize <- get_minimal_cnsegment_size(chrom = chrom, start = locus[1], end = locus[2], cnbreaks = cnbreaks)
  if (is.null(svs)) {
    SV <- data.frame(sv_dist = NA, sv_pos = NA)
  } else if (nrow(svs) > 0) {
    SV <- get_dist_nearest_sv(chr = chrom, startpos = locus[1], endpos = locus[2], svbreaks = svs, dist_only = F)
  } else {
    SV <- data.frame(sv_dist = Inf, sv_pos = NA)
  }
  
  mutspectrum <- paste0(get_focal_mutspectrum(mutations = focus, trinuc_freq = trinuc_freq), collapse = ",")
  
  outdf <- data.frame(sample = sample, histology = histology, chr = chrom, start = locus[1], end = locus[2],
                   nMaj = nMaj, nMin = nMin, nNoGain = siminfo[["nogain"]], nGain = siminfo[["gain"]], nGainLOH = siminfo[["gainloh"]], no.chrs.bearing.focus = no.chrs.bearing.focus, 
                   switches = switchcount, total = mutcount, informative = informative,
                   katccf = ccfinfo, subclccf = subclccf,
                   # timing = timing, timing_if_clonal,
                   # w_clonal_NA = w_classes_if_clonal[1], w_clonal_early = w_classes_if_clonal[2], w_clonal_late = w_classes_if_clonal[3],
                   # timing_chisq_p = timing_chisq_p,
                   # timing_binom_p = timing_binom_p,
                   # subcl_chisq_p = subcl_chisq_p, #subcl_binom_p = subcl_binom_p,
                   p_clonal = posterior[1], p_subclonal = posterior[2],
                   p_early = probs_earlylate["p_early"], p_late = probs_earlylate["p_late"], p_sub = probs_earlylate["p_sub"], p_na = probs_earlylate["p_na"],
                   SV, cn = CN, cnsize = cnsize, # active_sig = active_sigs[["signatures"]],
                   # var_expl = active_sigs[["varexpl"]], active_sigs2 = paste(active_sigs2[["signatures"]], collapse = ","),
                   # var_expl2 = active_sigs2[["varexpl"]],
                   # sig_cos1 = sigs_cosines[1], sig_cos1_name = names(sigs_cosines[1]), sig_cos2 = sigs_cosines[2], sig_cos2_name = names(sigs_cosines[2]),
                   sig_cosines = paste0(sigs_cosines, collapse = ","), sig_cosines_names = paste0(sub(pattern = "SBS", replacement = "", x = names(sigs_cosines)), collapse = ","),
                   # sig_mult1 = sigs_multi[1], sig_mult1_name = names(sigs_multi[1]), sig_mult2 = sigs_multi[2], sig_mult2_name = names(sigs_multi[2]),
                   sig_multi = paste0(sigs_multi, collapse = ","), sig_multi_names = paste0(sub(pattern = "SBS", replacement = "", x = names(sigs_multi)), collapse = ","),
                   sig1 = sigs_names[1], sig2 = sigs_names[2], sig3 = sigs_names[3],
                   p_nostrandbias = stranded, p_streak = pstreak, p_deaminated = deaminated, p_aid = aidprint,
                   no_phased_muts = phasing_stats[1], no_subclonal_muts = phasing_stats[2], no_antiphased_muts = phasing_stats[3],
                   mutspectrum = mutspectrum,
                   stringsAsFactors = F)
  mutannotdf <- cbind(focus[, c("pos", "ref", "alt", "trinuc")], outdf)
  
  return(mutannotdf)
}


get_probs_earlylate <- function(focus) {
  gainloh <- factor(ifelse(focus$major_cn < 2, "noGain", ifelse(focus$minor_cn > 0, "Gain", "GainLOH")), levels = c("noGain", "Gain", "GainLOH"))
  gainlohfrac <- c(table(gainloh)/nrow(focus))
  probs_early_late <- by(data = focus[, c("prob_clonal_early", "prob_clonal_late", "prob_subclonal")], INDICES = gainloh, FUN = function(x) colSums(log(x)))
  probs_early_late <- lapply(probs_early_late, FUN = function(x) if (is.null(x)) c(prob_clonal_early = 0, prob_clonal_late = 0, prob_subclonal = 0) else exp(x)/sum(exp(x)))
  outvect <- c(p_early = sum(0, gainlohfrac[c("Gain", "GainLOH")]*sapply(X = probs_early_late[c("Gain", "GainLOH")], FUN = "[[", "prob_clonal_early")),
               p_late = sum(0, 0, gainlohfrac["GainLOH"]*probs_early_late$GainLOH["prob_clonal_late"]),
               p_sub = sum(gainlohfrac*sapply(X = probs_early_late, FUN = "[[", simplify = T, "prob_subclonal")),
               p_na = sum(gainlohfrac["noGain"]*probs_early_late$noGain[c("prob_clonal_early", "prob_clonal_late")], gainlohfrac["Gain"]*probs_early_late$Gain["prob_clonal_late"], 0))
  return(outvect)
}


Mode <- function(x, na.rm = T) {
  if (na.rm)
    x <- x[!is.na(x)]
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# get the set of segments for the corresponding foci
get_segments <- function(mutations, foci, expand_by = 0) {
  segs <- do.call(rbind, by(data = mutations, INDICES = foci, 
                            FUN = function(x) data.frame(chr = x[1, "chromosome"],
                                                         start = min(x$pos) - expand_by, end = max(x$pos) + expand_by,
                                                         strand = ifelse(sum(x$ref %in% c("C", "T")) >= sum(x$ref %in% c("G", "A")), "+", "-"))))
  colnames(segs) <- c("chr", "start", "end", "strand")
  segs$chr <- as.character(segs$chr)
  segs$foci <- unique(foci)
  return(segs)
}


# get the size of the corresponding CN segments
get_cnsegment_size <- function(chrom, start, end, cnsegments) {
  cnsegment <- subset(cnsegments, subset = (chr == chrom & startpos <= start & endpos >= end), drop = F)
  if(nrow(cnsegment) == 0) 
    return(NA) 
  else 
    return(cnsegment$endpos - cnsegment$startpos)
}

# get the size of the minimal corresponding CN segment (based on consensus breakpoints)
get_minimal_cnsegment_size <- function(chrom, start, end, cnbreaks) {
  cnbreaks_on_chr <- cnbreaks[cnbreaks$chrom == chrom, ]
  cnsize <- cnbreaks_on_chr[which(cnbreaks_on_chr$pos >= end)[1], "pos"] - cnbreaks_on_chr[rev(which(cnbreaks_on_chr$pos <= start))[1], "pos"]
  return(cnsize)
}

# does the focus contain a CN breakpoint
has_cnbreak <- function(chr, startpos, endpos, cnbreaks) {
  cnbreak_in_focus <- ifelse(cnbreaks$chrom == chr,
                   ifelse(startpos <= cnbreaks$pos & endpos >= cnbreaks$pos, T, F),
                   F)
  return(any(cnbreak_in_focus))
}

# has the focus got a nearby SV
# has_sv_nearby <- function(chr, startpos, endpos, svbreaks, proximity_window = 1e4) {
#   nearby_breaks <- subset(svbreaks, subset = (chrom == chr & start >= startpos - proximity_window & start <= endpos + proximity_window), select = id, drop = T)
#   if(length(nearby_breaks) > 0) T else F
# }


# get the distance to the closest SV
get_dist_nearest_sv <- function(chr, startpos, endpos, svbreaks, dist_only = T) {
  svdist <- ifelse(svbreaks$chrom == chr,
                   ifelse(startpos <= svbreaks$start & endpos >= svbreaks$start,
                          0, min(abs(c(startpos - svbreaks$start, endpos - svbreaks$start)))),
                   Inf)
  mindist <- min(svdist, na.rm = T)
  if (!dist_only) {
    sv_pos <- svbreaks[which.min(svdist) , "start"]
    return(data.frame(sv_dist = mindist, sv_pos = ifelse(is.infinite(mindist), NA, sv_pos)))
  }
  return(mindist)
}


## helper functions

#' Robust removal of outliers
#'
#' \code{winsor} returns the non-outlying datapoints
#'
#' @param x Vector containing numerical data
#' @param multiple Number of standard deviations from the mean required to be outlying
#' @param na.rm Should NAs be removed?
#'
#' @return The indices of robustly identified outliers
winsor <- function (x, multiple = 3, na.rm = T) {
  med <- median(x, na.rm = na.rm)
  sc <- mad(x, center = med, na.rm = na.rm) * multiple
  if (sc == 0) 
    return(x)
  else
    return(x[abs(x - med) < sc])
}


# create the reverse complement of a DNA sequence
create_reverse_complement <- function(dna) {
  return(paste0(rev(unlist(strsplit(x = chartr(old = "ACGT", new = "TGCA", x = dna), split = ""))), collapse = ""))
}


# create the complement of a DNA sequence
create_complement <- function(dna) {
  return(x = chartr(old = "ACGT", new = "TGCA", x = dna))
}


# higher level function to compute mutation spectra
analyse_mutationspectra <- function(mutations, sample = "", histol = "", outdir = "", subset) {
  if (nrow(mutations[[subset]]) == 0 ) 
    return(NULL)
  
  if (subset %in% c("apobec", "nonstand", "cpcf"))
    trinuc_freq <- get_trinuc_normalisation_factors(regions = get_segments(mutations = mutations[[subset]], expand_by = 1, foci = mutations[[subset]]$foci))
  else
    trinuc_freq <- unlist(read.delim(file = WGTRINUCCONTFILE)[1,])
  
  plot_mutationspectrum(mutations = mutations[[subset]], trinuc_freq = trinuc_freq, sample = sample, histol = histol, outdir = outdir, suffix = subset)
}


# plot the mutation spectrum
plot_mutationspectrum <- function(mutations, trinuc_freq, sample = "", histol = "", outdir = "", suffix = "") {
  # generate all bases/trinucleotides/mutation types + factors
  bases <- c("A", "C", "G", "T")
  # bases_fact <- factor(bases, levels = bases)
  types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  types_fact <- factor(types, levels = types)
  trinucleotides <- paste0(rep(rep(bases, rep(4,4)), 6),
                           rep(c("C", "T"), c(48, 48)),
                           rep(bases, 24))
  trinucleotides_empty <- paste0(rep(rep(bases, rep(4,4)), 6),
                                 rep(c(" "), c(96)),
                                 rep(bases, 24))
  trinucleotides_mutations <- paste0(paste0(rep(types, rep(16,6))), "_", trinucleotides)
  trinucleotides_mutations_fact <- factor(trinucleotides_mutations, levels = trinucleotides_mutations)
  
  if (!"trinuc" %in% colnames(mutations))
    mutations$trinuc <- as.character(get_trinuc_context(mutations = mutations))
  
  # reverse complement and data augmentation
  revcomp <- data.frame(ref = create_complement(mutations$ref),
                        alt = create_complement(mutations$alt),
                        trinuc = sapply(X = mutations$trinuc, FUN = create_reverse_complement),
                        stringsAsFactors = F)
  mut_full <- factor(ifelse(mutations$ref %in% c("C", "T"), 
                                      paste0(mutations$ref, ">", mutations$alt, "_", mutations$trinuc),
                                      paste0(revcomp$ref, ">", revcomp$alt, "_", revcomp$trinuc)),
                               levels = trinucleotides_mutations)

  # compute frequencies and normalised probabilities
  muttype_freq <- as.vector(table(mut_full)) / nrow(mutations)
  
  counts_normalised <- (muttype_freq / trinuc_freq[trinucleotides]) / sum(muttype_freq / trinuc_freq[trinucleotides], na.rm = T)
  mutdata <- data.frame(type = trinucleotides_mutations_fact, change = paste0(rep(types, rep(16,6))), trinuc = trinucleotides, 
                        freq_plotted = counts_normalised)
  
  # plotting
  p1 <- ggplot(data = mutdata, mapping = aes(x = type)) + 
    geom_bar(mapping = aes(y = freq_plotted, fill = change), show.legend = F, stat="identity",  width = 0.5)
  p1 <- p1 + scale_x_discrete(drop = F, labels = NULL) + scale_fill_discrete(drop=F) +
    theme_bw() + scale_y_continuous(limits = c(-0.25, 1), name = "Mutation type probability") + 
    theme(axis.text.x = element_text(angle = 90), panel.grid.major.x = element_blank(),
          axis.title.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5))
  p1 <- p1 + geom_segment(data = data.frame(start = seq(0.65, 96, 16), end = seq(16.35, 97, 16),
                                            type = types_fact), mapping = aes(x = start, xend = end, y = -0.025, yend = -0.025, color = type),
                          size = 2, show.legend = F)
  p1 <- p1 + geom_text(data = data.frame(type = types_fact, pos = seq(8, 96, 16)),
                       mapping = aes(x = pos, y = -0.25, color = type, label = type),
                       show.legend = F, family = "mono")
  p1 <- p1 + geom_text(data = data.frame(trinuc = trinucleotides, pos = 1:96),
                       mapping = aes(x = pos, y = -0.125, label = trinuc), show.legend = F,
                       angle = 90, size = 3, family = "mono")
  p1 <- p1 + geom_text(data = data.frame(nuc = rep(c("C", "T"), c(48, 48)), pos = 1:96,
                                         col = rep(types_fact, rep(16, 6))),
                       mapping = aes(x = pos, y = -0.125, label = nuc, color = col),
                       show.legend = F, angle = 90, size = 3, fontface = "bold", family = "mono")
  p1 <- p1 + labs(title = paste0(sample, ": ", histol, " - ",nrow(mutations), " focal mutations")) + theme(legend.position="none")
  suppressWarnings(ggsave(filename = file.path(outdir, paste0(sample, "_mutationspectrum_", suffix, ".png")), plot = p1, dpi = 300, width = 10, height = 3))
  return(NULL)
}


# get the trinucleotide content of the whole genome
get_wg_trinuc_normalisation_factors <- function(bsgenome = genome) {
  regions <- data.frame(chr = seqnames(bsgenome), start = 1, end = seqlengths(bsgenome))
  regions <- regions[regions$chr %in% c(1:22, "X", "Y"), ]
  get_trinuc_normalisation_factors(regions = regions, bsgenome = bsgenome)
}


# get the trinucleotide content of specific regions in the genome
get_trinuc_normalisation_factors <- function(regions, bsgenome = genome, overall = T) {
  trinucleotides <- mkAllStrings(alphabet = c("A","C", "T", "G"), width = 3)
  reference_trinucleotides <- grep(pattern = "[ACGT][CT][ACTG]", x = trinucleotides, value = T)
  nonreference_trinucleotides <- as.character(reverseComplement(DNAStringSet(reference_trinucleotides)))

  sequences <- DNAStringSet(getSeq(x = bsgenome, names = regions$chr, start = regions$start, end = regions$end, strand = "+"))
  trinuc_counts_all <- trinucleotideFrequency(x = sequences, as.prob = F)
  trinuc_counts <- trinuc_counts_all[ , reference_trinucleotides, drop = F] + trinuc_counts_all[ , nonreference_trinucleotides, drop = F]

  if (overall) {
    trinuc_counts <- colSums(trinuc_counts)
    return(trinuc_counts / sum(trinuc_counts))
  } else {
    return(trinuc_counts / rowSums(trinuc_counts))
  }
}


# get the trinucleotide context of a set of point mutations
get_trinuc_context <- function(mutations, bsgenome = genome, size = 1) {
  sequences <- getSeq(x = genome, names = mutations$chr, start = mutations$pos - size, end = mutations$pos + size, strand = "+")
  return(sequences)
}


# plot the intermutation distance across the genome
# deprecated and replaced by more extended version
# plot_rainfall <- function(mutations, sample = "", outdir = "") {
#   mutations$chr <- factor(x = mutations$chr, levels = c(1:22, "X", "Y"))
#   mutations <- mutations[order(mutations$chr, mutations$pos),]
#   
#   mutations$limd <- log10(do.call(c, by(data = mutations, INDICES = mutations$chr, FUN = function(x) c(diff(x$pos), NA))))
#   
#   mutations$color <- ifelse(mutations$ref == "C", ifelse(mutations$alt == "A", "C>A", ifelse(mutations$alt == "G", "C>G", "C>T") ),
#                             ifelse(mutations$ref == "G", ifelse(mutations$alt == "T", "C>A", ifelse(mutations$alt == "C", "C>G", "C>T") ),
#                                    ifelse(mutations$ref == "T", ifelse(mutations$alt == "A", "T>A", ifelse(mutations$alt == "C", "T>C", "T>G") ),
#                                           ifelse(mutations$alt == "T", "T>A", ifelse(mutations$alt == "G", "T>C", "T>G")))))
#   
#   mutations$color <- factor(mutations$color, levels = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
#   mutations$strand <- factor(ifelse(mutations$ref %in% c("C", "T"), "+", "-"), levels = c("+", "-"))
#   breaks <- which(is.na(mutations$limd))
#   label_positions <- ( c(0,breaks[-length(breaks)]) + breaks )/2
#   annotations <- data.frame(breaks = breaks, lab_pos = label_positions, label = mutations[is.na(mutations$limd), "chr"])
#   
#   p1 <- ggplot(data = mutations, mapping = aes(x = 1:length(limd), y = limd)) +
#     geom_point(mapping = aes(colour = color, shape = strand), show.legend = F, alpha = .75, size = 1)
#   p1 <- p1 + scale_colour_manual(values = c("C>A" = "blue", "C>G" = "black", "C>T" = "red", "T>A" = "purple", "T>C" = "yellow", "T>G" = "green"))
#   p1 <- p1 + theme_minimal() + theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), 
#                                      axis.text.x = element_blank(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5), 
#                                      legend.position = "none", axis.ticks = element_blank())
#   p1 <- p1 + geom_vline(xintercept = breaks, color = "grey", alpha = 0.5)
#   p1 <- p1 + geom_text(data = annotations, mapping = aes(x = lab_pos, y = -0.5, label = label), size = 3)
#   p1 <- p1 + labs(title = paste0(sample, ": ", nrow(mutations), " mutations"), y = "log(intermutation distance)")
#   ggsave(filename = file.path(outdir, paste0(sample, "_rainfall.png")), plot = p1, dpi = 300, width = 10, height = 3)
#   return(NULL)
# }


# read all required PCAWG data for a sample
read_all_data <- function(sample_id, pcfdir, svdir, snv_mnvdir, cnbreaksdir, rhopsifile, cclustdir, cccfdir, ctimingdir, phasingdir) {
  
  pcf_file <- list.files(path = pcfdir, pattern = sample_id, full.names = T)
  muts_all_pcf <- read.delim(file = pcf_file, header = T, sep = "\t", as.is = T,
                             colClasses = c("character", "integer", rep("character", 2), "integer"))
  colnames(muts_all_pcf) <- c("chromosome", "pos", "ref", "alt", "foci")
  muts_all_pcf$chromosome <- as.character(muts_all_pcf$chromosome)
  muts_all_pcf$trinuc <- as.character(get_trinuc_context(mutations = muts_all_pcf, bsgenome = genome, size = 1))
  
  svfile <- list.files(path = svdir, pattern = paste0(sample_id, "[0-9A-Za-z_.]*", ".somatic.sv.bedpe.gz$"), full.names = T)
  svs <- tryCatch(
    {
      read.delim(file = gzfile(svfile), header = T, sep = "\t", as.is = T) 
    }, error = function(err) {
      print(paste("ERROR: sample ", sample_id, err))
      return(NULL)
    }
  )
  
  if (is.null(svs)) {
    svbreakpoints <- NULL
  } else if (nrow(svs) == 0) {
    svbreakpoints <- data.frame(id = character(), bp_side = numeric(), chrom = character(), start = integer())
  } else {
    colnames(svs) <- c("chrom.1", "start.1", "end.1", "chrom.2", "start.2", "end.2", "id", "pe_support", "strand.1", "strand.2", "svclas", "svmethod")
    svbreakpoints <- reshape(svs[ , c(1,2,4,5,7)], direction = "long", varying = 1:4, timevar = "bp_side")
    row.names(svbreakpoints) <- NULL
  }
  # 
  # if (!is.null(svs)) {
  #   colnames(svs) <- c("chrom.1", "start.1", "end.1", "chrom.2", "start.2", "end.2", "id", "pe_support", "strand.1", "strand.2", "svclas", "svmethod")
  # }
  # if (nrow(svs) == 0) {
  #   svbreakpoints <- data.frame(id = character(), bp_side = numeric(), chrom = character(), start = integer())
  # } else {
  #   svbreakpoints <- reshape(svs[ , c(1,2,4,5,7)], direction = "long", varying = 1:4, timevar = "bp_side")
  # }
  # row.names(svbreakpoints) <- NULL
  
  snvfile <- list.files(path = snv_mnvdir, pattern = paste0(sample_id, ".*somatic.snv_mnv.vcf.gz$"), full.names = T, recursive = T)
  allmuts <- tryCatch(
    {
      read.delim(file = gzfile(snvfile), header = F, sep = "\t", as.is = T, comment.char = "#")[, c(1,2,4,5)]
    }, error = function(err) {
      print(paste("ERROR: sample ", sample_id, err))
      return(data.frame(chr = character(), pos = integer(), ref = character(), alt = character()))
    }
  )
  colnames(allmuts) <- c("chr", "pos", "ref", "alt")
  
  cnbreaks_file <- list.files(path = cnbreaksdir, pattern = sample_id, full.names = T)
  cnbreakpoints <- read.delim(file = gzfile(cnbreaks_file), header = T, sep = "\t", as.is = T)
  
  rhopsiall <- read.delim(file = gzfile(rhopsifile), header = T, sep = "\t", as.is = T)
  rho <- rhopsiall[rhopsiall$samplename == sample_id, "purity"]
  
  # mcn_file <- list.files(path = mcndir, pattern = sample_id, full.names = T)
  # multiplicities <- read.delim(file = gzfile(mcn_file), header = T, sep = "\t", as.is = T)
  # multiplicities$ccf <- multiplicities$mutation.copy.number / ifelse(multiplicities$multiplicity == 0, multiplicities$mutation.copy.number, multiplicities$multiplicity)
  
  # clusters_file <- list.files(path = clustdir, pattern = paste0(sample_id, "_subclonal_structure.txt.gz"), full.names = T)
  # clusters <- read.delim(file = gzfile(clusters_file), header = T, sep = "\t", as.is = T)
  # clusters$ccf <- clusters$proportion / rho
  
  phasing_file <- file.path(phasingdir, paste0(sample_id, "_phasedmuts.txt.gz"))
  phasing <- tryCatch(
    {
      phasing <- read.delim(file = gzfile(phasing_file), header = T, sep = "\t", as.is = T)
      if (ncol(phasing) == 12) {
      phasing <- read.delim(file = gzfile(phasing_file), header = T, sep = "\t", as.is = T, colClasses = c("character", "integer", rep("character", 2),
                                                                                                           "integer", rep("character", 2), rep("integer", 4), 
                                                                                                           "character"))
      } else {
        phasing <- read.delim(file = gzfile(phasing_file), header = T, sep = "\t", as.is = T, colClasses = c("character", "integer", rep("character", 2),
                                                                                                             "integer", rep("character", 2), rep("numeric", 2), rep("integer", 4), 
                                                                                                             "character"))[, c(1:7, 10:14)]
      }
    }, error = function(err) {
      print(paste("ERROR: sample ", SAMPLE, err))
      phasingdf <- data.frame(Chr = character(), Pos1 = integer(), Ref1 = character(), Var1 = character(), Pos2 = integer(), Ref2 = character(), Var2 = character(),
                            Num_WT_WT = integer(), Num_Mut_Mut = integer(), Num_Mut_WT = integer(), Num_WT_Mut = integer(), phasing = character())
      return(phasingdf)
      }
    )

  # clustering consensus data
  cassignments <- read.delim(file = gzfile(file.path(cclustdir, paste0(sample_id, "_cluster_assignments.txt.gz"))), header = T, sep = "\t", as.is = T)
  cccfs <- read.delim(file = gzfile(file.path(cccfdir, paste0(sample_id, "_mutation_ccf.txt.gz"))), header = T, sep = "\t", as.is = T)
  ctiming <- read.delim(file = gzfile(file.path(ctimingdir, paste0(sample_id, "_prob_gained.txt"))), header = T, sep = "\t", as.is = T)
  colnames(ctiming) <- c("chromosome", "position", "mut_type", "timing", "chromosome2", "position2", "svid", "prob_clonal_early", "prob_clonal_late", "prob_subclonal")
  
  if (all(c(nrow(cassignments), nrow(cccfs), nrow(ctiming)) == nrow(cassignments)) ) {
    mutsdf <- cbind(cassignments[ , c("chromosome", "position", "mut_type", grep(pattern = "cluster_", x = colnames(cassignments), value = T))],
                    ctiming[ , c("timing", "prob_clonal_early", "prob_clonal_late", "prob_subclonal")],
                    cccfs[ , c("ccf", "major_cn", "minor_cn", "mcn", "mult")])
    mutsdf <- mutsdf[mutsdf$mut_type == "SNV", c("chromosome", "position", "ccf", "major_cn", "minor_cn", "mcn", "mult", 
                                                 grep(pattern = "cluster_", x = colnames(mutsdf), value = T), 
                                                 "timing", "prob_clonal_early", "prob_clonal_late", "prob_subclonal")]
    # if a mut pos is duplicated, see if there's one which has a ccf
    duppos <- mutsdf[duplicated(mutsdf$position), "position"]
    if (length(duppos) > 0) {
      dupmuts <- mutsdf[which(mutsdf$position %in% duppos), ]
      dedupmuts <- do.call(rbind, by(data = dupmuts, INDICES = dupmuts$position, FUN = function(x) if (any(!is.na(x$ccf))) x[!is.na(x$ccf), ][1, ] else x[1, ]))
      mutsdf <- rbind(mutsdf[-which(mutsdf$position %in% duppos), ], dedupmuts)
    }
    
    cmuts_all_pcf <- merge(x = muts_all_pcf[, c("chromosome", "pos", "ref", "alt", "trinuc", "foci")], #, "mut.count", "WT.count")], 
                           y = mutsdf, by.x = c("chromosome", "pos"), by.y = c("chromosome", "position"), all.x = T)
  } else {
    print("read_all_data: should not be here!")
    cmuts_all_pcf <- merge(x = muts_all_pcf[, c("chromosome", "pos", "ref", "alt", "trinuc", "foci")], #, "mut.count", "WT.count")], 
                           y = cccfs[cccfs$type == "SNV", c("chromosome", "position", "ccf", "major_cn", "minor_cn", "mcn", "mult")],
                           by.x = c("chromosome", "pos"), by.y = c("chromosome", "position"), all.x = T)
    cmuts_all_pcf <- merge(x = cmuts_all_pcf,
                           y = cassignments[cassignments$mut_type == "SNV", c("chromosome", "position", grep(pattern = "cluster_", x = colnames(cassignments), value = T))],
                           by.x = c("chromosome", "pos"), by.y = c("chromosome", "position"), all.x = T)
    cmuts_all_pcf <- merge(x = cmuts_all_pcf, 
                           y = ctiming[ctiming$mut_type == "SNV", c("chromosome", "position", "timing", "prob_clonal_early", "prob_clonal_late", "prob_subclonal")],
                           by.x = c("chromosome", "pos"), by.y = c("chromosome", "position"), all.x = T)
    cmuts_all_pcf <- cmuts_all_pcf[!duplicated(paste0(cmuts_all_pcf$chromosome, "_", cmuts_all_pcf$position, "_", cmuts_all_pcf$ref, "/", cmuts_all_pcf$alt)),]
  }
  
  # colnames(cmuts_all_pcf)[grep(colnames(cmuts_all_pcf), pattern = "cluster_1")] <- "p_clonal"

  csubcl <- read.delim(file = gzfile(file.path(cclustdir, paste0(sample_id, "_subclonal_structure.txt.gz"))), header = T, sep = "\t", as.is = T)
  
  # fixing negative/regularising pGain/pSingle/pSubclonal
  # values < 1e-15 (approx .Machine$double.eps) set to 1e-15 and renormalised such that pSingle + pGain + pSub = 1
  cmuts_all_pcf$prob_clonal_early <- ifelse(cmuts_all_pcf$major_cn > 1 & cmuts_all_pcf$prob_clonal_early < 1e-15, 1e-15, cmuts_all_pcf$prob_clonal_early)
  cmuts_all_pcf$prob_clonal_late <- ifelse(cmuts_all_pcf$prob_clonal_late < 1e-15, 1e-15, cmuts_all_pcf$prob_clonal_late)
  # if >= 1 subclone present, also regularise these probs
  if (nrow(csubcl) > 1) {
    cmuts_all_pcf$prob_subclonal <- ifelse(cmuts_all_pcf$prob_subclonal < 1e-15, 1e-15, cmuts_all_pcf$prob_subclonal)
  }
  # renorm
  cmuts_all_pcf[ , c("prob_clonal_early", "prob_clonal_late", "prob_subclonal")] <- cmuts_all_pcf[ , c("prob_clonal_early", "prob_clonal_late", "prob_subclonal")] / rowSums(cmuts_all_pcf[ , c("prob_clonal_early", "prob_clonal_late", "prob_subclonal")])

  timing_summary <- c(table(factor(ctiming$timing, levels = c("clonal [NA]", "clonal [early]", "clonal [late]", "subclonal"))))
  
  return(list(pcf = muts_all_pcf,
              cpcf = cmuts_all_pcf,
              sv = svbreakpoints,
              snv = allmuts,
              cn = cnbreakpoints,
              rho = rho,
              # mcn = multiplicities,
              # subcl = clusters,
              phasing = phasing,
              csubcl = csubcl,
              cass = cassignments,
              cccfs = cccfs,
              ctiming = ctiming,
              timing_summary = timing_summary))
}


# read and fix combined histology file
read_histology <- function(histologyfile, melafile = "/srv/shared/vanloo/ICGC_annotations/icgc_melanoma_new_label.txt") {
  histology_all <- read.delim(file = histologyfile, as.is = T)
  melaannots <- read.delim(file = melafile, as.is = T)
  melaannots[melaannots$subtype == "Cutaneous", "subtype"] <- "Cut"
  histology_all$histology_abbreviation <- ifelse(histology_all$histology_abbreviation == "Kidney-RCC",
                                                    ifelse(grepl("papillary", histology_all$histology_tier4), "Kidney-RCC-Pap","Kidney-RCC-Clear"), histology_all$histology_abbreviation)
  histology_all <- merge(x = histology_all, y = melaannots[, c("icgc_aliquot", "subtype")], by.x = "samplename", by.y = "icgc_aliquot", all.x = T)
  histology_all$histology_abbreviation <- ifelse(is.na(histology_all$subtype), histology_all$histology_abbreviation, paste0(histology_all$histology_abbreviation, "-", histology_all$subtype))
  histology_all <- histology_all[, !colnames(histology_all) %in% "subtype"]
  return(histology_all)
}


# combine PCAWG histology and release data files
get_clean_histology <- function(histologyfile, releasedatafile, exlude_greylisted = T) {
  pcawg_histology <- read.delim(file = histologyfile, header = T, sep = "\t", as.is = T)
  pcawg_histology <- pcawg_histology[pcawg_histology$specimen_library_strategy == "WGS", ]
  
  pcawg_release_data <- read.delim(file = releasedatafile, header = T, sep = "\t", as.is = T)
  pcawg_release_data_redup <- pcawg_release_data[rep(1:nrow(pcawg_release_data), pcawg_release_data$tumor_wgs_specimen_count),]
  pcawg_release_data_redup$tumor_wgs_icgc_specimen_id <- unlist(strsplit(pcawg_release_data$tumor_wgs_icgc_specimen_id, split = ","))
  pcawg_release_data_redup$tumor_wgs_aliquot_id <- unlist(strsplit(pcawg_release_data$tumor_wgs_aliquot_id, split = ","))
  
  pcawg_annot <- merge(x = pcawg_histology, y = pcawg_release_data_redup, by.x = "icgc_specimen_id",
                       by.y = "tumor_wgs_icgc_specimen_id")
  
  if (exlude_greylisted)
    pcawg_annot_trim <- pcawg_annot[pcawg_annot$donor_wgs_included_excluded != "Excluded", c("tumor_wgs_aliquot_id", "histology_abbreviation", "donor_wgs_included_excluded")]
  else 
    pcawg_annot_trim <- pcawg_annot[, c("tumor_wgs_aliquot_id", "histology_abbreviation", "donor_wgs_included_excluded")]
  
  # test: any(duplicated(pcawg_annot_trim))
  row.names(pcawg_annot_trim) <- pcawg_annot_trim$tumor_wgs_aliquot_id
  return(pcawg_annot_trim)
}


# creates the draw/no_draw regions for each chromosome
get_regions_to_draw <- function(segments_on_chrom) {
  segments_on_chrom <- segments_on_chrom[order(segments_on_chrom$start), ]
  outer_borders <- c(min(segments_on_chrom$start), max(segments_on_chrom$end))
  to_draw <- paste0("hs", segments_on_chrom[1, "chr"], ":", outer_borders[1] - 1000, "-", outer_borders[2] + 1000)
  shifted_df <- data.frame(gapstart = segments_on_chrom[-nrow(segments_on_chrom), "end"],
                           gapend = segments_on_chrom[-1, "start"])
  shifted_df <- shifted_df[shifted_df$gapend - shifted_df$gapstart >= 5000, ]
  not_to_draw <- paste0(apply(shifted_df, MARGIN = 1, FUN = function(x) paste0("-hs", segments_on_chrom[1, "chr"], ":", x[1] + 1000, "-", x[2] - 1000)), collapse = ";")
  return(list(draw = to_draw, no_draw = not_to_draw))
}


# make circos plot of a set of foci
plot_kataegis_circos <- function(mutations, summarized_foci, sample_id, circosdir, outdir, svbreaks, cnbreaks, suffix = "") {
  ## For plotting all of the identified segments
  mutations$color <- ifelse(mutations$ref == "C", ifelse(mutations$alt == "A", "blue", ifelse(mutations$alt == "G", "black", "red") ),
                               ifelse(mutations$ref == "G", ifelse(mutations$alt == "T", "blue", ifelse(mutations$alt == "C", "black", "red") ),
                                      ifelse(mutations$ref == "T", ifelse(mutations$alt == "A", "purple", ifelse(mutations$alt == "C", "yellow", "green") ),
                                             ifelse(mutations$alt == "T", "purple", ifelse(mutations$alt == "G", "yellow", "green")))))
  mutations$strand <- ifelse(mutations$ref %in% c("C", "T"), "circle", "triangle")

  ## data reshaping and PLOTTING
  # write mutations to be plotted + color and CCF, remove NAs
  circos_data <- data.frame(chrom = paste0("hs", mutations$chromosome),
                            start = mutations$pos, end = mutations$pos,
                            ccf = mutations$ccf, 
                            color = paste0("color=", mutations$color, ",glyph=", mutations$strand), stringsAsFactors = F)

  circos_data[is.na(circos_data$ccf), "ccf"] <- 0
  circos_katmut_file <- file.path(circosdir, "data", paste0(sample_id, "_katmut.txt"))
  write.table(x = circos_data, file = circos_katmut_file, sep = " ", col.names = F, row.names = F, quote = F)
  
  # write "fits" of the mean CCF per focus and color according to clonality
  circos_highlights <- data.frame(chr = paste0("hs", summarized_foci$chr), start = summarized_foci$start -250, 
                                  end = summarized_foci$end+250,
                                  ccf = paste0("offset=", round(summarized_foci$katccf, 3) * 200, "p,fill_color=",
                                               ifelse(is.na(summarized_foci$timing), "grey",
                                                      ifelse(summarized_foci$timing == "clonal [NA]", "blue",
                                                             ifelse(summarized_foci$timing == "clonal [early]", "green",
                                                                    ifelse(summarized_foci$timing == "clonal [late]", "purple",
                                                                           ifelse(summarized_foci$timing == "subclonal", "red", "grey")))))))
  circos_highlights_file <- file.path(circosdir, "data", paste0(sample_id, "_meanccfs.txt"))
  write.table(x = circos_highlights, file = circos_highlights_file, sep = " ", col.names = F, row.names = F, quote = F)
  
  # write focal signatures into text track
  circos_signatures <- data.frame(chr = paste0("hs", summarized_foci$chr), start = summarized_foci$start, 
                                  end = summarized_foci$end,
                                  signature = sub(pattern = "SBS", replacement = "", x = summarized_foci$active_sig))
  circos_signatures_file <- file.path(circosdir, "data", paste0(sample_id, "_signatures.txt"))
  write.table(x = circos_signatures, file = circos_signatures_file, sep = " ", col.names = F, row.names = F, quote = F)
  
  
  # determine which regions to draw and which to exclude ...
  circos_to_draw <- unlist(by(data = summarized_foci, INDICES = summarized_foci$chr, FUN = get_regions_to_draw))
  chromosomes_draw <- paste0(circos_to_draw[seq(1, length(circos_to_draw), 2)], collapse = ";")
  chromosomes_nodraw <- paste0(circos_to_draw[intersect(seq(2, length(circos_to_draw), 2), which(circos_to_draw != ""))], collapse = ";")
  
  # ... and modify circos conf accordingly
  circos_conf_template <- file.path(circosdir, "circos.conf")
  circos_conf <- file.path(circosdir, paste0(sample_id, "_circos.conf"))
  sedcommand <- paste0("sed -e \'s/CHROMOSOMESTODRAW/", chromosomes_draw, 
                       "/g; s/CHROMOSOMESNODRAW/", chromosomes_nodraw,
                       "/g; s/SAMPLEID/", sample_id, "/g\' ", circos_conf_template, " > ", circos_conf)
  system(sedcommand)
  
  # reformat matching SV data for circos
  if (nrow(svbreaks) > 0) {
    circos_svs <- data.frame(chr = paste0("hs", svbreaks$chrom), start = svbreaks$start, end = svbreaks$start)
  } else {
    circos_svs <- data.frame(chr = character(), start = integer(), end = integer())
  }
  circos_svs_file <- file.path(circosdir, "data", paste0(sample_id, "_svs.txt"))
  write.table(x = circos_svs, file = circos_svs_file, sep = " ", col.names = F, row.names = F, quote = F)
  
  # reformat matching CNV data for circos
  circos_cnvs <- data.frame(chr = paste0("hs", cnbreaks$chr), start = cnbreaks$pos, end = cnbreaks$pos)
  circos_cnvs_file <- file.path(circosdir, "data", paste0(sample_id, "_cnvs.txt"))
  write.table(x = circos_cnvs, file = circos_cnvs_file, sep = " ", col.names = F, row.names = F, quote = F)
  
  ## execute circos command and cleanup
  circoscmd <- paste0("circos -conf ", circos_conf, " -outputdir ", outdir, " -outputfile ", sample_id, "_circos_", suffix)#," -nosvg")
  
  system(command = circoscmd, ignore.stdout = T, ignore.stderr = T)
  unlink(x = c(circos_conf, circos_cnvs_file, circos_highlights_file, circos_svs_file, circos_katmut_file))
}


generate_bases_types_trinuc <- function() {
  # generate all bases/trinucleotides/mutation types + factors
  bases <- c("A", "C", "G", "T")
  # bases_fact <- factor(bases, levels = bases)
  types <- c("C.A", "C.G", "C.T", "T.A", "T.C", "T.G")
  types_gt <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  types_fact <- factor(types, levels = types)
  trinucleotides <- paste0(rep(rep(bases, rep(4,4)), 6),
                           rep(c("C", "T"), c(48, 48)),
                           rep(bases, 24))
  trinucleotides_empty <- paste0(rep(rep(bases, rep(4,4)), 6),
                                 rep(c(" "), c(96)),
                                 rep(bases, 24))
  trinucleotides_mutations <- paste0(paste0(rep(types, rep(16,6))), "_", trinucleotides)
  trinucleotides_mutations_gt <- paste0(paste0(rep(types_gt, rep(16,6))), "_", trinucleotides)
  trinucleotides_mutations_fact <- factor(trinucleotides_mutations, levels = trinucleotides_mutations)
  return(list(bases = bases, types = types, types_fact = types_fact, trinucleotides = trinucleotides,
              trinucleotides_empty = trinucleotides_empty, trinucleotides_mutations = trinucleotides_mutations,
              trinucleotides_mutations_gt = trinucleotides_mutations_gt,
              trinucleotides_mutations_fact = trinucleotides_mutations_fact))
}


get_focal_mutspectrum <- function(mutations, trinuc_freq) {
  base_type_trinuc_info <- generate_bases_types_trinuc()
  
  if (!"trinuc" %in% colnames(mutations))
    mutations$trinuc <- as.character(get_trinuc_context(mutations = mutations))
  
  # reverse complement and data augmentation
  revcomp <- data.frame(ref = create_complement(mutations$ref),
                        alt = create_complement(mutations$alt),
                        trinuc = sapply(X = mutations$trinuc, FUN = create_reverse_complement),
                        stringsAsFactors = F)
  mut_full <- factor(ifelse(mutations$ref %in% c("C", "T"), 
                            paste0(mutations$ref, ".", mutations$alt, "_", mutations$trinuc),
                            paste0(revcomp$ref, ".", revcomp$alt, "_", revcomp$trinuc)),
                     levels = base_type_trinuc_info[["trinucleotides_mutations"]])
  
  # compute frequencies and normalised probabilities
  muttype_freq <- as.vector(table(mut_full)) / nrow(mutations)
  
  counts_normalised <- (muttype_freq / trinuc_freq[base_type_trinuc_info[["trinucleotides"]]]) / sum(muttype_freq / trinuc_freq[base_type_trinuc_info[["trinucleotides"]]], na.rm = T)
  # mutdata <- data.frame(type = trinucleotides_mutations_fact, change = paste0(rep(types, rep(16,6))), trinuc = trinucleotides, 
                        # freq_plotted = counts_normalised)
  names(counts_normalised) <- base_type_trinuc_info[["trinucleotides_mutations"]]
  return(counts_normalised)
}
# 
# get_sample_mutspectrum <- function(sample, histology, resultsdir, trinuc_freq = trinuc_freq) {
#   apobec_resultsfile <- file.path(resultsdir, sample, paste0(sample, "_kataegis_apobec_muts.txt"))
#   nonstand_resultsfile <- file.path(resultsdir, sample, paste0(sample, "_kataegis_nonstand_muts.txt"))
#   if (!file.exists(apobec_resultsfile) && !file.exists(nonstand_resultsfile)) {
#     return(data.frame())
#   }
#   if (file.exists(apobec_resultsfile)) {
#     apobec_results_muts <- read.delim(file = apobec_resultsfile, header = T, sep = "\t", as.is = T)
#     apobec_results <- data.frame(do.call(rbind, by(data = apobec_results_muts, INDICES = apobec_results_muts$foci, FUN = get_focal_mutspectrum, trinuc_freq = trinuc_freq)))
#     apobec_results$type <- "APOBEC3A/B"
#   }
#   if (file.exists(nonstand_resultsfile)) {
#     nonstd_results_muts <- read.delim(file = nonstand_resultsfile, header = T, sep = "\t", as.is = T)
#     nonstd_results <- data.frame(do.call(rbind, by(data = nonstd_results_muts, INDICES = nonstd_results_muts$foci, FUN = get_focal_mutspectrum, trinuc_freq = trinuc_freq)))
#     nonstd_results$type <- "nonstd"
#   }
#   all_results$sample <- sample
#   all_results$histology <- histology
#   return(all_results)
# }

# plot the intermutation distance across the genome
plot_rainfall <- function(mutations, sample = "", outdir = "", svbreaks, plot_position = F, bsgenome = genome) {
  mutations$chr <- factor(x = mutations$chr, levels = c(1:22, "X", "Y"))
  mutations <- mutations[order(mutations$chr, mutations$pos), ]
  
  chrs_with_multiple_muts <- names(which(by(data = mutations, INDICES = mutations$chr, FUN = function(x) nrow(x) >= 2)))
  mutations <- mutations[mutations$chr %in% chrs_with_multiple_muts, ]
  
  cumulative_genomic_pos <- cumsum(c(0, as.numeric(seqlengths(genome)[c(1:22, "X")])))
  names(cumulative_genomic_pos) <- c(1:22, "X", "Y")
  
  if (!is.null(svbreaks) && nrow(svbreaks) > 0)
    svbreaks$cumulpos <- svbreaks$start + cumulative_genomic_pos[svbreaks$chrom]
  
  mutations$limd <- log10(do.call(c, by(data = mutations, INDICES = mutations$chr, FUN = function(x) c(diff(x$pos), NA))))
  mutations$color <- ifelse(mutations$ref == "C", ifelse(mutations$alt == "A", "C>A", ifelse(mutations$alt == "G", "C>G", "C>T") ),
                            ifelse(mutations$ref == "G", ifelse(mutations$alt == "T", "C>A", ifelse(mutations$alt == "C", "C>G", "C>T") ),
                                   ifelse(mutations$ref == "T", ifelse(mutations$alt == "A", "T>A", ifelse(mutations$alt == "C", "T>C", "T>G") ),
                                          ifelse(mutations$alt == "T", "T>A", ifelse(mutations$alt == "G", "T>C", "T>G")))))
  
  mutations$color <- factor(mutations$color, levels = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
  mutations$strand <- factor(ifelse(mutations$ref %in% c("C", "T"), "+", "-"), levels = c("+", "-"))
  
  if (!plot_position) {
    breaks <- which(is.na(mutations$limd))
    label_positions <- ( c(0,breaks[-length(breaks)]) + breaks )/2
    annotations <- data.frame(breaks = breaks, lab_pos = label_positions, label = mutations[is.na(mutations$limd), "chr"])
    
    p1 <- ggplot(data = mutations, mapping = aes(x = 1:length(limd), y = limd)) +
      geom_point(mapping = aes(colour = color, shape = strand), show.legend = F, alpha = .75, size = 1)
    p1 <- p1 + scale_colour_manual(values = c("C>A" = "blue", "C>G" = "black", "C>T" = "red", "T>A" = "purple", "T>C" = "yellow", "T>G" = "green"))
    p1 <- p1 + theme_minimal() + theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), 
                                       axis.text.x = element_blank(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5), 
                                       legend.position = "none", axis.ticks = element_blank())
    p1 <- p1 + geom_vline(xintercept = breaks, color = "grey", alpha = 0.5)
    p1 <- p1 + geom_text(data = annotations, mapping = aes(x = lab_pos, y = -0.5, label = label), size = 3)
    p1 <- p1 + labs(title = paste0(sample, ": ", nrow(mutations), " mutations"), y = "log(intermutation distance)")
    
  } else {
    mutations$cumulpos <- mutations$pos + cumulative_genomic_pos[mutations$chr]
    mutations$avcumpos <- do.call(c, by(data = mutations$cumulpos, INDICES = mutations$chr,
                                          FUN = function(x) filter(x, rep(0.5, 2), sides=2)))
    
    breaks <- cumsum(c(0, as.numeric(seqlengths(genome)[c(1:22, "X", "Y")])))
    label_positions <- filter(breaks, rep(.5, 2), sides = 2)[-length(breaks)]
    annotations <- data.frame(lab_pos = label_positions, label = c(1:22, "X", "Y"))
    
    p1 <- ggplot(data = mutations, mapping = aes(x = avcumpos, y = limd)) +
      geom_point(mapping = aes(colour = color, shape = strand), show.legend = F, alpha = .75, size = 1, na.rm=T)
    p1 <- p1 + scale_colour_manual(values = c("C>A" = "blue", "C>G" = "black", "C>T" = "red", "T>A" = "purple", "T>C" = "yellow", "T>G" = "green"))
    p1 <- p1 + theme_minimal() + theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), 
                                       axis.text.x = element_blank(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5), 
                                       legend.position = "none", axis.ticks = element_blank())
    p1 <- p1 + geom_vline(xintercept = breaks, color = "grey", alpha = 0.5)
    
    if (!is.null(svbreaks) && nrow(svbreaks) > 0)
      p1 <- p1 + geom_point(data = svbreaks, mapping = aes(x = cumulpos, y = -1), shape = "|")
    
    p1 <- p1 + geom_text(data = annotations, mapping = aes(x = lab_pos, y = -0.5, label = label), size = 3, na.rm = T)
    p1 <- p1 + labs(title = paste0(sample, ": ", nrow(mutations), " mutations"), y = "log(intermutation distance)")
  }
  ggsave(filename = file.path(outdir, paste0(sample, "_rainfall.png")), plot = p1, dpi = 300, width = 10, height = 3)
  # ggsave(filename = file.path(outdir, paste0(sample, "_rainfall.pdf")), plot = p1, width = 10, height = 3)
  return(NULL)
}


get_intermut_dist <- function(mutations, q = 0.01) {
  mutations$chr <- factor(x = mutations$chr, levels = c(1:22, "X", "Y"))
  mutations <- mutations[order(mutations$chr, mutations$pos), ]
  
  chrs_with_multiple_muts <- names(which(by(data = mutations, INDICES = mutations$chr, FUN = function(x) nrow(x) >= 2)))
  mutations <- mutations[mutations$chr %in% chrs_with_multiple_muts, ]
  
  imd <- do.call(c, by(data = mutations, INDICES = mutations$chr, FUN = function(x) diff(x$pos)))
  quant <- quantile(x = imd, na.rm = T, probs = q)
  return(quant)
}


get_active_sigs_per_sample <- function(sigs_per_sample_file = SIGSINSAMPLESFILE, mergesigs = NULL) {
  sigs_per_sample <- read.csv(file = sigs_per_sample_file, as.is = T)
  sigs_only <- sigs_per_sample[ , grep(pattern = "SBS", x = colnames(sigs_per_sample))]
  
  active_sigs_per_sample <- apply(X = sigs_only, MARGIN = 1, FUN = function(x, cols) cols[which(x > 0)], cols = colnames(sigs_only))
  names(active_sigs_per_sample) <- sigs_per_sample$Sample.Name
  
  if (!is.null(mergesigs)) {
    for (sig in mergesigs) {
      active_sigs_per_sample <- lapply(X = active_sigs_per_sample, FUN = function(x, sig) unique(sub(pattern = paste0(sig, "[abcd]"), x = x, replacement = sig)), sig = sig)
    }
  }
  
  return(active_sigs_per_sample)
}



load_signatures <- function(signatures_file = SIGFILE, pad = F, mergesigs = NULL) {
  
  # read data
  signatures <- read.csv(file = signatures_file, as.is = T)
  sig_names <- colnames(signatures)[-1]
  if (pad) {
    signatures[signatures == 0] <- 2.23E-16
    signatures[ , -1] <- signatures[ , -1]/colSums(signatures[ ,-1])
  }
  
  if (!is.null(mergesigs)) {
    for (sigid in mergesigs) {
      colids <- grepl(pattern = paste0("SBS", sigid, "[ab]"), x = colnames(signatures))
      signatures[, paste0("SBS", sigid)] <- rowSums(signatures[, colids])/sum(colids)
      signatures <- signatures[, !colids]
    }
     
  }

  # general statics: generate all bases/trinucleotides/mutation types + factors
  bases <- c("A", "C", "G", "T")
  # bases_fact <- factor(bases, levels = bases)
  types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  types_fact <- factor(types, levels = types)
  trinucleotides <- paste0(rep(rep(bases, rep(4,4)), 6),
                           rep(c("C", "T"), c(48, 48)),
                           rep(bases, 24))
  trinucleotides_empty <- paste0(rep(rep(bases, rep(4,4)), 6),
                                 rep(c(" "), c(96)),
                                 rep(bases, 24))
  trinucleotides_mutations <- paste0(paste0(rep(types, rep(16,6))), "_", trinucleotides)
  trinucleotides_mutations_fact <- factor(trinucleotides_mutations, levels = trinucleotides_mutations)
  
  signatures$Mutation.Type <- paste0(substr(x = signatures$Mutation.Type, start = 2, stop = 2), 
                                     substr(x = signatures$Mutation.Type, start = 4, stop = 5), "_",
                                     substr(x = signatures$Mutation.Type, start = 1, stop = 3))
  
  # reorder input according to ludmil plots
  signatures$trinuc_muts <- factor(signatures$Mutation.Type, levels = trinucleotides_mutations)
  signatures <- signatures[order(signatures$trinuc_muts), ]
  sig_list <- split(x = t(signatures[, -c(1, ncol(signatures))]), f = 1:(ncol(signatures)-2))
  sig_list <- lapply(X = sig_list, FUN = function(x) {
    names(x) <- trinucleotides_mutations
    return(x)
  })
  names(sig_list) <- colnames(signatures)[-c(1, ncol(signatures))]
  return(sig_list)
}


adjust_signatures_local_trinuc <- function(signatures, mutations, trinuc_genome = trinuc_freq) {
  # signatures <- load_signatures()
  # mutations <- alldata[["pcf"]]
  # trinuc_genome <- trinuc_freq
  
  base_type_trinuc_info <- generate_bases_types_trinuc()
  
  foci <- get_segments(mutations = mutations, foci = mutations$foci, expand_by = 1)
  trinuc_foci <- get_trinuc_normalisation_factors(regions = foci, overall = F)
  correction_factors <- t(t(trinuc_foci) / trinuc_genome)
  # correction_factors[is.infinite(correction_factors)] <- 0
  correction_factors_matched <- correction_factors[, base_type_trinuc_info[["trinucleotides"]]]
  # temp <- lapply(X = signatures, FUN = function(y) y*x / (sum(y*x)))
  signatures_renorm <- apply(X = correction_factors_matched, MARGIN = 1, 
                             FUN = function(x) lapply(X = signatures, 
                                                      FUN = function(y) y*x / (sum(y*x))))
  return(signatures_renorm)
}


compute_likelihoods <- function(signatures_renorm, mutations) {
  
  base_type_trinuc_info <- generate_bases_types_trinuc()
  
  revcomp <- data.frame(ref = as.character(reverseComplement(DNAStringSet(mutations$ref))),
                        alt = as.character(reverseComplement(DNAStringSet(mutations$alt))),
                        trinuc = as.character(reverseComplement(DNAStringSet(mutations$trinuc))),
                        stringsAsFactors = F)
  mut_full <- factor(ifelse(mutations$ref %in% c("C", "T"), 
                            paste0(mutations$ref, ">", mutations$alt, "_", mutations$trinuc),
                            paste0(revcomp$ref, ">", revcomp$alt, "_", revcomp$trinuc)),
                     levels = base_type_trinuc_info[["trinucleotides_mutations_gt"]])
  
  # signatures_foci <- do.call(rbind, by(data = mut_full, INDICES = mutations$foci, FUN = function(x) as.vector(table(x))))
  signatures_foci <- by(data = mut_full, INDICES = mutations$foci, FUN = function(x) as.vector(table(x)))
  # names(signatures_foci) <- base_type_trinuc_info[["trinucleotides_mutations"]]
  # out <- mapply(FUN = compute_likelihood_focus, counts = signatures_foci, signatures_renorm = signatures_renorm)
  
  likelihoods_foci <- t(mapply(FUN = function(counts, sigs) do.call(c, lapply(X = sigs, FUN = function(sig) prod(sig^counts))),
                               counts = signatures_foci, sigs = signatures_renorm))
  # likelihoods_foci <- likelihoods_foci/rowSums(likelihoods_foci)
  return(likelihoods_foci)
}


compute_cosine_sim <- function(signatures_renorm, mutations) {
  
  base_type_trinuc_info <- generate_bases_types_trinuc()
  
  revcomp <- data.frame(ref = as.character(reverseComplement(DNAStringSet(mutations$ref))),
                        alt = as.character(reverseComplement(DNAStringSet(mutations$alt))),
                        trinuc = as.character(reverseComplement(DNAStringSet(mutations$trinuc))),
                        stringsAsFactors = F)
  mut_full <- factor(ifelse(mutations$ref %in% c("C", "T"), 
                            paste0(mutations$ref, ">", mutations$alt, "_", mutations$trinuc),
                            paste0(revcomp$ref, ">", revcomp$alt, "_", revcomp$trinuc)),
                     levels = base_type_trinuc_info[["trinucleotides_mutations_gt"]])
  
  # signatures_foci <- do.call(rbind, by(data = mut_full, INDICES = mutations$foci, FUN = function(x) as.vector(table(x))))
  signatures_foci <- by(data = mut_full, INDICES = mutations$foci, FUN = function(x) as.vector(table(x)))
  # names(signatures_foci) <- base_type_trinuc_info[["trinucleotides_mutations"]]
  # out <- mapply(FUN = compute_likelihood_focus, counts = signatures_foci, signatures_renorm = signatures_renorm)
  
  likelihoods_foci <- t(mapply(FUN = function(counts, sigs) do.call(c, lapply(X = sigs, FUN = function(sig) sum(counts*sig)/(sqrt(sum(counts^2))*sqrt(sum(sig^2))) )),
                               counts = signatures_foci, sigs = signatures_renorm))
  # likelihoods_foci <- likelihoods_foci/rowSums(likelihoods_foci)
  return(likelihoods_foci)
}


get_cosine_sig_sims <- function(focal_mutations, signatures = pcawg_sigs, trinuc_genome = trinuc_freq, adjust_locally = T) {
  base_type_trinuc_info <- generate_bases_types_trinuc()
  
  revcomp <- data.frame(ref = as.character(reverseComplement(DNAStringSet(focal_mutations$ref))),
                        alt = as.character(reverseComplement(DNAStringSet(focal_mutations$alt))),
                        trinuc = as.character(reverseComplement(DNAStringSet(focal_mutations$trinuc))),
                        stringsAsFactors = F)
  mut_full <- factor(ifelse(focal_mutations$ref %in% c("C", "T"), 
                            paste0(focal_mutations$ref, ">", focal_mutations$alt, "_", focal_mutations$trinuc),
                            paste0(revcomp$ref, ">", revcomp$alt, "_", revcomp$trinuc)),
                     levels = base_type_trinuc_info[["trinucleotides_mutations_gt"]])
  signature_focus <- as.vector(table(mut_full))
  
  if (adjust_locally) {
  segment <- data.frame(chr = focal_mutations[1, "chromosome"], start = min(focal_mutations$pos) - 1,
                        end = max(focal_mutations$pos) + 1)
  
  # adjust signatures for local trinuc content
  signatures_renorm <- adjust_signatures_focus(locus = segment, signatures = signatures, trinuc_genome = trinuc_genome)
  } else {
    signatures_renorm <- signatures
  }
  
  # get local signatures via penalized non-negative least squares fit
  # active_signatures <- get_locus_signatures_pnnls(mutations = focal_mutations, signatures_renorm = signatures_renorm, n = n)
  cosine_sim_focus <- do.call(c, lapply(X = signatures_renorm, FUN = function(sig, counts) sum(counts*sig)/(sqrt(sum(counts^2))*sqrt(sum(sig^2))), counts = signature_focus ))
  # weighamds <- do.call(c, lapply(X = signatures_renorm, FUN = function(sig, counts) sum(sig*(1-sign(sig*counts)))/sum(sig) + sum(counts*(1-sign(counts*sig)))/sum(counts), counts = signature_focus ))
  
  dmulti <- do.call(c, lapply(X = signatures_renorm, FUN = function(sig, counts) dmultinom(x = counts, prob = sig, log = T), counts = signature_focus ))

  # dwcosine <- cosine_sim_focus/(weighamds^2 + cosine_sim_focus)  
  # names(cosine_sim_focus) <- names(pcawg_sigs)
  # return(dwcosine)
  return(list(cos = cosine_sim_focus, dmulti = dmulti))
}


get_signature_priors_typed <- function(sigsinsamplesfile = SIGSINSAMPLESFILE, histology_all = histology_all, pseudocounts = T, add_apobec = F) {
  sigsinsamples <- read.delim(file = gzfile(sigsinsamplesfile), as.is = T)
  if (add_apobec) {
    sigsinsamples$Signature.apobec <- ifelse(sigsinsamples$Signature.2 > 0 | sigsinsamples$Signature.13 > 0, 1, 0)
  }
  signaturenames <- grep(pattern = "SBS", x = colnames(sigsinsamples), value = T)
  
  if (pseudocounts) {
    pscounts <- rep(0, length(signaturenames))
    names(pscounts) <- signaturenames
    pscounts[c("SBS2", "SBS13")] <- 1
  }
  
  sigsinsamples_aug <- merge(x = histology_all[ , c("tumor_wgs_aliquot_id", "histology_abbreviation")], y = sigsinsamples,
                             by.x = "tumor_wgs_aliquot_id", by.y = "Sample.Name")

  priors <- by(data = sigsinsamples_aug[ , signaturenames],
               INDICES = sigsinsamples_aug$histology_abbreviation, FUN = function(x) (colSums(x > 0) + pscounts) / (sum(x > 0) + sum(pscounts)))
  
  return(priors)
}


get_signature_priors_sample <- function(sigsinsamplesfile = SIGSINSAMPLESFILE, pseudocounts = T, add_apobec = F) {
  sigsinsamples <- read.delim(file = gzfile(sigsinsamplesfile), as.is = T)
  if (add_apobec) {
    sigsinsamples$Signature.apobec <- ifelse(sigsinsamples$Signature.2 > 0 | sigsinsamples$Signature.13 > 0, 1, 0)
  }
  signaturenames <- grep(pattern = "SBS", x = colnames(sigsinsamples), value = T)
  
  if (pseudocounts) {
    sigsinsamples[sigsinsamples$Signature.2 == 0 , "SBS2"] <- 1
    sigsinsamples[sigsinsamples$Signature.13 == 0 , "SBS13"] <- 1
  }
  
  priors <- by(data = sigsinsamples[ , signaturenames],
               INDICES = sigsinsamples$Sample.Name, FUN = function(x) as.vector((x > 0) / sum(x > 0)) )
  
  return(priors)
}



## scoring by non-neg least-squares weighting


get_signatures_pnnls <- function(mutations, signatures_renorm, n = 1) {
  # library(penalized)
  
  base_type_trinuc_info <- generate_bases_types_trinuc()
  
  revcomp <- data.frame(ref = as.character(reverseComplement(DNAStringSet(mutations$ref))),
                        alt = as.character(reverseComplement(DNAStringSet(mutations$alt))),
                        trinuc = as.character(reverseComplement(DNAStringSet(mutations$trinuc))),
                        stringsAsFactors = F)
  mut_full <- factor(ifelse(mutations$ref %in% c("C", "T"), 
                            paste0(mutations$ref, ">", mutations$alt, "_", mutations$trinuc),
                            paste0(revcomp$ref, ">", revcomp$alt, "_", revcomp$trinuc)),
                     levels = base_type_trinuc_info[["trinucleotides_mutations_gt"]])
  
  focal_counts <- by(data = mut_full, INDICES = mutations$foci, FUN = function(x) as.vector(table(x)))
  
  focal_sigs <- mapply(FUN = function(cnts, sigs) {
    get_nonzero_coef_signatures(penalized(response = cnts/sum(cnts), penalized = t(do.call(rbind, sigs)),
                                          steps = 25, maxiter = 100, positive = T, unpenalized = ~0, trace = F), n = 3)
  },
  cnts = focal_counts, sigs = signatures_renorm)
  
  return(focal_sig)
}


adjust_signatures_focus <- function(locus, signatures = pcawg_sigs, trinuc_genome = trinuc_freq) {
  base_type_trinuc_info <- generate_bases_types_trinuc()
  
  trinuc_locus <- get_trinuc_normalisation_factors(regions = locus, overall = F)
  correction_factors <- trinuc_locus / trinuc_genome
  correction_factors <- correction_factors[, base_type_trinuc_info[["trinucleotides"]]]
  signatures_renorm <- lapply(X = signatures, FUN = function(x) x*correction_factors / (sum(x*correction_factors)))
  return(signatures_renorm)
}


get_locus_signatures_pnnls <- function(mutations, signatures_renorm, n = 1) {
  base_type_trinuc_info <- generate_bases_types_trinuc()
  
  revcomp <- data.frame(ref = as.character(reverseComplement(DNAStringSet(mutations$ref))),
                        alt = as.character(reverseComplement(DNAStringSet(mutations$alt))),
                        trinuc = as.character(reverseComplement(DNAStringSet(mutations$trinuc))),
                        stringsAsFactors = F)
  mut_full <- factor(ifelse(mutations$ref %in% c("C", "T"), 
                            paste0(mutations$ref, ">", mutations$alt, "_", mutations$trinuc),
                            paste0(revcomp$ref, ">", revcomp$alt, "_", revcomp$trinuc)),
                     levels = base_type_trinuc_info[["trinucleotides_mutations_gt"]])
  
  focal_counts <- as.vector(table(mut_full))
  focal_sigs <- tryCatch( 
    {
      penalized_fit <- penalized(response = focal_counts/sum(focal_counts), penalized = t(do.call(rbind, signatures_renorm)),
                steps = 25, maxiter = 100, positive = T, unpenalized = ~0, trace = F)
      get_nonzero_coef_signatures(penfit_list = penalized_fit, n = n)
    }, warning = function(warn) {
      print(paste("WARNING: sample ", SAMPLE, warn))
      return(list(signatures = NA, lambda = NA))
    }, error = function(err) {
      print(paste("ERROR: sample ", SAMPLE, err))
      return(list(signatures = NA, lambda = NA))
    }, finally = {}
  )
  # penalized_fit <- penalized(response = focal_counts/sum(focal_counts), penalized = t(do.call(rbind, signatures_renorm)),
  #                            steps = 25, maxiter = 100, positive = T, unpenalized = ~0, trace = F)

  if (anyNA(focal_sigs[["signatures"]])) {
    focal_sigs[["varexpl"]] <- NA 
    # focal_sigs[["coef_tstat"]] <- NA 
  } else if (length(focal_sigs[["signatures"]]) == 1) {
    lm.out <- lm(formula = focal_counts/sum(focal_counts) ~ signatures_renorm[[focal_sigs[["signatures"]] ]] - 1)
    af <- anova(lm.out)
    focal_sigs[["varexpl"]] <- ( af$`Sum Sq`/sum(af$`Sum Sq`)*100 )[1]
    # focal_sigs[["coef_tstat"]] <- coefficients(summary(lm.out))[4] 
  } else {
    nn.fit <- penalized(response = focal_counts/sum(focal_counts), penalized = do.call(cbind, signatures_renorm[focal_sigs[["signatures"]] ]),
                        lambda1 = 0, lambda2 = 0, maxiter = 100, positive = T, unpenalized = ~0, trace = F)
    focal_sigs[["varexpl"]] <- (1 - sum(residuals(nn.fit)^2)/sum((focal_counts/sum(focal_counts))^2))*100
    # focal_sigs[["coef_tstat"]] <- NA 
  }

  return(focal_sigs)
}


get_nonzero_coef_signatures <- function(penfit_list, n = 1) {
  coefs <- lapply(X = penfit_list, FUN = function(x) names(sort(coefficients(x), decreasing = T)))
  lambda1 <- lapply(X = penfit_list, FUN = function(x) x@lambda1)
  coef_idx <- which(do.call(c, lapply(X = coefs, FUN = function(x) length(x))) - n >= 0)[1]
  if (is.na(coef_idx)) {
    coef_idx <- which.min(abs(do.call(c, lapply(X = coefs, FUN = function(x) length(x))) - n))
  }
  return(list(signatures = coefs[[coef_idx]][1:n], lambda = lambda1[[coef_idx]]))
}



get_focal_signature <- function(focal_mutations, signatures = pcawg_sigs, trinuc_genome = trinuc_freq, n = 1) {
  # focal_mutations <- alldata[["pcf"]][alldata[["pcf"]]$foci == 6 , ]
  # signatures <- pcawg_sigs
  # trinuc_genome <- trinuc_freq
  # 
  segment <- data.frame(chr = focal_mutations[1, "chromosome"], start = min(focal_mutations$pos) - 1,
                        end = max(focal_mutations$pos) + 1)

  # adjust signatures for local trinuc content
  signatures_renorm <- adjust_signatures_focus(locus = segment, signatures = pcawg_sigs)
  # get local signatures via penalized non-negative least squares fit
  active_signatures <- get_locus_signatures_pnnls(mutations = focal_mutations, signatures_renorm = signatures_renorm, n = n)
  return(active_signatures)
}


is_stranded <- function(focal_mutations, trinuc_genome = trinuc_freq, bsgenome = genome) {
  # focal_mutations <- alldata[["pcf"]][alldata[["pcf"]]$foci == 2, ]
  # bsgenome = genome
  
  focal_mutations$ref <- factor(x = focal_mutations$ref, levels = c("A","C", "T", "G"))

  segment <- data.frame(chr = focal_mutations[1, "chromosome"], start = min(focal_mutations$pos),
                        end = max(focal_mutations$pos))
  
  sequences <- getSeq(x = bsgenome, names = segment$chr, start = segment$start, end = segment$end, strand = "+")
  nuc_counts_all <- oligonucleotideFrequency(x = sequences, width = 1, as.prob = F)
  nuc_ratio <- sum(nuc_counts_all[c("C", "T")]) / sum(nuc_counts_all)
  obs_counts <- table(focal_mutations$ref)
  pval <- binom.test(x = sum(obs_counts[c("C", "T")]), n = sum(obs_counts), p = nuc_ratio, alternative = "two.sided")$p.value
  return(pval)
}


get_prob_streak <- function(focal_mutations, trinuc_genome = trinuc_freq, bsgenome = genome) {
  # focal_mutations <- alldata[["pcf"]][alldata[["pcf"]]$foci == 2, ]
  # bsgenome = genome
  ntrials <- nrow(focal_mutations)
  
  refrle <- rle(focal_mutations$ref)
  streaklength <- max(refrle$lengths)
  streakbase <- refrle$values[which.max(refrle$lengths)]
  
  segment <- data.frame(chr = focal_mutations[1, "chromosome"], start = min(focal_mutations$pos),
                        end = max(focal_mutations$pos))
  sequences <- getSeq(x = bsgenome, names = segment$chr, start = segment$start, end = segment$end, strand = "+")
  pbase <- oligonucleotideFrequency(x = sequences, width = 1, as.prob = T)[[streakbase]]

  pval <- get_prob_streak_helper(ntrials = ntrials, streaklength = streaklength, psuccess = pbase)
  return(pval)
}


get_prob_streak_helper <- function(ntrials, streaklength, psuccess, saved = NULL) {
  if (is.null(saved)) 
    saved <- new.env(hash=TRUE)
  
  ID <- paste(ntrials, streaklength, psuccess, sep = "_")
  
  if (exists(x = ID, where = saved)) {
    return(saved[[ID]])
  } else {
    if (streaklength > ntrials || ntrials <= 0 ) {
      result <- 0
    } else {
      result <- exp(streaklength*log(psuccess))
      for (firsttail in 1:(streaklength+1) ) {
        pr <- get_prob_streak_helper(ntrials = ntrials - firsttail, streaklength = streaklength, psuccess = psuccess, saved = saved)
        result <- result + exp((firsttail-1)*log(psuccess) + log(1-psuccess) + log(pr))
      }
    }
    saved[[ID]] <- result
  }
  return(result)
}



get_focal_prob <- function(all_mutations, focal_mutations) {
  all_mutations$chr <- factor(x = all_mutations$chr, levels = c(1:22, "X", "Y"))
  all_mutations <- all_mutations[order(all_mutations$chr, all_mutations$pos), ]
  
  # chrs_with_multiple_muts <- names(which(by(data = all_mutations, INDICES = all_mutations$chr, FUN = function(x) nrow(x) >= 2)))
  # mutations <- mutations[mutations$chr %in% chrs_with_multiple_muts, ]

  imd <- do.call(c, by(data = all_mutations, INDICES = all_mutations$chr, FUN = function(x) diff(x$pos)))
  
  if (nrow(focal_mutations) > 1) {
    focal_imd_max <- max(diff(focal_mutations$pos))
    phat <- sum(imd <= focal_imd_max) / length(imd)
    focal_prob <- phat ^ nrow(focal_mutations)
  } else {
    focal_prob <- 1
  } 
  return(focal_prob)
}


get_phasing_graph_stats <- function(kat_phased, focus) {
  
  # kat_phased <- subset(x = phasing, subset = Chr == unique(focus$chromosome) & Pos1 %in% focus$pos & Pos2 %in% focus$pos)
  # if (nrow(kat_phased) == 0)
  #   return(c(no_phased_muts = 0, no_subclonal_muts = 0, no_antiphased_muts = 0))

  muts_on_same_allele <- kat_phased[kat_phased$phasing %in% c("phased", "clone-subclone", "subclone-clone") & kat_phased$Num_Mut_Mut > 1, ]
  no_phased_muts <- length(unique(c(muts_on_same_allele$Pos1, muts_on_same_allele$Pos2)))
  
  subclonal_muts <- kat_phased[kat_phased$phasing %in% c("clone-subclone", "subclone-clone") & (kat_phased$Num_WT_Mut > 1 | kat_phased$Num_WT_Mut > 1), ]
  no_subclonal_muts <- length(unique(ifelse(subclonal_muts$phasing == "clone-subclone", subclonal_muts$Pos2, subclonal_muts$Pos1)))
  
  muts_on_diff_allele <- kat_phased[kat_phased$phasing == "anti-phased" & kat_phased$Num_WT_Mut > 1 & kat_phased$Num_Mut_WT > 1, ]
  no_antiphased_muts <- length(unique(c(muts_on_diff_allele$Pos1, muts_on_diff_allele$Pos2)))
  
  return(c(no_phased_muts = no_phased_muts, no_subclonal_muts = no_subclonal_muts, no_antiphased_muts = no_antiphased_muts))
}


plot_phasing_graph <- function(kat_phased, focus, outdir = OUTDIR, sample_id = SAMPLE) {
  require(igraph)
 
  kat_vertices <- focus[, !colnames(focus) == "chromosome"]
  # kat_edges <- subset(x = phasing, subset = Chr == unique(focus$chromosome) & Pos1 %in% kat_vertices$pos & Pos2 %in% kat_vertices$pos)
  # if (nrow(kat_edges) == 0)
  #   return(NULL)
  
  # filtering
  kat_edges <- subset(x = kat_phased, subset = (phasing == "phased" & Num_Mut_Mut > 1) |
                                              (phasing %in% c("clone-subclone", "subclone-clone") & Num_Mut_Mut > 1 & (Num_WT_Mut > 1 | Num_WT_Mut > 1)) |
                                              (phasing == "anti-phased" & Num_WT_Mut > 1 & Num_Mut_WT > 1))
  reverse_edges <- kat_edges[kat_edges$phasing %in% c("phased", "subclone-clone"), 
                             c("Pos2", "Pos1", "Ref2", "Var2", "Ref1", "Var1", "Num_WT_WT", "Num_Mut_Mut", "Num_WT_Mut", "Num_Mut_WT", "phasing")]
  
  reverse_edges$phasing <- ifelse(reverse_edges$phasing == "subclone-clone", "clone-subclone", reverse_edges$phasing)
  colnames(reverse_edges) <- c("Pos1", "Pos2", "Ref1", "Var1", "Ref2", "Var2", "Num_WT_WT", "Num_Mut_Mut", "Num_Mut_WT", "Num_WT_Mut", "phasing")
  kat_edges <- rbind(kat_edges[!kat_edges$phasing == "subclone-clone", 
                               c("Pos1", "Pos2", "Ref1", "Var1", "Ref2", "Var2", "Num_WT_WT", "Num_Mut_Mut", "Num_Mut_WT", "Num_WT_Mut", "phasing")], reverse_edges)

  g <- graph_from_data_frame(d = kat_edges, directed=T, vertices = kat_vertices)
  E(g)$color <- ifelse(kat_edges$phasing == "phased", "green", ifelse(kat_edges$phasing == "clone-subclone", "orange", "red"))
  # V(g)$size <- kat_vertices$mut.count / (kat_vertices$mut.count + kat_vertices$WT.count) * 10
  
  png(filename = file.path(outdir, paste0(sample_id, "_graph_focus_chr", focus[1, "chromosome"], "_", focus[1, "pos"], "-", focus[nrow(focus), "pos"], ".png")))
  plot(g, layout = layout_nicely(graph = g), edge.arrow.size = .5)
  dev.off()
  return(NULL)
}


is_deaminated <- function(focal_mutations, trinuc_genome = trinuc_freq, bsgenome = genome) {

  focal_mutations$ref <- factor(x = focal_mutations$ref, levels = c("A", "C", "T", "G"))
  
  segment <- data.frame(chr = focal_mutations[1, "chromosome"], start = min(focal_mutations$pos),
                        end = max(focal_mutations$pos))
  
  sequences <- getSeq(x = bsgenome, names = segment$chr, start = segment$start, end = segment$end, strand = "+")
  nuc_counts_all <- oligonucleotideFrequency(x = sequences, width = 1, as.prob = F)
  nuc_ratio <- sum(nuc_counts_all[c("C", "G")]) / sum(nuc_counts_all)
  obs_counts <- table(focal_mutations$ref)
  pval <- binom.test(x = sum(obs_counts[c("C", "G")]), n = sum(obs_counts), p = nuc_ratio, alternative = "greater")$p.value
  return(pval)
}


has_AID_footprint <- function(focal_mutations, trinuc_genome = trinuc_freq, bsgenome = genome) {
  # count number of WRCY | RGYW mutated vs in segment
  focal_mutations$ref <- factor(x = focal_mutations$ref, levels = c("A", "C", "T", "G"))
  mutcounts <- sum(table(focal_mutations$ref)[c("C", "G")])
  if (mutcounts == 0) return(1)
  
  pentanuc <- get_trinuc_context(mutations = focal_mutations, bsgenome = bsgenome, size = 2)
  if (nrow(focal_mutations) < 2) {
    aidcounts <- sum(countPattern(pattern = "WRCYN", subject = pentanuc, fixed = F) + countPattern(pattern = "NRGYW", subject = pentanuc, fixed = F))
  } else {
    aidcounts <- sum(vcountPattern(pattern = "WRCYN", subject = pentanuc, fixed = F) + vcountPattern(pattern = "NRGYW", subject = pentanuc, fixed = F))
  }
  
  segment <- data.frame(chr = focal_mutations[1, "chromosome"], start = min(focal_mutations$pos),
                        end = max(focal_mutations$pos))
  sequences <- getSeq(x = bsgenome, names = segment$chr, start = segment$start, end = segment$end, strand = "+")
  pattern_counts_all <- sum(countPattern(pattern = "WRCY", subject = sequences, fixed = F) + countPattern(pattern = "RGYW", subject = sequences, fixed = F))
  gc_counts_all <- sum(oligonucleotideFrequency(x = sequences, width = 1, as.prob = F)[c("C", "G")])

  pval <- binom.test(x = aidcounts, n = mutcounts, p =  pattern_counts_all/gc_counts_all, alternative = "greater")$p.value
  return(pval)
}
