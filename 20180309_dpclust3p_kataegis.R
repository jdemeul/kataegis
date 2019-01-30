#' Identify Kataegis events in SNV mutation data
#'
#' @param samplename Samplename used when writing output files
#' @param dpInfile A DP input file
#' @param outdir Directory where output will be written
#' @param gamma_param Gamma parameter to be used for segmentation
#' @param kmin Kmin parameter to be used for segmentation
#' @param kataegis.threshold Intermutation distance, if NA will be set to: (a) 900000 > SNVs: 100, (b) 500000 > SNVs: 250, (c) 100000 > SNVs, (d) otherwise 1000 (Default NA)
#' @param minMuts Minimum number of mutations within kataegis.threshold to call Kataegis
#' @param logScale Transform intermutation distance to logscale
#' @param makePlots Make a figure for each Kataegis event
#' @param removeFromDPin Remove SNVs identified as part of a Kataegis event from the DP input file
#' @author dw9
#' @export
identifyKataegis <- function(samplename, snvs, outdir = ".", gamma_param = 25, kmin = 2, pstreak = 0.01, minmutsrange = 4:6, maxthresh = 1e3, bsgenome = genome){
  kataegis.thresholds <- get_kataegis_threshold(snvs = snvs, pstreak = pstreak, minmutsrange = minmutsrange, maxthresh = maxthresh)
  print(paste0(samplename, ": ml ", nrow(snvs), " / thresh ", kataegis.thresholds$thresh, " / minmuts ", kataegis.thresholds$minmuts))
  
  snvs_gr <- sort(GRanges(seqnames = snvs$chr, ranges = IRanges(start = snvs$pos, end = snvs$pos), mcols = snvs[, c("alt", "ref")], seqinfo = seqinfo(bsgenome)))
  
  is_dinuc <- unlist(by(data = start(snvs_gr), INDICES = seqnames(snvs_gr), FUN = function(x) c(F, diff(x) == 1)))

  katloci <- do.call(rbind, by(data = snvs_gr[!is_dinuc], INDICES = seqnames(snvs_gr[!is_dinuc]), FUN = run_pcf_chr,
                                    kmin = kmin, gamma_param = gamma_param, dthresh = kataegis.thresholds$thresh, nthresh = kataegis.thresholds$minmuts))
  
  if (nrow(katloci) == 0) return(NULL)
  
  katloci_gr <- GRanges(seqnames = katloci$chr, ranges = IRanges(start = katloci$start, end = katloci$end), seqinfo = seqinfo(bsgenome))
  
  snv_hits <- findOverlaps(query = snvs_gr, subject = katloci_gr, maxgap = kataegis.thresholds$thresh)
  kat_snv <- snvs_gr[queryHits(snv_hits)]
  mcols(kat_snv)$focus <- subjectHits(snv_hits)
  
  # # is_katlocus$values <- ifelse(is_katlocus$lengths >= kataegis.thresholds$minmuts, is_katlocus$values, F)
  # kat_snv <- snvs_gr[is_katlocus$is_kat]
  # mcols(kat_snv)$focus_id <- is_katlocus$focus_id[is_katlocus$is_kat]
  kat_snv <- as.data.frame(kat_snv)[ , c("seqnames", "start", "mcols.ref", "mcols.alt", "focus")]
  colnames(kat_snv) <- c("chr", "pos", "ref", "alt", "focus")
  
  print(paste0("     nfoci ", max(kat_snv$focus), " / nvar ", nrow(kat_snv)))
  
  write.table(kat_snv, file.path(outdir, paste0(samplename, "_kataegis_annotated.txt")), sep = "\t", quote = F, row.names = F)
  return(NULL)
}


run_pcf_chr <- function(snvs, kmin, gamma_param, dthresh, nthresh) {
  imds <- diff(start(snvs))

  # run pcf, omitting first mut in dinucs
  sdev <- getMad(imds, k=25)
  res <- exactPcf(imds, kmin, gamma_param*sdev,T)
  
  is_katlocus <- rle(res$yhat <= dthresh)
  is_katlocus$values <- ifelse(is_katlocus$lengths >= nthresh - 1, is_katlocus$values, F)
  
  if (sum(is_katlocus$values) == 0) {
    return(data.frame(chr = character(), start = integer(), end = integer()))
  }
  
  is_katlocus <- inverse.rle(is_katlocus)
  start_regions = start(snvs)[ (c(is_katlocus, F) & !c(F, is_katlocus)) ]
  end_regions = end(snvs)[ (!c(is_katlocus, F) & c(F, is_katlocus)) ]
  
  outdf <- data.frame(chr = rep(seqnames(snvs)[1], length(start_regions)),
                      start = start_regions, end = end_regions)
  return(outdf)
}


get_kataegis_threshold <- function(snvs, pstreak = 0.001, minmutsrange = 5:6, maxthresh = 1e3) {
  medimd <- median(unlist(by(data = snvs$pos, INDICES = snvs$chr, FUN = function(x) diff(sort(x)))))
  mutload <- nrow(snvs)
  
  exprate <- log(2)/medimd
  p <- (pstreak/mutload)^(1/(minmutsrange-1))
  d <- -log(1-p)/(exprate)
  
  d[d > maxthresh] <- maxthresh
  
  thresh <- d[which(d >= maxthresh)][1]
  minmuts <- minmutsrange[which(d >= maxthresh)][1]
  if (is.na(thresh)) {
    thresh <- d[length(minmutsrange)]
    minmuts <- max(minmutsrange)
  }
  
  return(list(thresh = thresh, minmuts = minmuts))
}