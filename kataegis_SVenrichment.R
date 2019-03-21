# View(katres_jmb_repsamples)

library(ggplot2)
library(GenomicRanges)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(parallel)

SVDIR <- "/srv/shared/vanloo/ICGC-structural-variants/"
SNVDIR <- "/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/"
CLEANHISTOLOGYFILE <- "/srv/shared/vanloo/ICGC_annotations/summary_table_combined_annotations_v2.txt"

source(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/kataegis_functions.R", local = T)

histology_all <- read_histology(histologyfile = CLEANHISTOLOGYFILE)


katres_jmb <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_Kataegis_calls_JD_allcolumns.txt", as.is = T)
katres_jmb_repsamples <- katres_jmb[katres_jmb$is_preferred, ]

## visualize distance to SVs for kat types
p1 <- ggplot(data = katres_jmb_repsamples, mapping = aes(x = sv_dist+1, fill = sigsum_SV)) + scale_x_log10() + geom_histogram(alpha = .4)
p1 <- p1 + theme_minimal()
p1

qplot(x = end-start, data = katres_jmb_repsamples, geom = "histogram")


# SAMPLE <- "075fc96d-6742-4ef3-9369-482592ad3a2f"

# debug(analyse_enrichment)
# debug(get_sv_enrichment_numbers)
enrichcounts <- mclapply(X = unique(katres_jmb_repsamples$sample), FUN = analyse_enrichment, svdir = SVDIR, snv_mnvdir = SNVDIR, all_katresults = katres_jmb_repsamples, genome = BSgenome.Hsapiens.1000genomes.hs37d5, mc.cores = 1, mc.preschedule = T)
# enrichcounts <- mclapply(X = unique(katres_jmb_repsamples$sample), FUN = analyse_enrichment, svdir = SVDIR, snv_mnvdir = SNVDIR, all_katresults = katres_jmb_repsamples, genome = BSgenome.Hsapiens.1000genomes.hs37d5, mc.cores = 12)

enrichcounts <- as.data.frame(do.call(rbind, enrichcounts))
enrichcounts$sample_id <- unique(katres_jmb_repsamples$sample) 
enrichcounts$histology <- histology_all[match(x = enrichcounts$sample_id, table = histology_all$samplename), "histology_abbreviation"]

p1 <- ggplot(data = enrichcounts, mapping = aes(x = histology, y = (n_APO_sv/(n_APO_nosv+n_APO_sv)) / (n_nonkat_sv/(n_nonkat_nosv+n_nonkat_sv)))) +
  geom_violin(mapping = aes(colour = histology), show.legend = F) + geom_jitter(mapping = aes(colour = histology), alpha = .5, show.legend = F) +
  theme_minimal() + theme(axis.text.x = element_text(angle = -90)) + scale_y_log10() + labs(x = "", y = "APOBEC SV association enrichment")
p1
ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190131_APO_SV_association.pdf", plot = p1, width = 20, height = 6)

save(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190131_kat_SV_enrichments.RData", enrichcounts)

#### functions

# read all required PCAWG data for a sample
read_svenrich_data <- function(sample_id, svdir, snv_mnvdir, all_katresults = katres_jmb_repsamples) {
  
  svfile <- list.files(path = svdir, pattern = paste0(sample_id, "[0-9A-Za-z_.]*", ".somatic.sv.bedpe.gz$"), full.names = T)
  svbreakpoints <- tryCatch(
    {
      svs <- read.delim(file = gzfile(svfile), header = T, sep = "\t", as.is = T)
      colnames(svs) <- c("chrom.1", "start.1", "end.1", "chrom.2", "start.2", "end.2", "id",
                         "pe_support", "strand.1", "strand.2", "svclas", "svmethod")
      if (nrow(svs) == 0) {
        svbreakpoints <- data.frame(id = character(), bp_side = numeric(), chrom = character(), start = integer())
      } else {
        svbreakpoints <- reshape(svs[ , c(1,2,4,5,7)], direction = "long", varying = 1:4, timevar = "bp_side")
      }
      row.names(svbreakpoints) <- NULL
      svbreakpoints
    }, error = function(err) {
      print(paste("ERROR: sample ", sample_id, err))
      return(NULL)
    }
  )
  
  snvfile <- list.files(path = snv_mnvdir, pattern = paste0(sample_id, ".consensus.20160830.somatic.snv_mnv.vcf.gz$"), full.names = T, recursive = T)
  allmuts <- tryCatch(
    {
      allmuts <- read.delim(file = gzfile(snvfile), header = F, sep = "\t", as.is = T, comment.char = "#")[, c(1,2,4,5)]
      colnames(allmuts) <- c("chr", "pos", "ref", "alt")
      allmuts
    }, error = function(err) {
      print(paste("ERROR: sample ", sample_id, err))
      return(NULL)
    }
  )
  
  return(list(sv = svbreakpoints,
              snv = allmuts,
              katout = katres_jmb_repsamples[katres_jmb_repsamples$sample == sample_id, ]))
}



# get the enrichment numbers
get_sv_enrichment_numbers <- function(alldata, genome = genome) {
  svs_gr <- GRanges(seqnames = alldata$sv$chrom, ranges = IRanges(start = alldata$sv$start, end = alldata$sv$start), seqinfo = genome@seqinfo)
  snvs_gr <- GRanges(seqnames = alldata$snv$chr, ranges = IRanges(start = alldata$snv$pos, end = alldata$snv$pos), seqinfo = genome@seqinfo)
  katfoci_gr <- GRanges(seqnames = alldata$katout$chr, ranges = IRanges(start = alldata$katout$start, alldata$katout$end), seqinfo = genome@seqinfo,
                        mcols = factor(alldata$katout$signature, levels = c("ALT", "APO", "CTT", "POLH", "uncertain")))
  
  katoverhits <- findOverlaps(query = snvs_gr, subject = katfoci_gr)
  katmuts <- snvs_gr[queryHits(katoverhits)]
  mcols(katmuts) <- mcols(katfoci_gr[subjectHits(katoverhits)])
  
  nonkatmuts <- subsetByOverlaps(x = snvs_gr, ranges = katfoci_gr, invert = T)

  katmuts_svassoc <- factor(overlapsAny(query = katmuts, subject = svs_gr, maxgap = 1e4), levels = c(F, T))
  nonkatmuts_svassoc <- overlapsAny(query = nonkatmuts, subject = svs_gr, maxgap = 1e4)
  
  katmuts_tab <- as.vector(table(data.frame(type = mcols(katmuts), sv = katmuts_svassoc)))
  names(katmuts_tab) <- c("n_ALT_nosv", "n_APO_nosv", "n_CTT_nosv", "n_POLH_nosv", "n_unc_nosv",
                          "n_ALT_sv", "n_APO_sv", "n_CTT_sv", "n_POLH_sv", "n_unc_sv")
  
  enrichcounts <- c(n_nonkat_nosv = sum(!nonkatmuts_svassoc), n_nonkat_sv = sum(nonkatmuts_svassoc), katmuts_tab)
  return(enrichcounts)
}


# wrapper function to apply to katsamples
analyse_enrichment <- function(sample_id, svdir, snv_mnvdir, all_katresults = katres_jmb_repsamples, genome = genome) {
  alldata <- tryCatch( 
    {
      read_svenrich_data(sample_id = sample_id, svdir = SVDIR, snv_mnvdir = SNVDIR, all_katresults = all_katresults)
    }, error = function(err) {
      print(paste("ERROR: sample ", sample_id, err))
      return(NULL)
    }, finally = {}
  )
  
  if (is.null(alldata$sv) || is.null(alldata$snv)) return(as.vector(rep(NA, 12), mode = "integer"))
  
  enrichcounts <- get_sv_enrichment_numbers(alldata = alldata, genome = genome)
  return(enrichcounts)
}


