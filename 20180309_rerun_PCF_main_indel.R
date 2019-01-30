## rerun PCF on SNVs
library(doParallel)
source(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/20180309_dpclust3p_kataegis.R")
source(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/20180309_dpclust3p_fastPCF.R")

# STATICS
CLEANHISTOLOGYFILE <- "/srv/shared/vanloo/ICGC_annotations/summary_table_combined_annotations_v2.txt"
SNVDIR <- "/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/"
OUTDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/20180524_pcf_rerun_indel"
NTHREADS <- 4


## functions
read_snv_file <- function(sample_id, snv_mnvdir) {
  snvfile <- list.files(path = snv_mnvdir, pattern = paste0(sample_id, ".*somatic.indel.vcf.gz$"), full.names = T, recursive = T)
  allmuts <- tryCatch(
    {
      mutsvcf <- read.delim(file = gzfile(snvfile), header = F, sep = "\t", as.is = T, comment.char = "#")
      mutsvcf[grepl(pattern = "t_alt_count=", x = mutsvcf$V8), c(1,2,4,5)]
    }, error = function(err) {
      print(paste("ERROR: sample ", sample_id, err))
      return(data.frame(chr = character(), pos = integer(), ref = character(), alt = character()))
    }
  )
  colnames(allmuts) <- c("chr", "pos", "ref", "alt")
  return(allmuts)
}



## end functions

histology_all <- read.delim(file = CLEANHISTOLOGYFILE, as.is = T)

# for every sample, run

clp <- makeCluster(NTHREADS)
registerDoParallel(clp)

foreach(i=1:nrow(histology_all), .packages = c("BSgenome", "BSgenome.Hsapiens.1000genomes.hs37d5", "GenomicRanges"), .verbose = T) %dopar% {
# for (i in 1:nrow(histology_all)) {
  # sample_id <- "00c27940-c623-11e3-bf01-24c6515278c0"
  sample_id <- histology_all$samplename[i]
  # debug(read_snv_file)
  allsnvs <- read_snv_file(sample_id = sample_id, snv_mnvdir = SNVDIR)
  # debug(identifyKataegis)
  # debug(run_pcf_chr)
  identifyKataegis(samplename = sample_id, snvs = allsnvs, outdir = OUTDIR, maxthresh = 1e4, pstreak = .01, minmutsrange = 3:5, bsgenome = BSgenome.Hsapiens.1000genomes.hs37d5)
}

stopCluster(clp)
