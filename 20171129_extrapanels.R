### try to take a look at local kataegis clusters
library(ggplot2)
str(katres_jmb_repsamples)

## sv enrichment
p1 <- ggplot(data = katres_jmb_repsamples, mapping = aes(x = log10(sv_dist+1), fill = sigsummary)) + geom_density(alpha = .3)
p1


#
tempresdf <- katresults_repsamples
# tempresdf <- katres_jmb_repsamples
tempresdf$chr <- factor(tempresdf$chr, levels = c(1:22,"X"))
tempresdf <- tempresdf[order(tempresdf$chr, tempresdf$start), ]
plot(2:nrow(tempresdf), log10(diff(tempresdf$start)))
selectv <- c(diff(tempresdf$start) < 1000, F) | c(F, diff(tempresdf$start) < 1000)
selectrle <- rle(selectv)
selectrle$values[selectrle$lengths < 3] <- F
selectdf <- tempresdf[inverse.rle(selectrle),]

# tempresdf2 <- tempresdf
tempresdf <- tempresdf2
# tempresdf <- tempresdf[tempresdf$sigsummary == "C[T>N]T",]

p2 <- ggplot(data = tempresdf, mapping = aes(x = 1:nrow(tempresdf), y = c(diff(tempresdf$start), NA))) +
               # geom_point(mapping = aes(colour = (tempresdf$chr %in% c(seq(1, 22, 2), "X"))), show.legend = F) +
  geom_point(mapping = aes(colour = sigsummary), alpha =.5, show.legend = T) +
  scale_y_log10() +  theme_minimal() +
  labs(x = "", y = "interfocal distance") + scale_x_continuous(breaks = which(tempresdf$chr[-nrow(tempresdf)] != tempresdf$chr[-1]), labels = NULL)
p2 <- p2 + theme(panel.grid.minor.x = element_blank())
p2
ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/20171129_kataegis_rainfall.pdf",
       plot = p2, width = 20, height = 6)




## on the raw cpcf output
# allpcfout$chr <- factor(allpcfout$chr, levels = c(1:22,"X"))
# allpcfout$histology <- newannots[allpcfout$sample, "histology_abbreviation"]
tempresdf <- allpcfout[order(allpcfout$chr, allpcfout$start), ]
tempresdf <- tempresdf[tempresdf$is_preferred, ]
# plot(2:nrow(tempresdf), log10(diff(tempresdf$start)))
selectv <- c(diff(tempresdf$start) < 1000, F) | c(F, diff(tempresdf$start) < 1000)
selectrle <- rle(selectv)
selectrle$values[selectrle$lengths < 3] <- F
selectdf <- tempresdf[inverse.rle(selectrle),]
# tempresdf <- tempresdf[grep(x = tempresdf$histology, pattern = "Lymph"), ]
# tempresdf <- tempresdf[tempresdf$sigsummary == "C[T>N]T",]

hist(log10(selectdf[grepl(pattern = "Lymph", x = selectdf$histology), "sv_dist"]+1 ))

cumsumdf <- as.data.frame(seqinfo(genome))
cumsumdf$cumseqlengths <- c(0,cumsum(as.numeric(cumsumdf$seqlengths))[-nrow(cumsumdf)])
tempresdf$cumstart <- tempresdf$start + cumsumdf[tempresdf$chr, "cumseqlengths"]
tempresdf$avpos <- do.call(c, by(data = tempresdf$cumstart, INDICES = tempresdf$chr,
              FUN = function(x) filter(x, rep(0.5, 2), sides=2)))
  
p2 <- ggplot(data = tempresdf, mapping = aes(x = avpos, y = c(diff(start), NA))) +
  # geom_point(mapping = aes(colour = (tempresdf$chr %in% c(seq(1, 22, 2), "X"))), show.legend = F) +
  geom_point(mapping = aes(colour = histology), alpha =.5, show.legend = F) +
  scale_y_log10() +  theme_minimal() +
  labs(x = "", y = "interfocal distance") + scale_x_continuous(breaks = cumsumdf[c(1:22, "X", "Y"), "cumseqlengths"], labels = NULL)
p2 <- p2 + theme(panel.grid.minor.x = element_blank())
p2
ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/20171219_kataegis_rainfall_unfiltered_repsamples.pdf",
       plot = p2, width = 20, height = 6)


##### just need to run from here
## rainfall formatting for Circos
tempresdf_circos <- allpcfout[allpcfout$is_preferred, ]
tempresdf_circos$avpos <- round(rowMeans(x = tempresdf_circos[, c("start", "end")]))
tempresdf_circos <- tempresdf_circos[order(tempresdf_circos$chr, tempresdf_circos$avpos), ]

# plot(2:nrow(tempresdf_circos), log10(diff(tempresdf_circos$avpos)))

tempresdf_circos$ifd <- c(tempresdf_circos[2:nrow(tempresdf_circos), "avpos"] - tempresdf_circos[1:(nrow(tempresdf_circos)-1), "avpos"], NA)
tempresdf_circos$ifd <- c(ifelse(tempresdf_circos[2:nrow(tempresdf_circos), "chr"] == tempresdf_circos[1:(nrow(tempresdf_circos)-1), "chr"], tempresdf_circos$ifd[-nrow(tempresdf_circos)], NA), NA)
tempresdf_circos$ifd <- log10(ifelse(tempresdf_circos$ifd > 0, tempresdf_circos$ifd, 1))
# tempresdf_circos$ifd <- log10(ifelse(tempresdf_circos$ifd <= 0, 1, tempresdf_circos$ifd))
# tempresdf_circos$avpos <- round(unlist(by(data = tempresdf_circos, INDICES = tempresdf_circos$chr, FUN = function(x) c(x[2:(nrow(x)), "start"] + x[1:(nrow(x)-1), "end"], NA)/2)))
# tempresdf_circos$avposplus <- tempresdf_circos$avpos + 1
tempresdf_circos$chr <- paste0("hs", tempresdf_circos$chr)

# unique(tempresdf_circos$tumour.subtype)[which(!unique(tempresdf_circos$tumour.subtype) %in% pcawg.colour.palette(scheme = "tumour.subtype", return.scheme = T)$levels)]
# pcawg.colour.palette(x = tempresdf_circos$tumour.subtype[1:5], scheme = "tumour.subtype")
# mycolours <- data.frame(t(col2rgb(pcawg.colour.palette(scheme = "tumour.subtype", return.scheme = T)$colours)))
# rownames(mycolours) <- sub(pattern = ".", replacement = "-", x = pcawg.colour.palette(scheme = "tumour.subtype", return.scheme = T)$levels, fixed = T)
# write.table(x = mycolours, file = "kataegis_recurrence/mycolours.conf", sep = ",", row.names = T, col.names = F, quote = F)

tempresdf_circos$colour <- paste0("fill_color=", tolower(tempresdf_circos$histology))

write.table(tempresdf_circos[!is.na(tempresdf_circos$ifd), c("chr", "avpos", "avpos", "ifd", "colour")], file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/kataegis_recurrence/data", "hypermut.foci.txt"), sep = " ", row.names = F, col.names = F, quote = F)



#### select Lymph-tumour regions with excess (>3 in repres samples) kataegis events
# tempresdf2 <- allpcfout[order(allpcfout$chr, allpcfout$start), ]
# tempresdf2 <- tempresdf2[tempresdf2$sample %in% newannots$samplename[newannots$is_preferred], ]
tempresdf_lymphonly <- tempresdf_circos[grepl(pattern = "Lymph-", x = tempresdf_circos$histology), ]
tempresdf_gr <- GRanges(seqnames = sub(pattern = "hs", replacement = "", x = tempresdf_lymphonly$chr), ranges = IRanges(start = tempresdf_lymphonly$start, end = tempresdf_lymphonly$end), mcols = DataFrame(sample = tempresdf_lymphonly$sample))
tempresdf_gr_reduced <- reduce(x = tempresdf_gr, min.gapwidth = 1001)
katoverlaps <- findOverlaps(query = tempresdf_gr_reduced, subject = tempresdf_gr)
katoverlapcounts <- countOverlaps(query = tempresdf_gr_reduced, subject = tempresdf_gr)
# which(katoverlapcounts >= 3)
metafoci <- tempresdf_gr_reduced[which(katoverlapcounts >= 3)]
mcols(metafoci) <- DataFrame(overlapcount <- katoverlapcounts[which(katoverlapcounts >= 3)])
metafoci_df <- as.data.frame(metafoci)
metafoci_events <- tempresdf_lymphonly[tempresdf_gr %over% metafoci, ]
metafoci_events_gr <- GRanges(seqnames = metafoci_events$chr, ranges = IRanges(start = metafoci_events$start, end = metafoci_events$end), mcols = DataFrame(sample = metafoci_events$sample))


## collect genes within/near regions
library(biomaRt)
## Get genes inside regions:
ensembl37 <- useMart(host = "grch37.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl37)

chromregions <- paste(seqnames(metafoci),
                      start(metafoci)-1e6,
                      end(metafoci)+1e6, sep = ":")

out <- getBM(attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name"),
             filters = c("chromosomal_region", "biotype", "with_hgnc"), values = list(chromregions, c("protein_coding", "lincRNA", "miRNA", "antisense", "IG_C_gene", "IG_D_gene", "IG_V_gene", "IG_J_gene"), T), mart = ensembl)

# out$external_gene_name <- ifelse(grepl(pattern = "^IGH.*", x = out$external_gene_name, perl = T), "IGH", 
                                 # ifelse(grepl(pattern = "^IGL.*", x = out$external_gene_name, perl = T), "IGL",
                                        # ifelse(grepl(pattern = "^IGK.*", x = out$external_gene_name, perl = T), "IGK", out$external_gene_name)))
out$external_gene_name <- ifelse(out$chromosome_name == "14" & out$start_position >= 106032614 & out$end_position <= 107288051, "IGH",
                                 ifelse(out$chromosome_name == "2" & out$start_position >= 89156674 & out$end_position <= 89630436, "IGK",
                                        ifelse(out$chromosome_name == "22" & out$start_position >= 22380474 & out$end_position <= 23265203, "IGL" ,out$external_gene_name)))
# out$histology <- "Lymph"

# unique_drivers <- setdiff(unique(kataegis_drivers$gene), out$external_gene_name)
# out2_otherdrivers <- getBM(attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name"),
             # filters = c("hgnc_symbol"), values = list(unique_drivers), mart = ensembl)

out_gr <- GRanges(seqnames = out$chromosome_name, ranges = IRanges(start = out$start_position, end = out$end_position))
# metafoci_df$target <- out_report[nearest(x = metafoci, subject = out_gr, ignore.strand = T), "external_gene_name"]
metafoci_nearest_gene_hits <- nearest(x = metafoci, subject = out_gr, ignore.strand = T, select = "all")
out <- out[subjectHits(metafoci_nearest_gene_hits), ]

out <- out[!duplicated(x = out$external_gene_name), ]
# out <- rbind(out, out2_otherdrivers)
out$chromosome_name <- paste0("hs", out$chromosome_name)

cancergenecensus <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/kataegis_recurrence/Census_allWed Jan  3 14_56_17 2018.tsv", as.is = T)
pcawgcancergenes <- read.delim(file = "/srv/shared/vanloo/ICGC_driver/TableS1_compendium_mutational_drivers_20180110.txt", as.is = T)
# out[out$external_gene_name %in% c(cancergenecensus$Gene.Symbol, pcawgcancergenes$gene), ]
out$colour = ifelse(out$external_gene_name %in% c(cancergenecensus$Gene.Symbol, pcawgcancergenes$Gene), "color=red,link_color=red", "color=black,link_color=black")

write.table(out, file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/kataegis_recurrence/data", "hit.genes.txt"), sep = " ", row.names = F, col.names = F, quote = F)



### try to make circos of kataegis rainfall with SVs in 

# CLEANHISTOLOGYFILE <- "/srv/shared/vanloo/data/vanloo/jdemeul/ICGC/kataegis/merged_histology_v6_release_v1.4.txt"
SVDIR <- "/srv/shared/vanloo/ICGC-structural-variants/"
CNDIR <- "/srv/data/vanloo/sdentro/ICGC/battenberg_rerun_consensusBP20160906/"
# CIRCOSDIR <- "/srv/data/vanloo/jdemeul/ICGC/quaid_circos/"
# setwd(CIRCOSDIR)


# histology <- read.delim(file = CLEANHISTOLOGYFILE, as.is = T)
lymphsamples <- unique(metafoci_events$sample)


# for (SAMPLENAME in lymphsamples) {

read_selected_svs <- function(sampleid, svdir, metafoci) {
  svfile <- grep(pattern = ".sv.bedpe.gz$", x = list.files(path = SVDIR, pattern = paste0(sampleid), full.names = T), value = T)
  svs <- read.delim(file = svfile, as.is = T)
  
  ## subset those that fall near (within 10kb) hypermutation metafocus
  svs_gr1 <- GRanges(seqnames = svs$chrom1, ranges = IRanges(start = svs$start1, end = svs$end1))
  svs_gr2 <- GRanges(seqnames = svs$chrom2, ranges = IRanges(start = svs$start2, end = svs$end2))
  # metafoci_events_insample <- metafoci_events[mcols(metafoci_events_gr)$mcols.sample == sampleid]
  # svs <- svs[overlapsAny(query = svs_gr1, subject = metafoci_events_insample, maxgap = 1e4) |
               # overlapsAny(query = svs_gr2, subject = metafoci_events_insample, maxgap = 1e4), ]
  svs <- svs[overlapsAny(query = svs_gr1, subject = metafoci, maxgap = 1e4) |
               overlapsAny(query = svs_gr2, subject = metafoci, maxgap = 1e4), ]
  
  if (nrow(svs) > 0) {
    svs$formatting <- ifelse(svs$svclass == "DUP", "color=red_a3",
                             ifelse(svs$svclass == "DEL", "color=blue_a3",
                                    ifelse(svs$svclass == "t2tINV", "color=green_a3",
                                           ifelse(svs$svclass == "h2hINV", "color=vvdyellow_a3",
                                                  "color=black_a2"))))
    
    ### Reformat structural variants bedpe file to a Circos input file
    ## Reformat the data and output to two tet files: rearr.ribs.txt and rearr.ribs.txt containing ribbons and links resp
    svs$chrom1 <- paste0("hs" ,svs$chrom1)
    svs$chrom2 <- paste0("hs" ,svs$chrom2)
    df.out.links <-  svs[,c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "formatting")]
  } else {
    df.out.links <- data.frame(chrom1 = character(), start1 = integer(), end1 = integer(),
                               chrom2 = character(), start2 = integer(), end2 = integer(), formatting = character())
  }
  return(df.out.links)
}
  # write.table(df.out.links, file.path(CIRCOSDIR, "data", "rearr.links.txt"), sep = " ", row.names = F, col.names = F, quote = F)
  
# temp <- read_selected_svs(sampleid = lymphsamples[1], svdir = SVDIR)
alllymphsvs <- lapply(X = lymphsamples, FUN = read_selected_svs, svdir = SVDIR, metafoci = metafoci)
# sort(sapply(X = alllymphsvs, FUN = nrow))
alllymphsvsdf <- do.call(rbind, alllymphsvs)

write.table(alllymphsvsdf, file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/kataegis_recurrence/data", "rearr.links.txt"), sep = " ", row.names = F, col.names = F, quote = F)





################# 
###### code to assign foci to nearest gene and report table for marker:
chromregions2 <- paste(seqnames(metafoci),
                      start(metafoci)-1e6,
                      end(metafoci)+1e6, sep = ":")

out_report <- getBM(attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name"),
             filters = c("chromosomal_region", "biotype", "with_hgnc"), values = list(chromregions2, c("protein_coding", "lincRNA", "miRNA", "antisense", "IG_C_gene", "IG_D_gene", "IG_V_gene", "IG_J_gene"), T), mart = ensembl)

# out$external_gene_name <- ifelse(grepl(pattern = "^IGH.*", x = out$external_gene_name, perl = T), "IGH", 
# ifelse(grepl(pattern = "^IGL.*", x = out$external_gene_name, perl = T), "IGL",
# ifelse(grepl(pattern = "^IGK.*", x = out$external_gene_name, perl = T), "IGK", out$external_gene_name)))
out_report$external_gene_name <- ifelse(out_report$chromosome_name == "14" & out_report$start_position >= 106032614 & out_report$end_position <= 107288051, "IGH",
                                 ifelse(out_report$chromosome_name == "2" & out_report$start_position >= 89156674 & out_report$end_position <= 89630436, "IGK",
                                        ifelse(out_report$chromosome_name == "22" & out_report$start_position >= 22380474 & out_report$end_position <= 23265203, "IGL" ,out_report$external_gene_name)))

out_gr <- GRanges(seqnames = out_report$chromosome_name, ranges = IRanges(start = out_report$start_position, end = out_report$end_position))
# metafoci_df$target <- out_report[nearest(x = metafoci, subject = out_gr, ignore.strand = T), "external_gene_name"]
metafoci_nearest_gene_hits <- nearest(x = metafoci, subject = out_gr, ignore.strand = T, select = "all")
metafoci_df$genes <- c(by(data = out_report[subjectHits(metafoci_nearest_gene_hits), "external_gene_name"], INDICES = queryHits(metafoci_nearest_gene_hits),
                           FUN = function(x) paste0(unique(x), collapse = ",")))
# metafoci_df <- metafoci_df[subjectHits(metafoci_nearest_gene_hits)]

colnames(metafoci_df) <- c("chr", "hotspot_start", "hotspot_end", "width", "strand", "No_foci", "genes")
write.table(metafoci_df[, c("chr", "hotspot_start", "hotspot_end", "No_foci", "genes")], file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/drivers_and_recurrence/20180319_focalhypermut_annotated.txt"), sep = "\t", row.names = F, col.names = T, quote = F)

  
  