## identifying commonalities in top 5% samples with extensive kataegis
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

BASEOUT <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/final_rerun_annotmuts/"
TOPHITSDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/tophits_rainfall/"

source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/pcawg.colour.palette.R")



katresults_jamboree <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_Kataegis_calls_JD_allcolumns.txt", as.is = T)
katresults_repsamples <- katresults_jamboree[katresults_jamboree$is_preferred & katresults_jamboree$is_punctuated, ]

### kataegis "intensity" plot - Ludmill style
sinaplotdf2 <- katresults_repsamples %>% group_by(sample, histology) %>% summarise(no_foci = length(total), no_sv_foci = sum(sv_assoc), no_apo_foci = sum(signature == "APO"))
### cervix adeno is basically only contributor to other
# sinaplotdf2 <- sinaplotdf2[!sinaplotdf2$histology %in% infreq_types, ]
sinaplotdf2$histology <- factor(x = sinaplotdf2$histology, levels = sort(unique(sinaplotdf2$histology)))
# sinaplotdf2 <- sinaplotdf2[order(sinaplotdf2$histology, sinaplotdf2$no_foci, decreasing = F), ]
# sinaplotdf2$rnk <- unlist(by(data = sinaplotdf2, INDICES = sinaplotdf2$histology, FUN = function(x) seq(0, .7, .7/(nrow(x)-1))))
# sinaplotdf2[is.nan(sinaplotdf2$rnk), "rnk"] <- .35
# sinaplotdf2$rnk <- sinaplotdf2$rnk + as.numeric(sinaplotdf2$histology)
# 
# write.table(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_kataegis_tumortypes_intensities.txt", x = sinaplotdf2, quote = F, sep = "\t", row.names = F)
# 
# meddf <- do.call(rbind, by(data = sinaplotdf2, INDICES = sinaplotdf2$histology, FUN = function(x) data.frame(xmin = min(x$rnk), xmax = max(x$rnk), y = median(x = x$no_foci))))
# meddf$histology <- rownames(meddf)
# meddf["CNS-PiloAstro", c("xmin", "xmax")] <- meddf["CNS-PiloAstro", c("xmin", "xmax")] + c(-.35,.35)
# 
# 
# cvect <- pcawg.colour.palette(x = tolower(sub(pattern = "-", replacement = ".", x = levels(sinaplotdf2$histology))), scheme = "tumour.subtype")
# names(cvect) <- levels(sinaplotdf2$histology)
# cvect[c("Kidney-RCC-Clear", "Other", "Kidney-RCC-Pap")] <- c('#FF4500', '#DDCDCD', '#FF4500')
# 
# psum4 <- ggplot(data = sinaplotdf2, mapping = aes(x = rnk, y = no_foci)) + geom_point(mapping = aes(fill = histology), stroke = .25, shape = 21, colour = "black", alpha = .75, size = .75) + scale_y_log10()
# # psum4 <- ggplot(data = sinaplotdf2, mapping = aes(x = rnk, y = no_foci)) + geom_point(colour = "black", alpha = .75, size = .25) + scale_y_log10()
# psum4 <- psum4 + geom_segment(data = meddf, mapping = aes(x = xmin, xend = xmax, y = y, yend = y), alpha = .5)
# psum4 <- psum4 + scale_fill_manual(values = cvect, guide = F) + labs(y = "# foci")
# psum4 <- psum4 + theme_minimal() + theme(axis.text.x = element_blank(), text = element_text(size = 6), axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
# psum4 <- psum4 + annotation_logticks(sides = "l", scaled = T, colour = "grey")
# psum4


# katfocidf <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_kataegis_tumortypes_intensities.txt", as.is = T)

# nrow(katfocidf)*.05
# quantile(katfocidf$no_foci, c(.05, .95))


####### plot of different classes of katsamples at top
katfoci_tophits <- sinaplotdf2[sinaplotdf2$no_foci > 30, ]
katfoci_tophits$frac_sv <- katfoci_tophits$no_sv_foci/katfoci_tophits$no_foci
katfoci_tophits$frac_apo <- katfoci_tophits$no_apo_foci/katfoci_tophits$no_foci

katfoci_tophits$class <- ifelse(katfoci_tophits$frac_apo < .5, "non-APOBEC", ifelse(katfoci_tophits$frac_sv <= .2, "Replication", ifelse(katfoci_tophits$frac_sv <= .45, "Mix", "SV-associated")))
katfoci_tophits$class <- factor(x = katfoci_tophits$class, levels = c("non-APOBEC", "SV-associated", "Mix", "Replication"))

sort(table(katfoci_tophits$histology))

# for (sampleid in katfoci_tophits$sample) {
#   file.copy(from = file.path(BASEOUT, sampleid, paste0(sampleid, "_rainfall.png")), to = TOPHITSDIR)
# }


cvect <- pcawg.colour.palette(x = tolower(sub(pattern = "-", replacement = ".", x = levels(sinaplotdf2$histology))), scheme = "tumour.subtype")
names(cvect) <- levels(sinaplotdf2$histology)
# cvect[c("Kidney-RCC-Clear", "Kidney-RCC-Pap")] <- c('#FF4500','#FF4500')
cvect[c("Skin-Melanoma-Acral", "Skin-Melanoma-Cut", "Kidney-RCC-Clear", "Kidney-RCC-Pap")] <- c("#000000", "#000000", '#FF4500', '#FF4500')


p1 <- ggplot(data = katfoci_tophits, mapping = aes(x = class, fill = histology)) + geom_bar() + scale_fill_manual(values = cvect, guide = F) + labs(y = "# samples")
p1 <- p1 + theme_minimal() + coord_flip() + labs(x = "", y = "# samples") + theme(panel.grid.major.x = element_blank())
p1

ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190204_kataegis_tophits_classification.pdf",
       plot = p1, width = 91.5, height = 40, units = "mm")


######## plot of distribution of SV distances
katresults_jamboree$signature <- factor(x = katresults_jamboree$signature, levels = c("APO", "ALT", "CTT", "POLH", "uncertain"))


p1 <- ggplot(data = katresults_jamboree[katresults_jamboree$is_preferred & katresults_jamboree$is_punctuated, ], mapping = aes(x = ifelse(sv_dist == 0, 100, sv_dist), fill = signature)) + geom_density(alpha = .8) + scale_x_log10(breaks = 10^(0:8), labels = 0:8) + annotation_logticks(base = 10, sides = "b")
# p1 <- ggplot(data = katresults_jamboree[katresults_jamboree$is_preferred & katresults_jamboree$is_punctuated, ], mapping = aes(x = ifelse(sv_dist == 0, 100, sv_dist), fill = signature)) + geom_histogram(alpha = .8, bins = 100) + scale_x_log10(breaks = 10^(0:8), labels = 0:8) + annotation_logticks(base = 10, sides = "bl")
p1 <- p1 + scale_colour_manual(values = c("#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4", "#E8E8E8", "#FFFFFF")) + scale_fill_manual(values = c("#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4", "#E8E8E8", "#FFFFFF")) +  theme_minimal() + geom_vline(mapping = aes(xintercept = 1000), alpha = .5, linetype = "dashed")  
p1 <- p1 + labs(x = "log10(distance to nearest SV (bp) )") + theme(panel.grid.minor = element_blank(), panel.grid = element_blank())
p1

ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190204_kataegis_sv_distance.pdf",
       plot = p1, width = 122, height = 91.5, units = "mm")




# 
# p1 <- ggplot(data = katfoci_tophits, mapping = aes(x = no_apo_foci/no_foci, y= no_sv_foci/no_foci, colour = histology)) + geom_point()
# p1
# 
# p1 <- ggplot(data = katfoci_tophits[katfoci_tophits$histology != "Lymph-BNHL",], mapping = aes(x= no_sv_foci/no_foci, fill = histology)) + geom_histogram(bins = 50)
# p1
# 
# p1 <- ggplot(data = katfoci_tophits[grep(pattern = "Breast", x = katfoci_tophits$histology),], mapping = aes(x= no_sv_foci/no_foci, fill = histology)) + geom_histogram(bins = 50)
# p1
# 
# p1 <- ggplot(data = sinaplotdf2[sinaplotdf2$no_foci > 10, ], mapping = aes(x= no_sv_foci/no_foci, fill = histology)) + geom_histogram(bins = 50)
# p1
# 
# 
# table(sinaplotdf2[sinaplotdf2$no_foci > 10 & sinaplotdf2$no_sv_foci/sinaplotdf2$no_foci < .25, "histology"])
# 
# 



### check bladder-TCC and Head-SCC samples for kataegis focus enrichment on lagging strand.
# Bladder_HNsamples <- katresults_repsamples[katresults_repsamples$histology %in% c("Skin-Melanoma", "SoftTissue-Liposarc", "SoftTissue-Leiomyo"), ]

# for every sample, get the annotated muts and check whether T or A in excess for each focus
get_focus_WC_annot <- function(sampleid, baseout) {
  samplemuts <- read.delim(file = file.path(baseout, sampleid, paste0(sampleid, "_all_muts_annot.txt")), as.is = T)
  samplemuts$ref <- factor(samplemuts$ref, levels = c("A", "C", "G", "T"))
  WCannots <- c(by(data = samplemuts$ref, INDICES = paste0(samplemuts$chr, ":", samplemuts$start, "-", samplemuts$end), FUN = function(x) {y <- table(x); if (y["C"] > y["G"]) "W" else if (y["G"] > y["C"]) "C" else NA} ))
  WCannotsdf <- data.frame(sampleid = sampleid, focus = names(WCannots), WC = WCannots, stringsAsFactors = F)
  return(WCannotsdf)
}

WCannots <- do.call(rbind, lapply(X = unique(katresults_jamboree$sample), FUN = get_focus_WC_annot, baseout = BASEOUT))
katresults_jamboree$WCannot <- WCannots[match(x = paste0(katresults_jamboree$chr, ":", katresults_jamboree$start, "-", katresults_jamboree$end), table = WCannots$focus), "WC"]
katresults_jamboree[is.na(katresults_jamboree$signature) | katresults_jamboree$signature != "APO", "WCannot"] <- NA




katresults_repsamples <- katresults_jamboree[katresults_jamboree$is_preferred & katresults_jamboree$is_punctuated & katresults_jamboree$signature == "APO", ]

#now match back to full dataframe
katresults_repsamples_gr <- GRanges(seqnames = katresults_repsamples$chr, ranges = IRanges(start = katresults_repsamples$start, end = katresults_repsamples$end))

## read in Okazaki seq RFD bed files
GMrep1 <- suppressWarnings(unlist(import.bedGraph(con = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/data/GM_rep1.bed")))
GMrep2 <- suppressWarnings(unlist(import.bedGraph(con = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/data/GM_rep2.bed")))
HeLarep1 <- suppressWarnings(unlist(import.bedGraph(con = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/data/Hela_rep1.bed")))
HeLarep2 <- suppressWarnings(unlist(import.bedGraph(con = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/data/Hela_rep2.bed")))

mcols(GMrep1)$avescore <- rowMeans(x = cbind(GMrep1$score,GMrep2$score))
mcols(HeLarep1)$avescore <- rowMeans(x = cbind(HeLarep1$score,HeLarep2$score))

seqlevelsStyle(GMrep1) <- "Ensembl"
seqlevelsStyle(HeLarep1) <- "Ensembl"


nearestOKidxs <- nearest(x = katresults_repsamples_gr, subject = GMrep1, select = "arbitrary")
katresults_repsamples$RFD_GM <- mcols(GMrep1)[nearestOKidxs, "avescore"]
katresults_repsamples[which(katresults_repsamples$RFD_GM < -1), "RFD_GM"] <- -1
katresults_repsamples$RFD_Hela <- mcols(HeLarep1)[nearestOKidxs, "avescore"]
katresults_repsamples[which(katresults_repsamples$RFD_Hela < -1), "RFD_Hela"] <- -1
# katresults_repsamples$sv_assoc_strict <- katresults_repsamples$sv_dist <= 1000


replitim <- import.bw(con = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/data/wgEncodeUwRepliSeqHelas3WaveSignalRep1.bigWig")
seqlevelsStyle(replitim) <- "Ensembl"
nearestRTidxs <- nearest(x = katresults_repsamples_gr, subject = replitim, select = "arbitrary")
katresults_repsamples$reptiming <- mcols(replitim)[nearestRTidxs, "score"]


GMrep1check <- GMrep1
mcols(GMrep1check)$reptiming <- mcols(replitim)[nearest(x = GMrep1check, subject = replitim, select = "arbitrary"), "score"]
mcols(GMrep1check)$avescore[mcols(GMrep1check)$avescore < -1] <- -1

# note: if deamination on Watson, then Okazaki == Crick (+)
katresults_repsamples$Lagging_GM <- ifelse(!is.na(katresults_repsamples$WCannot), ifelse(katresults_repsamples$WCannot == "W", katresults_repsamples$RFD_GM, -katresults_repsamples$RFD_GM), NA)
katresults_repsamples$Lagging_Hela <- ifelse(!is.na(katresults_repsamples$WCannot), ifelse(katresults_repsamples$WCannot == "W", katresults_repsamples$RFD_Hela, -katresults_repsamples$RFD_Hela), NA)


katresults_repsamples_intop <- katresults_repsamples[katresults_repsamples$sample %in% katfoci_tophits$sample,]

# p1 <- ggplot(data = katresults_repsamples[katresults_repsamples$sample %in% katfoci_tophits$sample,], mapping = aes(x = Lagging_Hela, fill = sv_assoc)) + geom_histogram(binwidth = .05, position = "identity", alpha = .3)
# p1
# 
# tabsvassoc <- table(lagging = katresults_repsamples_intop[katresults_repsamples_intop$sv_assoc, "Lagging_Hela"] >= 0, timing = katresults_repsamples_intop[katresults_repsamples_intop$sv_assoc, "reptiming"] >= 40)
# chisq.test(tabsvassoc)
# tabnonsvassoc <- table(lagging = katresults_repsamples_intop[!katresults_repsamples_intop$sv_assoc, "Lagging_Hela"] >= 0, timing = katresults_repsamples_intop[!katresults_repsamples_intop$sv_assoc, "reptiming"] >= 40)
# chisq.test(tabnonsvassoc)

# tabsvassoc
# tabnonsvassoc
# 
# (tabnonsvassoc/sum(tabnonsvassoc))/(tabsvassoc/sum(tabsvassoc))



# tabnonsvassoc <- table(sv_assoc = katresults_repsamples_intop$sv_assoc, lagtiming = (katresults_repsamples_intop$Lagging_Hela >= 0 | katresults_repsamples_intop$reptiming >= 40))
# tabnonsvassoc
# oddsratioWald.proc(n00 = 218, n01 = 2260, n10 = 249, n11 = 1163, alpha = .05)
# chisq.test(tabnonsvassoc)
# 
# tabnonsvassoc <- table(sv_assoc = katresults_repsamples_intop$sv_assoc, lag = katresults_repsamples_intop$Lagging_Hela >= 0)
# tabnonsvassoc
wilcox.test(x = katresults_repsamples_intop[!katresults_repsamples_intop$sv_assoc, "Lagging_Hela"],
            y = katresults_repsamples_intop[katresults_repsamples_intop$sv_assoc, "Lagging_Hela"], alternative = "greater")
wilcox.test(x = katresults_repsamples_intop[!katresults_repsamples_intop$sv_assoc, "reptiming"],
            y = katresults_repsamples_intop[katresults_repsamples_intop$sv_assoc, "reptiming"], alternative = "greater")



# by(data = katresults_repsamples_intop, INDICES = katresults_repsamples_intop$sv_assoc, FUN = function(x) mean(x$reptiming, na.rm = T))


# [grepl(pattern = "Bladder|Head|Skin|SoftTissue", x = katresults_repsamples$histology),]
p1 <- ggplot(data = katresults_repsamples[katresults_repsamples$sample %in% katfoci_tophits$sample,], mapping = aes(x = reptiming, fill = sv_assoc, color = sv_assoc)) + geom_histogram(binwidth = 1, position = "identity", alpha = .6, show.legend = F) + geom_density(mapping = aes(y = ..count..), alpha = .1)
p1 <- p1 + scale_fill_manual(values = c('TRUE' = "#fbb4ae", 'FALSE' ="#FDD2CE")) + scale_color_manual(values = c('TRUE' = "#fbb4ae", 'FALSE' ="#FDD2CE")) + theme_minimal() + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text.x = element_blank())
p1 <- p1 + labs(x = "<<- late replicating --- early replication ->> ")
p1

ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190204_kataegis_reptiming.pdf",
       plot = p1, width = 122, height = 91.5, units = "mm")


p1 <- ggplot(data = katresults_repsamples[katresults_repsamples$sample %in% katfoci_tophits$sample,], mapping = aes(x = Lagging_Hela, fill = sv_assoc, color = sv_assoc)) + geom_histogram(binwidth = .025, position = "identity", alpha = .6, show.legend = F) + geom_density(mapping = aes(y = ..count../40), alpha = .1)
p1 <- p1 + scale_fill_manual(values = c('TRUE' = "#fbb4ae", 'FALSE' ="#FDD2CE")) + scale_color_manual(values = c('TRUE' = "#fbb4ae", 'FALSE' ="#FDD2CE")) + theme_minimal() + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text.x = element_blank())
p1 <- p1 + labs(x = "<<- leading strand --- lagging strand ->> ")
p1

ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190204_kataegis_strandbias.pdf",
       plot = p1, width = 122, height = 91.5, units = "mm")


p1 <- ggplot(data = katresults_repsamples[katresults_repsamples$sample %in% katfoci_tophits$sample,], mapping = aes(x = reptiming, y = Lagging_Hela, colour = sv_assoc)) + geom_point()
p1

# hist(mcols(replitim)[, "score"])

# 
# table(sv_assoc = katresults_repsamples$sv_assoc_strict,
#                  lagging = katresults_repsamples$Lagging_Hela > 0)
# 
# 
oddsratioWald.proc <- function(n00, n01, n10, n11, alpha = 0.05){
  #
  #  Compute the odds ratio between two binary variables, x and y,
  #  as defined by the four numbers nij:
  #
  #    n00 = number of cases where x = 0 and y = 0
  #    n01 = number of cases where x = 0 and y = 1
  #    n10 = number of cases where x = 1 and y = 0
  #    n11 = number of cases where x = 1 and y = 1
  #
  OR <- (n00 * n11)/(n01 * n10)
  #
  #  Compute the Wald confidence intervals:
  #
  siglog <- sqrt((1/n00) + (1/n01) + (1/n10) + (1/n11))
  zalph <- qnorm(1 - alpha/2)
  logOR <- log(OR)
  loglo <- logOR - zalph * siglog
  loghi <- logOR + zalph * siglog
  #
  ORlo <- exp(loglo)
  ORhi <- exp(loghi)
  #
  oframe <- data.frame(LowerCI = ORlo, OR = OR, UpperCI = ORhi, alpha = alpha)
  oframe
}
# 
# 
# prop.test(x = 2536, n = 2536+3405, p = .5)
# prop.test(x = 2338, n = 2338+2493, p = .5)
# oddsratioWald.proc(n00 = 2536, n01 = 3405, n10 = 2338, n11 = 2493, alpha = .05)
# 
# 
# 
# # ori's
# 
# Oris <- HeLarep1
# Oris$timing <- mcols(replitim)[nearest(x = Oris, subject = replitim, select = "arbitrary"), "score"]
# Oris <- reduce(Oris[which(abs(Oris$avescore) < .05 & Oris$timing > 70)])
# # Oris
# nearestOriidxs <- distanceToNearest(x = katresults_repsamples_gr, subject = Oris, select = "arbitrary")
# katresults_repsamples[queryHits(nearestOriidxs), "ori_dist"] <- mcols(nearestOriidxs)$distance
# 
# p1 <- ggplot(data = katresults_repsamples, mapping = aes(x = ori_dist+1, fill = sv_assoc_strict)) + geom_histogram(position = "identity", bins = 100, alpha = .5) + scale_x_log10() + annotation_logticks(base = 10, sides = "b")
# p1
# 
# 


### check exon enrichment
# # library(BSgenome.Hsapiens.UCSC.hg19)
# 
hstxdb <- makeTxDbFromGFF(file = "/srv/shared/vanloo/pipeline-files/human/references/annotation/GENCODE/gencode.v23lift37.annotation.gtf", organism = "Homo sapiens")
# seqlevels(hstxdb) <- sub(pattern = "chr", replacement = "", x = seqlevels(seqinfo(hstxdb)))
hsexondb <- exons(x = hstxdb, columns = c("gene_id"))
seqlevelsStyle(hsexondb) <- "Ensembl"

katresults_repsamples$overexon <- overlapsAny(query = katresults_repsamples_gr, subject = hsexondb)
table(over_exon = katresults_repsamples[katresults_repsamples$sample %in% katfoci_tophits$sample,"overexon"],
      sv_assoc = katresults_repsamples[katresults_repsamples$sample %in% katfoci_tophits$sample,"sv_assoc"])
chisq.test(table(over_exon = katresults_repsamples[katresults_repsamples$sample %in% katfoci_tophits$sample,"overexon"],
                 sv_assoc = katresults_repsamples[katresults_repsamples$sample %in% katfoci_tophits$sample,"sv_assoc"]))

sum(width(reduce(hsexondb)))
sum(seqlengths(seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5)))
1313/(1313+4675)
839/(839+4015)

oddsratioWald.proc(n00 = 1853, n01 = 1164, n10 = 635, n11 = 254)


# ###### additional exon enrich checks
# 
# snvfiles <- list.files(path = "/srv/shared/vanloo/ICGC_snv_mnv/final_consensus_12oct_passonly/", pattern = ".consensus.20160830.somatic.snv_mnv.vcf.gz$", full.names = T, recursive = T)
# 
# katsnvfiles <- list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/final_rerun_annotmuts/", pattern = "_all_muts.txt$", full.names = T, recursive = T)
# samples <- gsub(pattern = "_all_muts.txt", replacement = "", x = basename(katsnvfiles))
# 
# katsnvs <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_Kataegis_SNVs_JD.txt", as.is = T)
# 
# library(VariantAnnotation)
# library(BSgenome.Hsapiens.1000genomes.hs37d5)
# 
# checkexonoverlaps <- function(sampleid, katsnvfiles, snvfiles, hsexondb) {
#   print(sampleid)
#   katsnvs <- read.delim(file = grep(pattern = sampleid, x = katsnvfiles, value = T), as.is = T)
#   katsnvs <- GRanges(seqnames = katsnvs$chromosome, ranges = IRanges(start = katsnvs$pos, end = katsnvs$pos), seqinfo = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
#   
#   allsnvs <- readVcfAsVRanges(x = grep(pattern = sampleid, x = snvfiles, value = T), genome = seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5))
#   
#   allsnvs$inkat <- overlapsAny(query = allsnvs, subject = katsnvs, type = "equal")
#   allsnvs$exonic <- overlapsAny(query = allsnvs, subject = hsexondb, type = "within")
#   return(table(inkat = allsnvs$inkat, exonic = allsnvs$exonic))
# }
# # sampleid <- samples[1]
# 
# # debug(checkexonoverlaps)
# exonenrichout <- mclapply(X = samples, FUN = checkexonoverlaps, katsnvfiles = katsnvfiles, snvfiles = snvfiles, hsexondb = hsexondb, mc.preschedule = T, mc.cores = 16)
# exonenrichoutpchi <- sapply(X = exonenrichout, FUN = function(x) chisq.test(x)$p.value)
# exonenrichpool <- colSums(do.call(rbind, lapply(X = exonenrichout, FUN = 'c')))
# exonenrichsep <- as.data.frame(do.call(rbind, lapply(X = exonenrichout, FUN = 'c')))
# exonenrichsep$pval <- exonenrichoutpchi
# exonenrichsep$sample <- samples
# 
# exonenrichsep <- merge(x = exonenrichsep, y = sinaplotdf2, all = T)
# View(exonenrichsep[which(exonenrichsep$no_foci>=10), ])

# 
# 
# #### GSEA style test for checking associations with heavy non-SV assoc kataegis
# # create gmx file with columns are driver genes and rows samples containing that driver
# drivers <- read.delim(file = "/srv/shared/vanloo/ICGC_driver/TableS2_driver_point_mutations_annotation_20180110.txt", as.is = T)
# drivers$gene <- sub(pattern = "::.*$", replacement = "", x = sub(pattern = "^.*?::.*?::", replacement = "", x = drivers$gene_id))
# 
# head(drivers)
# samplesbydriv <- by(data = drivers, INDICES = drivers$gene, FUN = function(x) unique(x$sample))
# # emptyvec <- rep("",  max(lengths(samplesbydriv)))
# 
# samplesbydriv <- do.call(cbind, by(data = drivers, INDICES = drivers$gene, FUN = function(x) {y <- unique(x$sample); c(y, rep("", 900-length(y)))}))
# write.table(x = samplesbydriv, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/data/samplesbydriver.gmx", quote = F, sep = "\t", row.names = F)
# 
# outdf <- sinaplotdf2[grepl(pattern = "Bladder", x = sinaplotdf2$histology), ]
# outdf$logfoci <- log(outdf$no_foci - outdf$no_sv_foci + 1) - mean(log(outdf$no_foci - outdf$no_sv_foci + 1))
# write.table(x = outdf[order(outdf$logfoci, decreasing = T), c("sample", "logfoci")], file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/focirankedlist.rnk", quote = F, sep = "\t", row.names = F)
