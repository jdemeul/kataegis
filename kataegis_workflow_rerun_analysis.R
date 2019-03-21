## summarize initial signature assingments:
library(ggplot2)
library(reshape2)
library(dplyr)

source(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/kataegis_functions.R", local = T)
source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/pcawg.colour.palette.R")

RESULTSBASE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/final_rerun_annotmuts/"
SIGFILE <- "/srv/shared/vanloo/ICGC_signatures/20180322_release/sigProfiler_SBS_signatures.csv"
KATSIGFILE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20180322_kataegis_signature_patterns.csv"
CLEANHISTOLOGYFILE <- "/srv/shared/vanloo/ICGC_annotations/summary_table_combined_annotations_v2.txt"

histology_all <- read_histology(histologyfile = CLEANHISTOLOGYFILE)
all_sigs <- c(load_signatures(SIGFILE, mergesigs = c(7,10,17)), load_signatures(KATSIGFILE, pad = T))

### checking cosine simil sigs
pcffiles <- list.files(path = RESULTSBASE, pattern = "_kataegis_cpcf.txt", full.names = T, recursive = T)
allpcfout <- do.call(rbind, lapply(X = pcffiles, FUN = read.delim, header = T, sep = "\t", as.is = T))
allpcfout$histology <- histology_all[match(x = allpcfout$sample, table = histology_all$samplename), "histology_abbreviation"]
allpcfout[, c("nrpcc", "is_preferred")] <- histology_all[match(allpcfout$sample, histology_all$samplename), c("nrpcc", "is_preferred")]
# allpcfout$active_sig <- factor(allpcfout$active_sig, levels = names(pcawg_sigs))
allpcfout$sig1 <- factor(allpcfout$sig1, levels = names(all_sigs), labels = names(all_sigs))
allpcfout$sig2 <- factor(allpcfout$sig2, levels = names(all_sigs), labels = names(all_sigs))
allpcfout$sig3 <- factor(allpcfout$sig3, levels = names(all_sigs), labels = names(all_sigs))

# p1 <- ggplot(data = allpcfout, mapping = aes(x = sig1, y = histology, colour = histology)) + geom_jitter(alpha = .25) + theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90))
# p1


allpcfout$p_streak_adj <- p.adjust(ifelse(allpcfout$p_streak > 1, 1, allpcfout$p_streak), method = "fdr")
allpcfout$p_nostrandbias_adj <- p.adjust(allpcfout$p_nostrandbias, method = "fdr")
allpcfout$p_aid_adj <- NA
allpcfout[allpcfout$histology %in% c("Lymph-BNHL", "Lymph-CLL"), "p_aid_adj"] <- p.adjust(allpcfout[allpcfout$histology %in% c("Lymph-BNHL", "Lymph-CLL"), "p_aid"], method = "fdr")

# p1 <- ggplot(data = allpcfout, mapping = aes(x = -log10(p_streak), y = no_phased_muts)) + geom_point()
# p1

# allpcfout_similarcosine <- allpcfout[allpcfout$sig_cos1 - allpcfout$sig_cos2 < .1, ]
# p1 <- ggplot(data = allpcfout_similarcosine, mapping = aes(x = sig_cos1_name, y = sig_cos2_name, colour = histology)) + geom_jitter(alpha = .25) + theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90))
# p1
# 

# allpcfout_spectra <- as.data.frame(do.call(rbind, lapply(X = allpcfout$mutspectrum, FUN = function(x) as.numeric(unlist(strsplit(x = x, split = ","))))))
# colnames(allpcfout_spectra) <- generate_bases_types_trinuc()[["trinucleotides_mutations"]]
# allpcfout_spectra <- data.frame(allpcfout[, c("sig1", "sig2", "sample", "histology")], allpcfout_spectra)
# pr.out <- prcomp(allpcfout_spectra[, -c(1:4)] , scale=T)
# 
# p1 <- ggplot(data = allpcfout_spectra, mapping = aes(x = sig1, y = clust)) + geom_jitter(alpha = .25) + theme(axis.text.x = element_text(angle = 90))
# p1
# # 
# pr_plotdf <- data.frame(allpcfout_spectra[, 1:4], pr.out$x)
# p1 <- ggplot(data = pr_plotdf[,], mapping = aes(x = PC1, y = PC2, colour = sig1)) + geom_jitter(alpha = .25)
# p1
# 
# p1 <- ggplot(data = pr_plotdf[pr_plotdf$sig_cos1_name %in% c("Signature.2", "Signature.13", "Signature.1", "Signature.7a", "Signature.7b", "Signature.7c", "Signature.17a", "Signature.17b", "Signature.10a", "Signature.10b"),], mapping = aes(x = PC2, colour = sig_cos1_name)) + geom_density()
# p1

# pr_plotrot <- as.data.frame(pr.out$rotation)
# pr_plotrot$muttype <- generate_bases_types_trinuc()[["trinucleotides_mutations"]]
# # pr_plotrot_melt <- melt(data = pr_plotrot, id.vars = "muttype")
# p2 <- ggplot(data = pr_plotrot, mapping = aes(x = muttype, y = PC4)) + geom_col() + theme(axis.text.x = element_text(angle = 90))
# p2

# allpcfout$upweigth_apobec <- ifelse(allpcfout$active_sig %in% factor(c("Signature.2", "Signature.13"), levels = c("Signature.2", "Signature.13")), 
# allpcfout$active_sig)
# 
# allpcfout_cosines <- apply(X = allpcfout[, c("sig_cosines", "sig_cosines_names")], MARGIN = 1, FUN = function(x) setNames(object = as.numeric(unlist(strsplit(x = x[1], split = ","))), nm = unlist(strsplit(x = x[2], split = ","))))
# allpcfout_cosines <- as.data.frame(do.call(rbind, lapply(X = allpcfout_cosines, FUN = function(x) setNames(object = x[names(pcawg_sigs)], nm = names(pcawg_sigs)))))
# allpcfout_cosines <- data.frame(allpcfout[, c("sig_cos1_name", "sig_cos2_name", "sample", "histology")], allpcfout_cosines)
# 
# p1 <- ggplot(data = allpcfout_cosines, mapping = aes(x = Signature.9, y = Signature.5, colour = sig_cos1_name)) + geom_jitter(alpha = .25) + theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90))
# p1
# 
# allpcfout_cosines[is.na(allpcfout_cosines)] <- 0
# pr.out <- prcomp(allpcfout_cosines[, -c(1:4,45:49,52)] , scale=F)
# 
# pr_plotdf <- data.frame(allpcfout_cosines[, 1:4], pr.out$x)
# p1 <- ggplot(data = pr_plotdf[,], mapping = aes(x = PC3, y = PC4, colour = pr_clust_cut)) + geom_jitter(alpha = .25)
# p1
# 
# # p1 <- ggplot(data = pr_plotdf[pr_plotdf$sig_cos1_name %in% c("Signature.2", "Signature.13", "Signature.1", "Signature.7a", "Signature.7b", "Signature.7c", "Signature.17a", "Signature.17b", "Signature.10a", "Signature.10b"),], mapping = aes(x = PC2, colour = sig_cos1_name)) + geom_density()
# # p1
# 
# pr_plotrot <- as.data.frame(pr.out$rotation)
# pr_plotrot$sig <- rownames(pr_plotrot)
# # pr_plotrot_melt <- melt(data = pr_plotrot, id.vars = "muttype")
# p2 <- ggplot(data = pr_plotrot, mapping = aes(x = sig, y = PC4)) + geom_col() + theme(axis.text.x = element_text(angle = 90))
# p2
# 
# 
# pr_clust <- hclust(dist(pr_plotdf[, c("PC1", "PC2", "PC3", "PC4")], method = "euclidean"))
# kmeans()
# pr_clust_cut <- cutree(tree = pr_clust, k = 6)


sigs_per_sample <- by(data = allpcfout$sig1, INDICES = allpcfout$sample, unique)
sigs_per_sample_df <- data.frame(sample = rep(names(sigs_per_sample), lapply(X = sigs_per_sample, FUN = length)), signature = unlist(sigs_per_sample), stringsAsFactors = F)
# histology_all[match(sigs_per_sample_df$sample, histology_all$samplename), "histology_abbreviation"]
sigs_per_sample_df$histology <- histology_all[match(sigs_per_sample_df$sample, histology_all$samplename), "histology_abbreviation"]


## summary stats/plots:
p1 <- ggplot(data = allpcfout) + geom_bar(mapping = aes(x = sig1, fill = sig1), show.legend = F)
p1 <- p1 + theme_minimal() + theme(axis.text.x = element_text(size = 6, angle = 90), legend.position = "none")
p1 <- p1 + scale_x_discrete(labels = sub(pattern = "SBS", replacement = "", x = names(all_sigs)), drop = F)
p1 <- p1 + facet_wrap(~ histology, scales = "free_y")
p1

ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_focal_signatures_overall.png",
       plot = p1, width = 25, height = 8)


# p1 <- ggplot(data = allpcfout[allpcfout$total >= 5, ]) + geom_jitter(mapping = aes(x = sig_cos1_name, y = sig_cos1, colour = sig_cos1_name), show.legend = F, shape = ".", alpha = .85)
# p1 <- p1 + geom_boxplot(mapping = aes(x = sig_cos1_name, y = sig_cos1, colour = sig_cos1_name), show.legend = F, alpha = .25)
# p1 <- p1 + theme_minimal() + theme(axis.text.x = element_text(size = 6, angle = 90), legend.position = "none")
# p1 <- p1 + scale_x_discrete(labels = sub(pattern = "Signature.", replacement = "", x = names(all_sigs)), drop = F)
# p1

p1 <- ggplot(data = allpcfout) + geom_jitter(mapping = aes(x = sig1, y = no_phased_muts/total, colour = sig1), show.legend = F, shape = ".", alpha = .85)
p1 <- p1 + geom_boxplot(mapping = aes(x = sig1, y = no_phased_muts/total, colour = sig1), show.legend = F, alpha = .25)
p1 <- p1 + theme_minimal() + theme(axis.text.x = element_text(size = 6, angle = 90), legend.position = "none")
p1 <- p1 + scale_x_discrete(labels = sub(pattern = "SBS", replacement = "", x = names(all_sigs)), drop = F)
p1

p1 <- ggplot(data = allpcfout) + geom_histogram(mapping = aes(x = no_phased_muts/total, fill = sig1), show.legend = F, binwidth = .05)
p1 <- p1 + theme_minimal() + theme(axis.text.x = element_text(size = 6, angle = 90), legend.position = "none")
p1 <- p1 + facet_wrap( ~ sig1, scales = "free_y")
p1

p1 <- ggplot(data = sigs_per_sample_df) + geom_bar(mapping = aes(x = histology, fill = histology), show.legend = F)
p1 <- p1 + theme_minimal() + theme(axis.text.x = element_text(size = 6, angle = 90), legend.position = "none")
p1 <- p1 + scale_x_discrete(drop = F)
p1 <- p1 + facet_wrap( ~ signature, scales = "free_y")
p1

p2 <- ggplot(data = sigs_per_sample_df) + geom_bar(mapping = aes(x = signature, fill = signature), show.legend = F)
p2 <- p2 + theme_minimal() + theme(axis.text.x = element_text(size = 6, angle = 90), legend.position = "none")
p2 <- p2 + scale_x_discrete(labels = sub(pattern = "SBS", replacement = "", x = names(all_sigs)), drop = F)
p2 <- p2 + facet_wrap(~ histology, scales = "free_y")
p2

ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_focal_signatures_allhistol_persample.pdf",
       plot = p2, width = 25, height = 8)


p1 <- ggplot(data = allpcfout, mapping = aes(y = -log10(p_nostrandbias_adj), x = total)) + geom_jitter(show.legend = F, shape = ".")
p1 <- p1 + geom_density2d() 
p1 <- p1 + coord_cartesian(xlim = c(0, 15), ylim = c(0, 2.5))
p1

# p1 <- ggplot(data = allpcfout, mapping = aes(x = sigs1-sig_cos2)) + geom_histogram(bins = 250)
# p1


allpcfout_clean <- allpcfout[( allpcfout$p_streak_adj <= .1 | (allpcfout$no_phased_muts >= 4 & allpcfout$no_phased_muts/allpcfout$total >= .75)) &
                               !(allpcfout$no_antiphased_muts/allpcfout$total >= .1 | allpcfout$no_subclonal_muts/allpcfout$total >= .1), ]

# write main list
write.table(x = allpcfout_clean, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_Kataegis_Results_all.txt", sep = "\t", col.names = T, row.names = F, quote = F)



####

katresults <- allpcfout_clean

### for upload to jamboree
acceptsigs <- paste0("SBS", c(2, 9, 13, "17", 28, "APO", "CTT", "ALT"))
aposigs <- paste0("SBS", c(2, 13, "APO"))
cttsigs <- paste0("SBS", c("17", 28, "CTT"))
altsigs <- "SBSALT"
polnsigs <- "SBS9"

# final assingment 
# katresults$signature <- katresults$sig1
katresults$signature <- ifelse(katresults$sig1 %in% acceptsigs,
                               ifelse(katresults$sig1 %in% aposigs, "APO", 
                                      ifelse(katresults$sig1 %in% cttsigs, "CTT",
                                             ifelse(katresults$sig1 %in% altsigs, 
                                                    ifelse(katresults$sig2 %in% aposigs, "APO", "ALT"),
                                                    ifelse(katresults$sig1 %in% polnsigs, "POLH", "should not be here")))),
                               ifelse(katresults$sig2 %in% acceptsigs,
                                      ifelse(katresults$sig2 %in% aposigs, "APO", 
                                             ifelse(katresults$sig2 %in% cttsigs, "CTT",
                                                    ifelse(katresults$sig2 %in% altsigs,
                                                           ifelse(katresults$sig3 %in% aposigs, "APO", "ALT"),
                                                           ifelse(katresults$sig2 %in% polnsigs, "POLH", "should not be here")))),
                                      ifelse(katresults$sig3 %in% aposigs, "APO", 
                                             ifelse(katresults$sig3 %in% cttsigs, "CTT",
                                                    ifelse(katresults$sig3 %in% altsigs, "ALT",
                                                           ifelse(katresults$sig3 %in% polnsigs, "POLH", "uncertain"))))))

# & !katresults$histology %in% c("Lymph-BNHL", "Lymph-CLL")
katresults$sigsum_SV <- paste0(katresults$signature, ifelse(katresults$sv_dist <= 1e3, "_SV", ""))
katresults$sigsum_SV <- factor(x = katresults$sigsum_SV, levels = c("APO_SV", "APO", "ALT_SV", "ALT", "CTT_SV", "CTT", "POLH_SV", "POLH", "uncertain_SV", "uncertain"))

katresults_repsamples <- katresults[katresults$is_preferred, ]

# plotdf[katres_jmb_repsamples$histology_abbreviation %in% infreq_types, "histology_abbreviation"] <- "other"
# 
# plotdf <- as.data.frame(do.call(rbind, by(data = katresults_repsamples, INDICES = katresults_repsamples$sample, FUN = function(x) table(x$sigsum_SV)/sum(table(x$sigsum_SV)) )))
# # applplotdf[is.nan(plotdf)] <- 0
# plotdf$sample <- rownames(plotdf)
# plotdf$histology <- histology_all[match(plotdf$sample, histology_all$samplename), "histology_abbreviation"]
# # define infreq types
# infreq_types <- names(which(table(histology_all$histology_abbreviation) < 10))
# # infreq types
# plotdf[plotdf$histology %in% infreq_types, "histology"] <- "other"
# # plotdf[is.na(plotdf$APOBEC_SV), 1:10] <- 0
# 
# plotdf_tt <- as.data.frame(do.call(rbind, by(data = plotdf[, -c(11,12)], INDICES = plotdf$histology, FUN = colSums)))
# plotdf_tt[,1:10] <- plotdf_tt[,1:10] / rowSums(plotdf_tt[,1:10])
# plotdf_tt$histology <- rownames(plotdf_tt)
# plotdf_tt_melt <- melt(data = plotdf_tt, id.vars = 11, measure.vars = 1:10)
# 
# pl2 <- ggplot(data = plotdf_tt_melt) + geom_bar(stat = "identity", mapping = aes(x = histology, y = value, fill = variable))
# pl2 <- pl2 + theme_minimal() + theme(axis.text.x = element_text(angle = 90)) + labs(x = "", y = "fraction")
# pl2 <- pl2 + scale_fill_manual(values = c("#fbb4ae", "#FDD2CE", "#b3cde3", "#D1E1EE","#ccebc5", "#E0F3DC", "#decbe4", "#EBE0EF", "#E8E8E8", "#F1F1F1"))
# pl2
# ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/20170802_kataegis_tumortypes_signatures_v2.pdf",
#        plot = pl2, width = 20, height = 6)


####FINAL FIGURE MARKER
### redo with all samples from histol
plotdf <- as.data.frame(do.call(rbind, by(data = katresults_repsamples, INDICES = katresults_repsamples$sample, FUN = function(x) table(x$sigsum_SV)/sum(table(x$sigsum_SV)) )))
plotdf$sample <- rownames(plotdf)
plotdf <- merge(x = plotdf, y = histology_all[histology_all$is_preferred, "samplename", drop = F], by.x = "sample", by.y = "samplename", all.x = F, all.y = T)
plotdf$no_kat <- is.na(plotdf$APO_SV)
plotdf$histology <- histology_all[match(plotdf$sample, histology_all$samplename), "histology_abbreviation"]
# redefined infreq types
samplecounts <- c(table(histology_all[histology_all$is_preferred, "histology_abbreviation"]))
infreq_types <- names(which(samplecounts < 10))
# infreq types
plotdf <- plotdf[!plotdf$histology %in% infreq_types, ]
# plotdf[plotdf$histology %in% infreq_types, "histology"] <- "Other"
plotdf[is.na(plotdf$APO_SV), 2:11] <- 0

plotdf_tt <- as.data.frame(do.call(rbind, by(data = plotdf[, -c(1,13)], INDICES = plotdf$histology, FUN = colSums)))
plotdf_tt <- plotdf_tt / rowSums(plotdf_tt)
plotdf_tt$histology <- rownames(plotdf_tt)
plotdf_tt_melt <- melt(data = plotdf_tt, id.vars = 12, measure.vars = 1:11)
flevels <- plotdf_tt$histology[order(plotdf_tt$no_kat, decreasing = F)]
plotdf_tt_melt$histology <- factor(plotdf_tt_melt$histology, levels = flevels, labels = paste0(flevels, " (", samplecounts[flevels], ")"))

# samplecountsdf2 <- do.call(rbind, by(data = plotdf, INDICES = plotdf$histology, FUN = function(x) data.frame(histology = unique(x$histology), ntot = nrow(x))))
# samplecountsdf2$histology <- factor(samplecountsdf2$histology, levels = levels(plotdf_tt_melt$histology), labels = paste0(levels(plotdf_tt_melt$histology), " (", samplecountsdf2[levels(plotdf_tt_melt$histology), "ntot"], ")"))

pl2 <- ggplot(data = plotdf_tt_melt) + geom_bar(stat = "identity", mapping = aes(x = histology, y = value, fill = variable), position = position_stack(reverse = T))
pl2 <- pl2 + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 6), panel.grid.major.x = element_blank(), legend.position = "none") + labs(x = "", y = "fraction")
pl2 <- pl2 + scale_fill_manual(values = c("#fbb4ae", "#FDD2CE", "#b3cde3", "#D1E1EE","#ccebc5", "#E0F3DC", "#decbe4", "#EBE0EF", "#E8E8E8", "#F1F1F1", "#FFFFFF"))

# pl2 <- pl2 + geom_text(data = samplecountsdf2, mapping = aes(x = samplecountsdf2$histology, y = 1.05, label = ntot), size = 2.11667)

pl2

ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_kataegis_tumortypes_signatures.pdf",
       plot = pl2, width = 183, height = 40, units = "mm")


sum(katresults[katresults$histology %in% c("Lymph-BNHL"), "p_aid_adj"] <= .1)
sum(katresults[katresults$histology %in% c("Lymph-BNHL"), "p_aid_adj"] > .1)



### kataegis "intensity" plot - Ludmill style
sinaplotdf2 <- katresults_repsamples %>% group_by(sample, histology) %>% summarise(no_foci = length(total))
### cervix adeno is basically only contributor to other
# sinaplotdf2[sinaplotdf2$histology == "Cervix-AdenoCA", "histology"] <- "Other"
sinaplotdf2 <- sinaplotdf2[!sinaplotdf2$histology %in% infreq_types, ]
sinaplotdf2$histology <- factor(x = sinaplotdf2$histology, levels = flevels)
sinaplotdf2 <- sinaplotdf2[order(sinaplotdf2$histology, sinaplotdf2$no_foci, decreasing = F), ]
sinaplotdf2$rnk <- unlist(by(data = sinaplotdf2, INDICES = sinaplotdf2$histology, FUN = function(x) seq(0, .7, .7/(nrow(x)-1))))
sinaplotdf2[is.nan(sinaplotdf2$rnk), "rnk"] <- .35
sinaplotdf2$rnk <- sinaplotdf2$rnk + as.numeric(sinaplotdf2$histology)

write.table(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_kataegis_tumortypes_intensities.txt", x = sinaplotdf2, quote = F, sep = "\t", row.names = F)

meddf <- do.call(rbind, by(data = sinaplotdf2, INDICES = sinaplotdf2$histology, FUN = function(x) data.frame(xmin = min(x$rnk), xmax = max(x$rnk), y = median(x = x$no_foci))))
meddf$histology <- rownames(meddf)
meddf["CNS-PiloAstro", c("xmin", "xmax")] <- meddf["CNS-PiloAstro", c("xmin", "xmax")] + c(-.35,.35)


cvect <- pcawg.colour.palette(x = tolower(sub(pattern = "-", replacement = ".", x = levels(sinaplotdf2$histology))), scheme = "tumour.subtype")
names(cvect) <- levels(sinaplotdf2$histology)
cvect[c("Skin-Melanoma-Acral", "Skin-Melanoma-Cut", "Kidney-RCC-Clear", "Other", "Kidney-RCC-Pap")] <- c("#000000", "#000000", '#FF4500', '#DDCDCD', '#FF4500')

psum4 <- ggplot(data = sinaplotdf2, mapping = aes(x = rnk, y = no_foci)) + geom_point(mapping = aes(fill = histology), stroke = .25, shape = 21, colour = "black", alpha = .75, size = .75) + scale_y_log10()
# psum4 <- ggplot(data = sinaplotdf2, mapping = aes(x = rnk, y = no_foci)) + geom_point(colour = "black", alpha = .75, size = .25) + scale_y_log10()
psum4 <- psum4 + geom_segment(data = meddf, mapping = aes(x = xmin, xend = xmax, y = y, yend = y), alpha = .5)
psum4 <- psum4 + scale_fill_manual(values = cvect, guide = F) + labs(y = "# foci")
psum4 <- psum4 + theme_minimal() + theme(axis.text.x = element_blank(), text = element_text(size = 6), axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
psum4 <- psum4 + annotation_logticks(sides = "l", scaled = T, colour = "grey")
psum4

ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_kataegis_tumortypes_intensities.pdf",
       plot = psum4, width = 183, height = 20, units = "mm")





#### writing final output for jamboree
katresults_jamboree <- merge(x = allpcfout, y = katresults, all = T)
katresults_jamboree$is_punctuated <- !is.na(katresults_jamboree$sigsum_SV)
katresults_jamboree$sv_assoc <- F
katresults_jamboree[which(katresults_jamboree$sv_dist <= 1e3), "sv_assoc"] <- T
# katresults_jamboree$has_AID_sig <- "p_aid_adj"

write.table(x = katresults_jamboree[, c("sample", "histology", "is_preferred", "chr", "start", "end", "total", "p_clonal", "p_subclonal", "p_early", "p_late", "p_na", "p_streak_adj", "no_phased_muts", "no_subclonal_muts", "no_antiphased_muts", "p_aid_adj", "sv_assoc", "signature")],
            file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_Kataegis_calls_JD.txt", quote = F, sep = "\t", row.names = F, col.names = T)

write.table(x = katresults_jamboree, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_Kataegis_calls_JD_allcolumns.txt", quote = F, sep = "\t", row.names = F, col.names = T)


### merging annotated mutation calls for jamboree
allmutannotfiles <- list.files(path = RESULTSBASE, pattern = "_all_muts_annot.txt", full.names = T, recursive = T)
allmutannot <- do.call(rbind, lapply(X = allmutannotfiles, FUN = read.delim, header = T, sep = "\t", as.is = T, colClasses = c("integer", rep("character", 6), rep("integer", 2))))
write.table(x = allmutannot[, c("chr", "pos", "ref", "alt", "trinuc", "sample", "histology", "start", "end")],
            file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_Kataegis_SNVs_JD.txt", quote = F, sep = "\t", row.names = F, col.names = T)




#### intersecting with driver calls
katresults_jamboree <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_Kataegis_calls_JD_allcolumns.txt", as.is = T)

allmutsfiles <- list.files(path = RESULTSBASE, pattern = "_all_muts.txt", full.names = T, recursive = T)
allmutsout <- lapply(X = allmutsfiles, FUN = read.delim, header = T, sep = "\t", as.is = T)
allmutsout <- data.frame(sample = rep(x = sub(pattern = "_all_muts.txt", replacement = "", x = basename(allmutsfiles)), times = sapply(X = allmutsout, FUN = nrow)), 
                    do.call(rbind, lapply(allmutsout, FUN = subset, select = c(chromosome, pos, ref, alt, trinuc, foci, ccf, major_cn, minor_cn, mcn, mult, timing, prob_clonal_early, prob_clonal_late, prob_subclonal))))

## quick check for overlap:
drivers <- read.delim(file = "/srv/shared/vanloo/ICGC_driver/TableS2_driver_point_mutations_annotation_20180110.txt", as.is = T)
drivers$gene <- sub(pattern = "::.*$", replacement = "", x = sub(pattern = "^.*?::.*?::", replacement = "", x = drivers$gene_id))
hypermut_drivers <- merge(x = drivers, y = allmutsout, by.x = c("sample", "chr", "pos", "ref", "alt"), by.y = c("sample", "chromosome", "pos", "ref", "alt"))

annot_driver <- function(df, focus) {
  katannot <- df[df$sample == focus$sample & df$chr == focus$chr &
                                      df$start <= focus$pos & df$end >= focus$pos, ]
  if (nrow(katannot) == 0)
    return(NULL)
  return(cbind(focus, katannot))
}

kataegis_drivers <- do.call(rbind, lapply(split(x = hypermut_drivers, f = 1:nrow(hypermut_drivers)), FUN = annot_driver, df = katresults_jamboree))

write.table(x = kataegis_drivers, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_drivers_and_recurrence_kataegis_drivers.txt", quote = F, sep = "\t", col.names = T, row.names = F)

kataegis_drivers_out <- kataegis_drivers[, c("sample", "histology", "gene", "p_aid_adj", "signature", "is_punctuated", "sv_assoc")]
kataegis_drivers_out_collapse <- do.call(rbind, by(data = kataegis_drivers_out, INDICES = paste0(kataegis_drivers_out$sample, "_", kataegis_drivers_out$gene), 
                           FUN = function(x) {y <- x[1, , drop = F]; y$p_aid_adj <- any(x$p_aid_adj <= .1); y$is_punctuated <- names(sort(table(x$signature), decreasing = T))[1]; y$is_punctuated <- all(x$is_punctuated); y$sv_assoc <- any(x$sv_assoc); return(y) }))
write.table(x = kataegis_drivers_out_collapse, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_drivers_and_recurrence_kataegis_drivers_collapsed.txt", quote = F, sep = "\t", col.names = T, row.names = F)

# ####
# 
# kataegis_drivers_uniq <- kataegis_drivers[!duplicated(kataegis_drivers$sample) & !duplicated(kataegis_drivers$start), ]
# # kataegis_drivers_uniq[, c("signature", "timing")] <- katres_jmb_repsamples_redotiming[match(paste0(kataegis_drivers_uniq$sample, "_", kataegis_drivers_uniq$focus_start), paste0(katres_jmb_repsamples_redotiming$sample, "_", katres_jmb_repsamples_redotiming$start)), c("sigsummary", "timing_fin")]
# # write.table(x = kataegis_drivers_uniq, file = "drivers_and_recurrence/20171120_kataegis_drivers_unique.txt", quote = F, sep = "\t", col.names = T, row.names = F)
# # plotdriverdf <- as.data.frame(table(kataegis_drivers_uniq$histology, kataegis_drivers_uniq$gene, kataegis_drivers_uniq$signature))
# # # plotdriverdf$Freq <- ifelse(plotdriverdf$Freq == 0, NA, plotdriverdf$Freq)
# # plotdriverdf <- plotdriverdf[plotdriverdf$Freq != 0, ]
# # 
# # drivp1 <- ggplot(data = plotdriverdf, mapping = aes(x = Var2, y = Var1)) + geom_point(aes(size = Freq, colour = Var3), alpha = .5)
# # drivp1 <- drivp1 + theme_minimal() + theme(axis.text.x = element_text(angle = 90))
# # drivp1
# 
# kataegis_drivers_uniq$focaltiming <- apply(kataegis_drivers_uniq[, c("p_sub", "p_early", "p_late", "p_na")], MARGIN = 1, FUN = function(x) c("p_sub", "p_early", "p_late", "p_na")[which.max(x)])
# # khypermut_drivers$focaltiming <- apply(hypermut_drivers[, c("p_sub", "p_early", "p_late", "p_na")], MARGIN = 1, FUN = function(x) c("p_sub", "p_early", "p_late", "p_na")[which.max(x)])
# drivp2 <- ggplot(data = kataegis_drivers_uniq[, -c(1,2)], mapping = aes(x = gene, y = histology)) + geom_jitter(mapping = aes(fill = signature, shape = focaltiming), colour = "grey", size = 3, alpha = .6, stroke = 1, width = .25, height = .25)
# # drivp2 <- drivp2 + geom_jitter(data = kataegis_drivers_uniq[kataegis_drivers_uniq$gene == "MYC", ], mapping = aes(shape = focaltiming, fill = signature), colour = "grey", size = 3, alpha = .6, stroke = 1)
# drivp2 <- drivp2 + theme_minimal() + theme(axis.text.x = element_text(angle = 90))
# drivp2 <- drivp2 + scale_shape_manual(values = c(24,25,23,21)) + scale_fill_manual(values = c("#b3cde3", "#fbb4ae", "#decbe4", "grey"))
# drivp2
# ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_drivers_and_recurrence_kataegis_drivers_plot.pdf", plot = drivp2, width = 7, height = 3.25)
# 
# 
# 
# 
# 
# 
# # ## checking correlation between sig1-sig2
# # sigsplit <- strsplit(x = katresults$active_sigs2, split = ",")
# # sigsplit <- as.data.frame(do.call(rbind, sigsplit[sapply(X = sigsplit, FUN = length) == 2]))
# # heatmap(table(sigsplit$V1, sigsplit$V2))
# # sigsplit_sub <- sigsplit[sigsplit$V1 %in% acceptsigs & sigsplit$V2 %in% acceptsigs,]
# # 
# # pl3 <- ggplot(data = rbind(sigsplit_sub, sigsplit_sub), mapping = aes(x = V1, y = V2)) + geom_jitter(alpha = .4, shape = ".", position = position_jitter(width = .1, height = .1))
# # pl3 <- pl3 + theme(axis.text.x = element_text(angle = 90))
# # pl3
# 
# ## generate Ludmill plots of trinuc freq + change for every type, merge if too close
# ## based on correlations would group 17a/b/28, 2/13 and 9/34
# 
# 
# 
# ## get Ludmill-type plots of "kataegis-signatures"
# katres_gr <- GRanges(seqnames = katres_jmb_repsamples$chr, ranges = IRanges(start = katres_jmb_repsamples$start, end = katres_jmb_repsamples$end), 
#                      mcols = katres_jmb_repsamples[, -c(3:5)], seqinfo = genome@seqinfo)
# # mcols(katres_gr)$mcols.signature <- ifelse(mcols(katres_gr)$mcols.var_expl2 >= 25, mcols(katres_gr)$mcols.signature, "uncertain")
# 
# sample_ids <- unique(katres_jmb_repsamples$sample)
# sample_muts <- lapply(X = sample_ids, FUN = function(x) read.table(file = file.path("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/analysed/", x, paste0(x, "_all_muts.txt")), 
#                                                                    header = T, sep = "\t", as.is = T))
# names(sample_muts) <- sample_ids
# sample_muts <- lapply(X = sample_muts, FUN = function(df) GRanges(seqnames = df$chromosome, IRanges(start = df$pos, end = df$pos), mcols = df[, -c(1,2)], seqinfo = genome@seqinfo))
# 
# annotate_muts_with_signature <- function(mutlist, focilist, sample_id) {
#   mutlist_sample <- mutlist[[sample_id]]
#   focilist_sample <- focilist[mcols(focilist)$mcols.sample == sample_id, ]
#   mut_foci_overlaps <- findOverlaps(query = mutlist_sample, subject = focilist_sample)
#   mutlist_sample_annotated <- mutlist_sample[queryHits(mut_foci_overlaps)]
#   mcols(mutlist_sample_annotated) <- DataFrame(mcols(mutlist_sample_annotated)[, 1:3],
#                                                mcols(focilist_sample[subjectHits(mut_foci_overlaps)])[, c("mcols.sample", "mcols.histology", "mcols.signature")],
#                                                start(focilist_sample[subjectHits(mut_foci_overlaps)]),
#                                                end(focilist_sample[subjectHits(mut_foci_overlaps)]) )
#   colnames(mcols(mutlist_sample_annotated)) <- c("ref", "alt", "trinuc", "sample", "histology", "signature", "focus_start", "focus_end")
#   return(mutlist_sample_annotated)
# }
# 
# # debug(annotate_muts_with_signature)
# # annotate_muts_with_signature(mutlist = sample_muts, focilist = katres_gr, sample_id = sample_ids[1])
# 
# sample_muts_antd <- lapply(X = sample_ids, FUN = annotate_muts_with_signature, mutlist = sample_muts, focilist = katres_gr)
# sample_muts_antdf <- do.call("c", sample_muts_antd)
# sample_muts_antdf <-  data.frame(chr = seqnames(sample_muts_antdf), pos = start(sample_muts_antdf), mcols(sample_muts_antdf))
# # sample_muts_antdf <- as.data.frame(unlist(sample_muts_antd, recursive = TRUE, use.names = TRUE)
# # sample_muts_antdf <- sample_muts_antdf[!grepl(pattern = "Lymph", x = sample_muts_antdf$histology, ignore.case = T), -c(4:5)]
# # colnames(sample_muts_antdf) <- "chr"
# 
# 
# # focalregions_sig2 <- katres_jmb_repsamples[katres_jmb_repsamples$signature == "Signature.2", ]
# trinuc_norms <- lapply(X = paste0("Signature.", c(2, 13, 19, "17a", "17b", 28, 9, 34)), FUN = function(x) get_trinuc_normalisation_factors(regions = katres_jmb_repsamples[katres_jmb_repsamples$signature == x, ], bsgenome = genome, overall = T))
# names(trinuc_norms) <- paste0("Signature.", c(2, 13, 19, "17a", "17b", 28, 9, 34))
# # trinuc_norms[1]
# 
# 
# lapply(X = paste0("Signature.", c(2, 13, 19, "17a", "17b", 28, 9, 34)), 
#        FUN = function(x) plot_mutationspectrum(mutations = sample_muts_antdf[sample_muts_antdf$signature == x, ], trinuc_freq = trinuc_norms[[x]], sample = "allsamples", histol = x, outdir = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis", suffix = paste0("kataegis_",x)))
# 
# trinuc_norm_2_13 <- get_trinuc_normalisation_factors(regions = katres_jmb_repsamples[katres_jmb_repsamples$signature %in% c("Signature.2", "Signature.13"), ], bsgenome = genome, overall = T)
# plot_mutationspectrum(mutations = sample_muts_antdf[sample_muts_antdf$signature %in% c("Signature.2", "Signature.13"), ], trinuc_freq = trinuc_norm_2_13, sample = "allsamples", histol = "APOBEC", outdir = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis", suffix = paste0("kataegis_" , "APOBEC"))
# 
# trinuc_norm_17_28 <- get_trinuc_normalisation_factors(regions = katres_jmb_repsamples[katres_jmb_repsamples$signature %in% c("Signature.17a", "Signature.17b") & !grepl(x = katres_jmb_repsamples$histology, pattern = "Lymph", ignore.case = T), ], bsgenome = genome, overall = T)
# plot_mutationspectrum(mutations = sample_muts_antdf[sample_muts_antdf$signature %in% c("Signature.17a", "Signature.17b")  & !grepl(x = sample_muts_antdf$histology, pattern = "Lymph", ignore.case = T), ], trinuc_freq = trinuc_norm_17_28, sample = "allsamples", histol = "TT", outdir = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis", suffix = paste0("kataegis_" , "TT_nolymph"))
# 
# trinuc_norm_17_28 <- get_trinuc_normalisation_factors(regions = katres_jmb_repsamples[katres_jmb_repsamples$signature %in% c("Signature.17a", "Signature.17b", "Signature.28", "Signature.34"), ], bsgenome = genome, overall = T)
# plot_mutationspectrum(mutations = sample_muts_antdf[sample_muts_antdf$signature %in% c("Signature.17a", "Signature.17b", "Signature.28", "Signature.34"), ], trinuc_freq = trinuc_norm_17_28, sample = "allsamples", histol = "TT", outdir = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis", suffix = paste0("kataegis_" , "TT_clean"))
# 
# 
# # plot_mutationspectrum(mutations = sample_muts_antdf[sample_muts_antdf$signature == "Signature.2", ], trinuc_freq = trinuc_norms[["Signature.2"]], sample = "allsamples", histol = "", outdir = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis", suffix = "kataegis_sig2_type")
# # undebug(plot_mutationspectrum)
# 
# 
# # bases <- c("A", "C", "G", "T")
# # # bases_fact <- factor(bases, levels = bases)
# # types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
# # types_fact <- factor(types, levels = types)
# # trinucleotides <- paste0(rep(rep(bases, rep(4,4)), 6),
# #                          rep(c("C", "T"), c(48, 48)),
# #                          rep(bases, 24))
# # trinucleotides_empty <- paste0(rep(rep(bases, rep(4,4)), 6),
# #                                rep(c(" "), c(96)),
# #                                rep(bases, 24))
# # trinucleotides_mutations <- paste0(paste0(rep(types, rep(16,6))), "_", trinucleotides)
# # trinucleotides_mutations_fact <- factor(trinucleotides_mutations, levels = trinucleotides_mutations)
# # 
# # if (!"trinuc" %in% colnames(mutations))
# 
# 
# # sample_muts_antdf$context <- 
# colnames(sample_muts_antdf)[1:2] <- c("chr", "pos")
# sample_muts_antdf$context <- as.character(get_trinuc_context(mutations = sample_muts_antdf, size = 10))
# colnames(sample_muts_antdf)[1:2] <- c("chr", "start")
# sample_muts_antdf[sample_muts_antdf$ref == "TRUE", "ref"] <- "T"
# sample_muts_antdf[sample_muts_antdf$alt == "TRUE", "alt"] <- "T"
# 
# # reverse complement and data augmentation
# revcomp <- data.frame(ref = create_complement(sample_muts_antdf$ref),
#                       alt = create_complement(sample_muts_antdf$alt),
#                       context = sapply(X = sample_muts_antdf$context, FUN = create_reverse_complement),
#                       trinuc = sapply(X = sample_muts_antdf$trinuc, FUN = create_reverse_complement),
#                       stringsAsFactors = F)
# sample_muts_antdf$trinuc_wc <- ifelse(sample_muts_antdf$ref %in% c("C", "T"), sample_muts_antdf$trinuc, revcomp$trinuc)
# write.table(x = sample_muts_antdf[sample_muts_antdf$signature %in% c("Signature.17a", "Signature.17b", "Signature.28", "Signature.34") , "context_wc"], file = "20170809_kataegis_context_TT.txt", quote = F, row.names = F, col.names = F)
# write.table(x = sample_muts_antdf[sample_muts_antdf$signature %in% c("Signature.19") , "context_wc"], file = "20170809_kataegis_context_altAPOBEC.txt", quote = F, row.names = F, col.names = F)
# write.table(x = sample_muts_antdf, file = "20170809_kataegis_allmuts.txt", quote = F, row.names = F, col.names = T, sep = "\t")
# 
# 
# get_locus_signatures_pnnls(mutations = sample_muts_antdf, signatures_renorm = pcawg_sigs_all, n = 1)
# debug(get_locus_signatures_pnnls)
# 
# 
# gordenin <- read.delim(file = "/srv/data/vanloo/jdemeul/Gordenin_APOBEC/MAF_Aug31_2016_sorted_A3A_A3B_comparePlus.txt")
# 
# p1 <- ggplot(data = katres_jamboree, mapping = aes(x = histology, y = switches/total, colour = histology)) + geom_point(alpha = .25, position = position_jitter(width = .1, height = 0)) + geom_boxplot(alpha = .5) 
# p1 <- p1 + theme(axis.text.x = element_text(angle = 90))
# p1
# 
# p1 <- ggplot(data = katres_jamboree, mapping = aes(x = var_expl, y = switches/total, colour = histology)) + geom_jitter(alpha = .25, position = position_jitter(width = .1, height = 0.1))
# p1 <- p1 + theme(axis.text.x = element_text(angle = 90))
# p1
# 
# # library(ggplot2)
# 
# 
# 
# 
# ggplot_df1 <- histology_all
# ggplot_df1$apobec_counts <- apobec_kataegis_counts[ggplot_df1$tumor_wgs_aliquot_id]
# ggplot_df1[is.na(ggplot_df1$apobec_counts), "apobec_counts"] <- 0
# ggplot_df1$class <- ifelse(ggplot_df1$apobec_counts == 0, 0,
#                            ifelse(ggplot_df1$apobec_counts < 10, "< 10", ">= 10"))
# 
# colnames(ggplot_df1)[1:3] <- c("sample", "histology", "whitelist")
# infreq_types <- names(which(table(ggplot_df1$histology) < 10))
# ggplot_df1[ggplot_df1$histology %in% infreq_types, "histology"] <- "other"
# 
# tt_order <- by(data = ggplot_df1$class, INDICES = ggplot_df1$histology, FUN = function(x) sum(x != "0")/nrow(x), simplify = T)
# ggplot_df1$histology <- factor(ggplot_df1$histology, levels = names(tt_order)[order(as.vector(tt_order), decreasing = T)])
# 
# 
# p5 <- ggplot(data = ggplot_df1, mapping = aes(x = histology)) + geom_bar(mapping = aes(fill = class), position = "fill")
# p5 <- p5 + theme_minimal() + theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank()) + ylab("Fraction")
# p5 <- p5 + geom_text(data = as.data.frame(table(ggplot_df1$histology)), mapping = aes(x = Var1, y = 1.05, label = Freq), angle = 90, size = 3)
# p5
# 
# ggsave(filename = "/srv/data/vanloo/jdemeul/ICGC/kataegis/20170516_apobec_tumortypes.pdf",
#        plot = p5, width = 20, height = 6)
# 
# 
# allpcfout_clean_mod <- allpcfout_clean[!is.na(allpcfout_clean$histology), ]
# allpcfout_clean_mod[allpcfout_clean_mod$histology %in% infreq_types, "histology"] <- "other"
# allpcfout_clean_mod$histology <- factor(allpcfout_clean_mod$histology, levels = levels(ggplot_df1$histology))
# 
# p7 <- ggplot(data = allpcfout_clean_mod, mapping = aes(x = histology)) + geom_bar(mapping = aes(fill = timing_fin), position = "fill")
# p7 <- p7 + theme_minimal() + theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank()) + ylab("Fraction") + scale_x_discrete(drop = F)
# p7 <- p7 + scale_fill_manual(values = c('#4daf4a','#984ea3','#377eb8', '#e41a1c', "grey"))
# p7 <- p7 + geom_text(data = as.data.frame(table(allpcfout_clean_mod$histology)), mapping = aes(x = Var1, y = 1.05, label = Freq), angle = 90, size = 3)
# p7
# 
# 
# ggsave(filename = "/srv/data/vanloo/jdemeul/ICGC/kataegis/20170516_apobec_tumortypes_timing.pdf",
#        plot = p7, width = 20, height = 6)
# 
# # 
# # 
# # timingcounts <- as.data.frame(unlist(lapply(by(data = allpcfout_clean_mod$timing_fin, INDICES = allpcfout_clean_mod$histology, table), function(x) x/sum(x))))
# # timingcounts <- cbind(timingcounts, do.call(rbind, strsplit(rownames(timingcounts), split = ".", fixed = T)))
# # colnames(timingcounts) <- c("frequency", "histology", "timing")
# # 
# # p8 <- ggplot(data = timingcounts, aes(x = histology, y = timing, fill = frequency)) + geom_raster() + scale_fill_continuous(low = "#2166ac", high = "#b2182b")
# # p8
# 
# 
# ## finding interesting cases (e.g. with both clonal early, late and subclonal events - but not too many)
# eventsdf <- as.data.frame(do.call(rbind, by(data = allpcfout_clean$timing_fin, INDICES = allpcfout_clean$sample, FUN = table)))
# # samplesnames <- names(eventsdf)
# represent_samples <- rownames(eventsdf)[(eventsdf$`clonal [early]` >= 1 & eventsdf$`clonal [late]` >= 1 & eventsdf$`subclonal` >= 1 & rowSums(eventsdf) <= 20)]
# allpcfout_clean_int <- allpcfout_clean[allpcfout_clean$sample %in% represent_samples, ]
# 
# ## selecting p-value cutoffs
# pcutoff <- 10^-(1:200)
# ratio <- vector(mode = "numeric", length = 200)
# ratio2 <- vector(mode = "numeric", length = 200)
# for (i in 1:200) {
#   no_subclonal <- sum(!is.na(allpcfout_clean$subcl_chisq_p) & allpcfout_clean$p_clonal_adj <= pcutoff[i], na.rm = T)
#   no_clonal <- sum(!is.na(allpcfout_clean$subcl_chisq_p) & allpcfout_clean$p_subclonal_adj <= pcutoff[i], na.rm = T)
#   ratio[i] <- (no_subclonal/ (no_subclonal+no_clonal))
#   
#   no_subclonal2 <- sum(!is.na(allpcfout_clean$subcl_chisq_p) & allpcfout_clean$p_clonal <= pcutoff[i], na.rm = T)
#   no_clonal2 <- sum(!is.na(allpcfout_clean$subcl_chisq_p) & allpcfout_clean$p_subclonal <= pcutoff[i], na.rm = T)
#   ratio2[i] <- (no_subclonal2/ (no_subclonal2+no_clonal2))
# }
# plot(1:200, ratio, ylim = c(0,1), col = "red")
# points(1:200, ratio2, col = "green")
# 
# allpcfout_clean$p_clonal_adj <- p.adjust(p = allpcfout_clean$p_clonal, method = "fdr")
# allpcfout_clean$p_subclonal_adj <- p.adjust(p = allpcfout_clean$p_subclonal, method = "fdr")
# 
# sum(!is.na(allpcfout_clean$subcl_chisq_p) & (allpcfout_clean$p_clonal_adj <= .0001 | allpcfout_clean$p_subclonal_adj <= .0001), na.rm = T)
# sum(!is.na(allpcfout_clean$subcl_chisq_p), na.rm = T)
# 
# 
# 
# ## annotation and checking overlaps with TSG exons
# library(GenomicRanges)
# library(biomaRt)
# 
# census <- read.delim("/srv/data/vanloo/jdemeul/refdata/Census_allMon Feb 13 09-36-09 2017.tsv", as.is = T)
# census <- census[grepl("Rec", census$Molecular.Genetics) | grepl("TSG", census$Role.in.Cancer),]
# 
# ensembl37 <- useMart(host = "grch37.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL")
# ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl37)
# 
# clean_foci_range <- GRanges(seqnames = allpcfout_clean$chr, ranges = IRanges(start = allpcfout_clean$start,
#                                                                              end = allpcfout_clean$end),
#                             mcols = DataFrame(allpcfout_clean[, -c(3:5)]))
# 
# ## kataegis ranges
# out <- getBM(attributes = c("chromosome_name", "exon_chrom_start", "exon_chrom_end", "external_gene_name", "ensembl_exon_id")
#              , filters = c("external_gene_name"), values = list(census$Gene.Symbol), mart = ensembl)
# genes_ranges <- GRanges(seqnames = out$chromosome_name, IRanges(start = out$exon_chrom_start,
#                                                                 end = out$exon_chrom_end),
#                         mcols = DataFrame(geneID = out$external_gene_name))
# # names(genes_ranges) <- out$external_gene_name
# seqlevels(genes_ranges, force=TRUE) <- c(1:22,"X")
# 
# hits <- findOverlaps(query = clean_foci_range, subject = genes_ranges)
# hitting_foci <- clean_foci_range[queryHits(hits)]
# mcols(hitting_foci) <- DataFrame(mcols(hitting_foci), tsg = mcols(genes_ranges)$mcols.geneID[subjectHits(hits)])
# hitting_foci <- hitting_foci[!duplicated(hitting_foci)]
# hitting_foci_df <- as.data.frame(hitting_foci)
# write.table(x = hitting_foci_df, file = "/Users/demeulj/Documents/Work/2017_PCAWG-kataegis/kataegis_fun_samples/hitTSGs.txt", quote = F, sep = "\t", row.names = F)
# 
# 
# ## phasing info analysis
# RESULTSFILE <- "/srv/data/vanloo/jdemeul/ICGC/kataegis/20170516_Kataegis_Results_all_finalpvals.txt"
# 
# allpcfout <- read.delim(file = RESULTSFILE, as.is = T)
# allpcfout_subs <- allpcfout[allpcfout$p_clonal_fin_adj <= .05 & allpcfout$p_subclonal_fin_adj <= .05 & !is.na(allpcfout$p_subclonal_fin_adj), ]
# 
# # ggd.qqplot(pvector = katresults$p_subclonal_fin_adj[!is.na(katresults$p_subclonal_fin_adj) & katresults$p_subclonal_fin_adj > 0])
# # 
# # plot(-log10(katresults$p_subclonal_fin), -log10(katresults$p_clonal_fin))
# # cor(katresults$p_subclonal_fin, katresults$p_clonal_fin, use = "pairwise.complete.obs")
# # pl1 <- ggplot(data = katresults, mapping = aes(x = -log10(p_subclonal_fin_adj + 1e-7), y = -log10(p_clonal_fin_adj+ 1e-7))) + geom_density2d(mapping = aes(x = -log10(p_subclonal_fin_adj), y = -log10(p_clonal_fin_adj))) + geom_point(alpha =.25, mapping = aes(colour = abs(log10((p_clonal_fin+1e-7)/(p_subclonal_fin+1e-7))) > 2))
# # pl1 <- pl1 + geom_vline(xintercept = -log10(.05)) + geom_hline(yintercept = -log10(.05))
# # pl1 <- pl1 + geom_vline(xintercept = -log10(.01), color = "red") + geom_hline(yintercept = -log10(.1), color = "red")
# # pl1 <- pl1 + geom_vline(xintercept = -log10(.1), color = "green") + geom_hline(yintercept = -log10(.01), color = "green")
# # # pl1 <- pl1 + geom_abline(slope = 1, intercept = -log10(.05))
# # # pl1 <- pl1 + geom_density2d()
# # pl1
# # 
# # qplot(x = -log10(p_subclonal_fin), y = -log10(p_clonal_fin), data = katresults, geom = "density2d")
# # 
# # ggd.qqplot = function(pvector, main=NULL, ...) {
# #   o = -log10(sort(pvector,decreasing=F))
# #   e = -log10( 1:length(o)/length(o) )
# #   plot(e,o,pch=19,cex=1, main=main, ...,
# #        xlab=expression(Expected~~-log[10](italic(p))),
# #        ylab=expression(Observed~~-log[10](italic(p))),
# #        xlim=c(0,max(e)), ylim=c(0,max(o)))
# #   lines(e,e,col="red")
# # }