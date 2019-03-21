### trying a linead mixed-effects model on to predict kataegis intensity

# Load dataset, inspect size and additional info

library(dplyr)
library(ggplot2)
# library(lme4)
library(blme)
# library(optimx)
# library(lmerTest)
# library(glmmTMB)
# library(nlme)
# library(parallel)


BASEOUT <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/final_rerun_annotmuts/"
TOPHITSDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/tophits_rainfall/"
CLEANHISTOLOGYFILE <- "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/gene_conversion/results/summary_table_combined_annotations_v2_JD.txt"
# ICGC_annotations/summary_table_combined_annotations_v2.txt"

source(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/kataegis_functions.R", local = T)
source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/code_kataegis/pcawg.colour.palette.R")

histology_all <- read_histology(histologyfile = CLEANHISTOLOGYFILE)

katresults_jamboree <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_Kataegis_calls_JD_allcolumns.txt", as.is = T)
katresults_punct <- katresults_jamboree[katresults_jamboree$is_punctuated, ]

### kataegis "intensity" plot - Ludmill style
sinaplotdf2 <- katresults_punct %>% group_by(sample, histology) %>% summarise(no_foci = length(total), no_sv_foci = sum(sv_assoc), no_sv_foci_apo = sum(signature == "APO" & sv_assoc), no_non_sv_foci_apo = sum(signature == "APO" & !sv_assoc), no_apo_foci = sum(signature == "APO"), no_apo_foci = sum(signature == "APO"), no_ctt_foci = sum(signature == "CTT"))

sinaplotdf2 <- merge(x = sinaplotdf2, y = histology_all[, c("samplename", "histology_abbreviation", "inferred_sex", "is_preferred", "donor_age_at_diagnosis", "histology_tier4")], by.x = c("sample", "histology"), by.y = c("samplename", "histology_abbreviation"), all = T)
rare_types <- names(which(sort(table(sinaplotdf2$histology)) < 10))
sinaplotdf2 <- sinaplotdf2[sinaplotdf2$is_preferred & !sinaplotdf2$histology %in% rare_types, ]
sinaplotdf2[is.na(sinaplotdf2$no_foci), "no_foci"] <- 0
sinaplotdf2[is.na(sinaplotdf2$no_sv_foci), "no_sv_foci"] <- 0
sinaplotdf2[is.na(sinaplotdf2$no_apo_foci), "no_apo_foci"] <- 0
sinaplotdf2[is.na(sinaplotdf2$no_non_sv_foci_apo), "no_non_sv_foci_apo"] <- 0
sinaplotdf2[is.na(sinaplotdf2$no_sv_foci_apo), "no_sv_foci_apo"] <- 0
### cervix adeno is basically only contributor to other
# sinaplotdf2 <- sinaplotdf2[!sinaplotdf2$histology %in% infreq_types, ]
sinaplotdf2$histology <- factor(x = sinaplotdf2$histology, levels = sort(unique(sinaplotdf2$histology)))
# sinaplotdf2$donor_age_at_diagnosis_grouped <- round(sinaplotdf2$donor_age_at_diagnosis/10)
sinaplotdf2$donor_age_at_diagnosis <- scale(log(sinaplotdf2$donor_age_at_diagnosis), center = T, scale = T)


####### adding in number of SVs/chromothripsis

ctcalls <- read.delim(file = "/srv/shared/vanloo/home/mtarabichi/PCAWG/chromothripsis/tableCT.Step6.022019.txt", as.is = T)
sinaplotdf2$hasCT <- sinaplotdf2$sample %in% ctcalls$samplename
sinaplotdf2$numSVs <- rowSums(histology_all[match(x = sinaplotdf2$sample, table = histology_all$samplename), c("num_clonal_svs", "num_subclonal_svs")])
sinaplotdf2[is.na(sinaplotdf2$numSVs), "numSVs"] <- 0
sinaplotdf2$numSVssc <- scale(log(sinaplotdf2$numSVs+1), center = T, scale = T)


drivers <- read.delim(file = "/srv/shared/vanloo/ICGC_driver/TableS3_panorama_drivers_pcawg_20180110.txt", as.is = T)
# drivers$gene <- sub(pattern = "::.*$", replacement = "", x = sub(pattern = "^.*?::.*?::", replacement = "", x = drivers$gene_id))
# commondriv <- names(which(table(drivers$gene) >= 1))
commondriv <- names(which(table(drivers$gene[!duplicated(paste0(drivers$sample_id, "_", drivers$gene))]) >= 1))


Xmat <- as.data.frame(do.call(cbind, lapply(commondriv, function(x, driv, katsam) katsam %in% driv[drivers$gene == x, "sample_id"], driv = drivers, katsam = sinaplotdf2$sample)))
commondriv <- paste0("x", gsub(pattern = "-", replacement = "", x = gsub(pattern = "[\\._\\-\\:]", replacement = "", x = commondriv)))

colnames(Xmat) <- commondriv



sinaplotdf2 <- cbind(sinaplotdf2, Xmat)
# sinaplotdf2$no_apo_foci <- log(sinaplotdf2$no_apo_foci+1)
# sinaplotdf2$no_foci <- log(sinaplotdf2$no_foci+1)
# sinaplotdf2$no_sv_foci <- log(sinaplotdf2$no_sv_foci+1)


### pull in APOBEC/AID expression
apoexpr <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/data/RNASeq_joint_fpkm_uq_APOBECs.tsv", as.is = T)
apoexpr$feature <- sapply(X = strsplit(x = apoexpr$feature, split = "::"), FUN = '[', 2)
# apoexpr <- t(apoexpr)


reltab <- read.delim(file = "/srv/shared/vanloo/ICGC_annotations/release_may2016.v1.4.tsv", as.is = T)
# nchar(reltab$tumor_rna_seq_aliquot_id) > 36
reltab <- reltab[reltab$tumor_rna_seq_aliquot_id != "", ]
# check for potential sample number mismatches
# reltab[reltab$tumor_wgs_specimen_count > 1, "tumor_rna_seq_aliquot_id"]
wgs_aliq_ids <- strsplit(reltab$tumor_wgs_aliquot_id, split = ",")

reltab_redup <- reltab[rep(1:nrow(reltab), lengths(wgs_aliq_ids)),]
# reltab_redup$tumor_wgs_icgc_specimen_id <- unlist(strsplit(reltab_redup$tumor_wgs_icgc_specimen_id, split = ","))
reltab_redup$tumor_wgs_aliquot_id <- unlist(wgs_aliq_ids)
reltab_redup$tumor_rna_seq_aliquot_id <- sapply(X = strsplit(reltab_redup$tumor_rna_seq_aliquot_id, split = ","), '[', 1)
reltab <- reltab_redup[, c("tumor_wgs_aliquot_id", "tumor_rna_seq_aliquot_id")]


remoddf <- do.call(rbind, lapply(X = gsub(pattern = "-", replacement = ".", x = reltab$tumor_rna_seq_aliquot_id), FUN = function(smpl, df) df[, grep(pattern = smpl, x = colnames(df))], df = apoexpr))
# sapply(1:13, function(x, df) min(df[,x][df[,x]>0]), df = remoddf)
remoddf <- as.data.frame(scale(log2(remoddf + 1e-3), center = T, scale = T))
colnames(remoddf) <- apoexpr$feature
remoddf$sample <- reltab$tumor_wgs_aliquot_id

# colnames(remoddf)[-c(6,14)]

sinaplotdf3 <- merge(x = sinaplotdf2, y = remoddf, by.x = "sample", by.y = "sample", all.x = T)





# lm1 <- lm(data = sinaplotdf2[, -c(1,3,4)], formula = no_apo_foci ~ histology*inferred_sex)
# summary(lm1)
# par(mfrow = c(2,2))
# plot(lm1)
# 
# 
# GLM <- gls(data = sinaplotdf2[, -c(1,3,4)], model = no_apo_foci ~ histology*TP53,
#            method = "ML")
# summary(GLM)
# 
# lme()
# 
# lmm1 <- lme(data = sinaplotdf2[, -c(1,3,4)], fixed = no_apo_foci ~ PTEN,
#             random = ~1|histology, method = "ML")

tophits <- c(colnames(Xmat)[order(colSums(Xmat), decreasing = T)], colnames(remoddf)[-c(5,14)])
# tophits
# paste0(tophits, collapse = "+")
# tophits <- tophits[-which(tophits %in% c("xTP53"))]
testvar <- function(driver, df) {
  newrow <- data.frame(Estimate = numeric(1L),
                       Stderr = numeric(1L),
                       Zvalue = numeric(1L),
                       Pvalue = numeric(1L),
                       Err = character(1L),
                       Warn = character(1L))
  print(driver)
  # lmmo1 <- glmer(formula = paste0("no_apo_foci ~ xTP53 + ", driver, " + (1 + xTP53 + |histology) + (1|donor_age_at_diagnosis) + (1|inferred_sex)"), data = df, family = "poisson")
  # lmmo2 <- glmer(formula = paste0("no_apo_foci ~ xTP53 + (1 + xTP53 + |histology) + (1|donor_age_at_diagnosis) + (1|inferred_sex)"), data = df, family = "poisson")

  # anovaout <- tryCatch({
  #   # blmmo_fix2 <- update(blmmo_fix, as.formula(paste0(". ~ . + ", driver)))
  #   blmmo_fix <- glmmTMB(data = df, formula = as.formula(paste0("no_apo_foci ~ numSVssc + donor_age_at_diagnosis + (1+",driver,"|histology)")),
  #                        family = "nbinom2", control = glmmTMBControl(profile=TRUE))
  #   blmmo_fix2 <- update(blmmo_fix, as.formula(paste0(". ~ . + ", driver)))
  #   bgout <- list(anva = anova(blmmo_fix, blmmo_fix2), sum = summary(blmmo_fix2), warn = character(1L))
  # }, warning = function(warn) {
  #   blmmo_fix <- glmmTMB(data = df, formula = as.formula(paste0("no_apo_foci ~ numSVssc + donor_age_at_diagnosis + (1+",driver,"|histology)")),
  #                        family = "nbinom2", control = glmmTMBControl(profile=TRUE))
  #   blmmo_fix2 <- update(blmmo_fix, as.formula(paste0(". ~ . + ", driver)))
  #   bgout <- list(anva = anova(blmmo_fix, blmmo_fix2), sum = summary(blmmo_fix2), warn = paste0("WARN: ", warn))
  #   return(bgout)
  # }, error = function(err) {
  #   # bgout <- list(anva = NULL)
  #   return(NULL)
  # })
  
  anovaout <- tryCatch({
    # blmmo_fix2 <- update(blmmo_fix, as.formula(paste0(". ~ . + ", driver)))
    blmmo_fix <- bglmer(data = df, formula = as.formula(paste0("no_apo_foci ~ numSVssc + donor_age_at_diagnosis + (1+",driver,"|histology)")),
                        family = "poisson", fixef.prior = normal, control = glmerControl(optimizer = "Nelder_Mead"))
    blmmo_fix2 <- update(blmmo_fix, as.formula(paste0(". ~ . + ", driver)))
    bgout <- list(anva = anova(blmmo_fix, blmmo_fix2), sum = summary(blmmo_fix2), warn = character(1L))
  }, warning = function(warn) {
    blmmo_fix <- bglmer(data = df, formula = as.formula(paste0("no_apo_foci ~ numSVssc + donor_age_at_diagnosis + (1+",driver,"|histology)")),
                        family = "poisson", fixef.prior = normal, control = glmerControl(optimizer = "Nelder_Mead"))
    blmmo_fix2 <- update(blmmo_fix, as.formula(paste0(". ~ . + ", driver)))
    bgout <- list(anva = anova(blmmo_fix, blmmo_fix2), sum = summary(blmmo_fix2), warn = paste0("WARN: ", warn))
    return(bgout)
  }, error = function(err) {
    # bgout <- list(anva = NULL)
    return(NULL)
  })
  return(anovaout)
  # return(anova(blmmo1,blmmo2))
}

# df <- sinaplotdf3
# blmmo_fix <- glmer.nb(formula = no_apo_foci ~ xTP53 + APOBEC3B + numSVssc + donor_age_at_diagnosis + (1+xTP53|histology),
                    # data = df, glmerControl(optimizer = "Nelder_Mead"))
# blmmo_fix


# # blmmo_fix2 <- bglmer(formula = no_apo_foci ~ xTP53 + numSVssc + (1|histology),
# #                           data = df, family = "poisson", fixef.prior = normal(cov = diag(1,3)), glmerControl(optimizer = "Nelder_Mead"))
# blmmo_fix <- bglmer(formula = no_apo_foci ~ numSVssc + APOBEC3B + (1|histology),
#                     data = df, family = "poisson", fixef.prior = normal(cov = diag(1,3)), glmerControl(optimizer = "Nelder_Mead"))
# temp <- summary(blmmo_fix)
# temp <- anova(blmmo_fix, blmmo_fix2)
# coef(blmmo_fix)
# allFit(blmmo_fix)

# 
# check zero-inflation
# p1 <- ggplot(data = sinaplotdf3, mapping = aes(x = no_apo_foci)) + geom_histogram(binwidth = 1) + facet_wrap(~histology, scales = "free_y") + xlim(c(0,30))
# p1
# 
# p1 <- ggplot(data = sinaplotdf3, mapping = aes(x = numSVssc, y = log(no_apo_foci+1))) + geom_point() + stat_smooth(method = "glm", method.args = list(family = "poisson"))# + facet_wrap(~histology, scales = "free")
# p1
# 
# 
# blmmo_fix <- glmmTMB(data = sinaplotdf3, formula = no_apo_foci ~ xTP53 + APOBEC3B + numSVssc + donor_age_at_diagnosis + (1|histology) + (1|inferred_sex),
#                       family = "nbinom2", control = glmmTMBControl(profile=TRUE))
# blmmo_fix <- glmmTMB(data = sinaplotdf3, formula = no_apo_foci ~ APOBEC3B + numSVssc + donor_age_at_diagnosis + (1|histology) + (1|inferred_sex),
#                      family = "nbinom2", control = glmmTMBControl(profile=TRUE))
# anova(blmmo_fix,blmmo_fix2)
# 
# summary(blmmo_fix)
# coef(blmmo_fix)
#   # glmer(formula = no_apo_foci ~ xTP53 + (1+xTP53|histology) + (1|inferred_sex) + (1|donor_age_at_diagnosis_grouped) + (1|numSVssc),
#   #                     data = sinaplotdf3, family = "poisson", glmerControl(optimizer = "Nelder_Mead"))
# 
# 
# blmmo_fix <- glmer.nb(formula = no_apo_foci ~ xTP53 + APOBEC3B + numSVssc + donor_age_at_diagnosis + (1+xTP53|histology) + (1|inferred_sex),
#                     data = sinaplotdf3, glmerControl(optimizer = "Nelder_Mead"))
# blmmo_fix
# 
# summary(blmmo_fix)
# blmmo_fix2 <- bglmer(formula = no_apo_foci ~ donor_age_at_diagnosis + numSVssc + (1+xTP53|histology) + (1|inferred_sex) + (1|hasCT),
#                     data = sinaplotdf2, family = "poisson", fixef.prior = normal(cov = diag(1,3)), glmerControl(optimizer = "Nelder_Mead"))
# 
# glmer(formula = no_apo_foci ~ xTP53 + donor_age_at_diagnosis + numSVssc + (1+xTP53|histology) + (1|inferred_sex),
#       data = sinaplotdf2, family = "poisson", glmerControl(optimizer = "Nelder_Mead"))

# summary(blmmo_fix)
# anova(blmmo_fix, blmmo_fix2)

# anovaresults <- lapply(X = tophits, FUN = testvar, df = sinaplotdf2, blmmo2 = blmmo2)
# undebug(testvar)
anovaresults <- lapply(X = tophits, FUN = testvar, df = sinaplotdf3)
# saveRDS(object = anovaresults, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/drivercorranovaresults.RDS")
# anovaresults <- readRDS(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/drivercorranovaresults.RDS")
# anovaresults <- mclapply(X = tophits[1:20], FUN = testvar, df = sinaplotdf3, blmmo_fix = blmmo_fix, mc.preschedule = T, mc.cores = 10)
fails <- sapply(anovaresults, FUN = is.null)
# anovaresults_df <- do.call(rbind, lapply(anovaresults, FUN = function(x) data.frame(AIC = x$anva$AIC[2], BIC = x$anva$BIC[2], Chisq = x$anva$Chisq[2], coef = x$sum$coefficients$cond[4,1], coefz = x$sum$coefficients$cond[4,3], pvalsum = x$sum$coefficients$cond[4,4], warn = x$warn)))
anovaresults_df <- do.call(rbind, lapply(anovaresults, FUN = function(x) data.frame(AIC = x$anva$AIC[2], BIC = x$anva$BIC[2], Chisq = x$anva$Chisq[2], pvalanva = x$anva$`Pr(>Chisq)`[2], coef = x$sum$coefficients[4,1], coefz = x$sum$coefficients[4,3], pvalsum = x$sum$coefficients[4,4], warn = x$warn)))

# anovaresults_df2$warn <- sapply(anovaresults, FUN = function(x) x$warn)
# anovaresults_df2$err <- sapply(anovaresults, FUN = function(x) x$err)

# anovaresults_df <- as.data.frame(do.call(rbind, anovaresults))

# anovaresults <- mclapply(X = tophits, FUN = testvar, df = sinaplotdf2, mc.preschedule = T, mc.cores = 15)
# basicmodel <- data.frame(AIC = anovaresults[[1]]$AIC[1], BIC = anovaresults[[1]]$BIC[1], logLik = anovaresults[[1]]$logLik[1], pval = anovaresults[[1]]$`Pr(>Chisq)`[2])
# anovaresults_df <- as.data.frame(do.call(rbind, anovaresults))
# resultsdf <- do.call(rbind, lapply(X = anovaresults, FUN = function(x) data.frame(AIC = x$AIC[2], BIC = x$BIC[2], logLik = x$logLik[2], pval = x$`Pr(>Chisq)`[2])))
# anovaresults_df2$driver <- tophits[1:20]


anovaresults_df$driver <- tophits[!fails]
# anovaresults_df <- anovaresults_df[anovaresults_df$warn == "", ]


# anovaresults_df[anovaresults_df$pvalsum > 0, "pvalsumadj"] <- p.adjust(anovaresults_df[anovaresults_df$pvalsum > 0, "pvalsum"], method = "fdr")
# anovaresults_df[anovaresults_df$pvalsum > 0, "pvalanvaadj"] <- p.adjust(anovaresults_df[anovaresults_df$pvalsum > 0, "pvalanva"], method = "fdr")
# resultsv_adj[resultsv_adj < .05]
# resultsv_adj[resultsv_adj < .1]
# blmmo1final <- bglmer(formula = no_apo_foci ~ APOBEC3B + numSVssc + donor_age_at_diagnosis + (1|histology),
#                data = sinaplotdf3, family = "poisson", fixef.prior = normal(cov = diag(1,4)), glmerControl(optimizer = "Nelder_Mead"))
# 

# blmmo_fix <- glmmTMB(data = sinaplotdf3, formula = no_apo_foci ~ numSVssc + donor_age_at_diagnosis + (1|histology),
#                      family = "nbinom2", control = glmmTMBControl(profile=TRUE))
blmmo_fix <- bglmer(data = sinaplotdf3, formula = no_apo_foci ~ numSVssc + donor_age_at_diagnosis + (1|histology),
                    family = "poisson", fixef.prior = normal, control = glmerControl(optimizer = "Nelder_Mead"))



# blmmo_fix
# coef(blmmo_fix)
# summary(blmmo_fix)
blmmo1finalsum <- summary(blmmo_fix)

# plot(blmmo1final)
# head(anovaresults_df)
# blmmo1final
# blmmo1finalsum$coefficients
# apo <- data.frame(AIC = numeric(2L), BIC = numeric(2L), Chisq = numeric(2L),
#                   coef = blmmo1finalsum$coefficients$cond[2:3,1], coefz = blmmo1finalsum$coefficients$cond[2:3,3],
#                   pvalsum = blmmo1finalsum$coefficients$cond[2:3,4], warn = character(2L),
#                   driver = rownames(blmmo1finalsum$coefficients$cond)[2:3])
apo <- data.frame(AIC = numeric(2L), BIC = numeric(2L), Chisq = numeric(2L), pvalanva = numeric(2L),
                  coef = blmmo1finalsum$coefficients[2:3,1], coefz = blmmo1finalsum$coefficients[2:3,3],
                  pvalsum = blmmo1finalsum$coefficients[2:3,4], warn = character(2L),
                  driver = rownames(blmmo1finalsum$coefficients)[2:3])

anovaresults_df <- rbind(anovaresults_df, apo)

anovaresults_df$pvalsumadj <- p.adjust(anovaresults_df$pvalsum, method = "fdr")
anovaresults_df$pvalanvaadj <- p.adjust(anovaresults_df$pvalanva, method = "fdr")


write.table(x = anovaresults_df, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190225_LMEM_NB2_no_apo_foci_vs_drivers_cytdeam_APO3B_nsvs_age_historand_BAYES.txt", quote = T, sep = "\t", col.names = T, row.names = F)
anovaresults_df <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190225_LMEM_NB2_no_apo_foci_vs_drivers_cytdeam_APO3B_nsvs_age_historand_BAYES.txt", stringsAsFactors = F)

p1 <- ggplot(data = anovaresults_df[anovaresults_df$warn == "", ], mapping = aes(x = coef, y = abs(coefz), colour = pvalanvaadj < .05)) + geom_point(data = anovaresults_df[anovaresults_df$driver %in% colnames(remoddf),], colour = "blue", size = 2, show.legend = F) + geom_point(show.legend = F)
p1 <- p1 + geom_text(mapping = aes(label = gsub(pattern = "^x", replacement = "", x = ifelse(pvalanvaadj < .05, driver, "")), color = driver %in% c("APOBEC3B", "donor_age_at_diagnosis", "numSVssc")), nudge_y = 1, show.legend = F) + theme_minimal() + scale_colour_manual(values = c('TRUE' = "red", 'FALSE' = "grey"))
p1 <- p1 + labs(x = "Coefficient estimate", y = "|Z-value|") + xlim(c(-2.5, 2.5))
p1

ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190225_LMEM_NB2_no_apo_foci_vs_drivers_cytdeam_APO3B_nsvs_age_historand_BAYES.pdf",
       plot = p1, width = 122, height = 91.5, units = "mm", useDingbats=FALSE)


# p1 <- ggplot(data = anovaresults_df[anovaresults_df$warn == "", ], mapping = aes(x = coef, y = abs(coefz), colour = pvalsumadj < .05)) + geom_point(data = anovaresults_df[anovaresults_df$driver %in% colnames(remoddf),], colour = "blue", size = 2, show.legend = F) + geom_point(show.legend = F)
# p1 <- p1 + geom_text(mapping = aes(label = gsub(pattern = "^x", replacement = "", x = ifelse(pvalsumadj < .05, driver, "")), color = driver %in% c("APOBEC3B", "donor_age_at_diagnosis", "numSVssc")), nudge_y = 0.15, show.legend = F) + theme_minimal() + scale_colour_manual(values = c('TRUE' = "red", 'FALSE' = "grey"))
# p1 <- p1 + labs(x = "Coefficient estimate", y = "|Z-value|")
# p1
# 
# ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190225_LMEM_NB2_no_apo_foci_vs_drivers_cytdeam_APO3B_nsvs_age_historand_nowarn.pdf",
#        plot = p1, width = 122, height = 91.5, units = "mm", useDingbats=FALSE)


blmmo_fix <- bglmer(data = sinaplotdf3, formula = no_apo_foci ~ APOBEC3B + numSVssc + donor_age_at_diagnosis + (1+xATM|histology),
                    family = "poisson", fixef.prior = normal, control = glmerControl(optimizer = "Nelder_Mead"))
blmmo_fix2 <- bglmer(data = sinaplotdf3, formula = no_apo_foci ~ xATM + APOBEC3B + numSVssc + donor_age_at_diagnosis + (1+xATM|histology),
                    family = "poisson", fixef.prior = normal, control = glmerControl(optimizer = "Nelder_Mead"))



blmmo_fix2
coef(blmmo_fix2)
summary(blmmo_fix2)
anova(blmmo_fix, blmmo_fix2)


# 
# blmmo_fix <- glmmTMB(data = sinaplotdf3, formula = no_apo_foci ~ APOBEC3B + numSVssc + donor_age_at_diagnosis + (1+xBRCA2|histology),
#                      family = "nbinom2", control = glmmTMBControl(profile=TRUE))
# 
# blmmo_fix2 <- glmmTMB(data = sinaplotdf3, formula = no_apo_foci ~ xBRCA2 + APOBEC3B + numSVssc + donor_age_at_diagnosis + (1+xBRCA2|histology),
#                      family = "nbinom2", control = glmmTMBControl(profile=TRUE))
# 
# summary(blmmo_fix2)
# anova(blmmo_fix, blmmo_fix2)




# 
# 
# 
# ### pull in APOBEC/AID expression
# apoexpr <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/data/RNASeq_joint_fpkm_uq_APOBECs.tsv", as.is = T)
# apoexpr$feature <- sapply(X = strsplit(x = apoexpr$feature, split = "::"), FUN = '[', 2)
# # apoexpr <- t(apoexpr)
# 
# 
# reltab <- read.delim(file = "/srv/shared/vanloo/ICGC_annotations/release_may2016.v1.4.tsv", as.is = T)
# # nchar(reltab$tumor_rna_seq_aliquot_id) > 36
# reltab <- reltab[reltab$tumor_rna_seq_aliquot_id != "", ]
# # check for potential sample number mismatches
# # reltab[reltab$tumor_wgs_specimen_count > 1, "tumor_rna_seq_aliquot_id"]
# wgs_aliq_ids <- strsplit(reltab$tumor_wgs_aliquot_id, split = ",")
# 
# reltab_redup <- reltab[rep(1:nrow(reltab), lengths(wgs_aliq_ids)),]
# # reltab_redup$tumor_wgs_icgc_specimen_id <- unlist(strsplit(reltab_redup$tumor_wgs_icgc_specimen_id, split = ","))
# reltab_redup$tumor_wgs_aliquot_id <- unlist(wgs_aliq_ids)
# reltab_redup$tumor_rna_seq_aliquot_id <- sapply(X = strsplit(reltab_redup$tumor_rna_seq_aliquot_id, split = ","), '[', 1)
# reltab <- reltab_redup[, c("tumor_wgs_aliquot_id", "tumor_rna_seq_aliquot_id")]
# 
# 
# remoddf <- do.call(rbind, lapply(X = gsub(pattern = "-", replacement = ".", x = reltab$tumor_rna_seq_aliquot_id), FUN = function(smpl, df) df[, grep(pattern = smpl, x = colnames(df))], df = apoexpr))
# # sapply(1:13, function(x, df) min(df[,x][df[,x]>0]), df = remoddf)
# remoddf <- as.data.frame(t(t(log2(remoddf + 1e-3)) - colMeans(log2(remoddf + 1e-3))))
# colnames(remoddf) <- apoexpr$feature
# remoddf$sample <- reltab$tumor_wgs_aliquot_id
# 
# 
# sinaplotdf3 <- merge(x = sinaplotdf2, y = remoddf, by.x = "sample", by.y = "sample", all.x = T)
# 
# 
# df <- sinaplotdf3
# blmmo_fix <- bglmer(formula = no_apo_foci ~ numSVssc + (1|histology),
#                     data = df, family = "poisson", fixef.prior = normal(cov = diag(1,2)), glmerControl(optimizer = "Nelder_Mead"))
# 
# 
# anovaresults2 <- lapply(X = colnames(remoddf)[-14], FUN = testvar, df = sinaplotdf3, blmmo_fix = blmmo_fix)
# anovaresults2_df <- as.data.frame(do.call(rbind, anovaresults2))
# 
# # anovaresults <- mclapply(X = tophits, FUN = testvar, df = sinaplotdf2, mc.preschedule = T, mc.cores = 15)
# # basicmodel <- data.frame(AIC = anovaresults[[1]]$AIC[1], BIC = anovaresults[[1]]$BIC[1], logLik = anovaresults[[1]]$logLik[1], pval = anovaresults[[1]]$`Pr(>Chisq)`[2])
# # anovaresults_df <- as.data.frame(do.call(rbind, anovaresults))
# # resultsdf <- do.call(rbind, lapply(X = anovaresults, FUN = function(x) data.frame(AIC = x$AIC[2], BIC = x$BIC[2], logLik = x$logLik[2], pval = x$`Pr(>Chisq)`[2])))
# anovaresults2_df$driver <- colnames(remoddf)[-14][!sapply(anovaresults2, is.null)]
# # resultsv
# anovaresults2_df[anovaresults2_df$Pvalue > 0, "Padj"] <- p.adjust(anovaresults2_df[anovaresults2_df$Pvalue > 0, "Pvalue"], method = "fdr")
# rownames(anovaresults2_df) <- NULL
# # resultsv_adj[resultsv_adj < .05]
# # resultsv_adj[resultsv_adj < .1]
# 
# write.table(x = anovaresults2_df, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190220_LMEMPois_no_apo_foci_vs_CytDeam_nsvs_historand.txt", quote = T, sep = "\t", col.names = T, row.names = F)
# 
# 
# res1 <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190220_LMEMPois_no_apo_foci_vs_CytDeam_nsvs_historand.txt", as.is = T)
# res2 <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190220_LMEMPois_no_apo_foci_vs_drivers_nsvs_historand.txt", as.is = T)
# res12 <- rbind(res2, res1)
# res12[res12$Pvalue > 0, "Padj"] <- p.adjust(res12[res12$Pvalue > 0, "Pvalue"], method = "fdr")
# 
# p1 <- ggplot(data = res12[res12$Pvalue > 0, ], mapping = aes(x = Estimate, y = abs(Zvalue), colour = Padj < .05)) + geom_point(data = res1[res1$Pvalue > 0, ], colour = "blue", size = 2, show.legend = F) + geom_point(show.legend = F)
# p1 <- p1 + geom_text(mapping = aes(label = gsub(pattern = "^x", replacement = "", x = ifelse(Padj < .05, driver, "")), color = driver %in% c("APOBEC", "xTP53")), nudge_y = 0.15, show.legend = F) + theme_minimal() + scale_colour_manual(values = c('TRUE' = "red", 'FALSE' = "grey"))
# p1 <- p1 + labs(x = "Coefficient estimate", y = "|Z-value|")
# p1
# 
# ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190220_LMEMPois_no_apo_foci_vs_drivers_nsvs_historand.pdf",
#        plot = p1, width = 122, height = 91.5, units = "mm", useDingbats=FALSE)
# 
# 
# 
# 
# #checking what best model is
# 
# blmmo3 <- bglmer(formula = no_apo_foci ~ xTP53 + APOBEC3B + donor_age_at_diagnosis + (1+xTP53|histology) + (1|inferred_sex),
#                  data = sinaplotdf3, family = "poisson", fixef.prior = normal(cov = diag(1,4))) #, glmerControl(optimizer = "optimx", optCtrl=list(method="nlminb")))
# 
# allFit(blmmo3)
# blmmo2 <- bglmer(formula = no_apo_foci ~ xTP53 + APOBEC3B + donor_age_at_diagnosis + (1+xTP53+xEPHA2|histology) + (1|inferred_sex),
#                  data = sinaplotdf3, family = "poisson", fixef.prior = normal(cov = diag(1,4)), glmerControl(optimizer = "Nelder_Mead"))
# 
# blmmo3
# blmmo2
# summary(blmmo3)
# coef(blmmo3)
# anova(blmmo3, blmmo2)
# View(sinaplotdf2[, c("sample", "histology", "no_apo_foci", "xTP53")])
# 
# 
# p1 <- ggplot(data = sinaplotdf3, mapping = aes(x = histology, y = no_apo_foci + 1, colour = xCDKN1B)) + geom_boxplot() + scale_y_log10()
# p1 <- p1 + theme(axis.text.x = element_text(angle = 90))
# p1
# 
# p1 <- ggplot(data = sinaplotdf3[sinaplotdf3$histology == "Kidney-RCC-Pap",], mapping = aes(x = xPBRM1, y = no_apo_foci + 1)) + scale_y_log10() + geom_jitter()
# p1 <- p1 + theme(axis.text.x = element_text(angle = 90))
# p1
# 
# p1 <- ggplot(data = sinaplotdf3, mapping = aes(x = donor_age_at_diagnosis, y = no_apo_foci + 1)) + scale_y_log10() + geom_density_2d() #+ facet_wrap(~ histology)
# # p1 <- p1 + theme(axis.text.x = element_text(angle = 90))
# p1
# 
# 
# 
# blmmo3 <- bglmer(formula = no_apo_foci ~ xTP53 + APOBEC3B + donor_age_at_diagnosis + hasCT + numSVssc + (1+xTP53|histology) + (1|inferred_sex),
#                  data = sinaplotdf3, family = "poisson", fixef.prior = normal(cov = diag(1,6)))
# blmmo2 <- bglmer(formula = no_apo_foci ~ APOBEC3B + donor_age_at_diagnosis + hasCT + numSVssc + (1+xTP53|histology) + (1|inferred_sex),
#                  data = sinaplotdf3, family = "poisson", fixef.prior = normal(cov = diag(1,5)), glmerControl(optimizer = "Nelder_Mead"))
# 
# allFit(blmmo3)
# 
# ### checks for CTT type:
# ctttypes <- c("Eso-AdenoCA", "Stomach-AdenoCA", "Biliary-AdenoCA", "Liver-HCC", "ColoRect-AdenoCA")
# sinaplotdf_ctt <- sinaplotdf3[which(sinaplotdf3$histology %in% ctttypes), ]
# cttdriv <- colnames(sinaplotdf_ctt[, -c(1:12)])[which(rowSums(sinaplotdf_ctt[, -c(1:12)]) > 0)]
#   
# 
# anovaresults3 <- lapply(X = cttdriv, FUN = testvar, df = sinaplotdf_ctt)
# anovaresults3_df <- as.data.frame(do.call(rbind, anovaresults3))
# 
# # anovaresults <- mclapply(X = tophits, FUN = testvar, df = sinaplotdf2, mc.preschedule = T, mc.cores = 15)
# # basicmodel <- data.frame(AIC = anovaresults[[1]]$AIC[1], BIC = anovaresults[[1]]$BIC[1], logLik = anovaresults[[1]]$logLik[1], pval = anovaresults[[1]]$`Pr(>Chisq)`[2])
# # anovaresults_df <- as.data.frame(do.call(rbind, anovaresults))
# # resultsdf <- do.call(rbind, lapply(X = anovaresults, FUN = function(x) data.frame(AIC = x$AIC[2], BIC = x$BIC[2], logLik = x$logLik[2], pval = x$`Pr(>Chisq)`[2])))
# anovaresults3_df$driver <- cttdriv[!sapply(anovaresults3, is.null)]
# # resultsv
# anovaresults3_df[anovaresults3_df$Pvalue > 0, "Padj"] <- p.adjust(anovaresults2_df[anovaresults3_df$Pvalue > 0, "Pvalue"], method = "fdr")
# rownames(anovaresults3_df) <- NULL
# # resultsv_adj[resultsv_adj < .05]
# # resultsv_adj[resultsv_adj < .1]
# 
# write.table(x = anovaresults3_df, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_LMEMPois_no_apo_foci_vs_drivers_agegroup_sex_histol_CTT.txt", quote = T, sep = "\t", col.names = T, row.names = F)
# 
# 
# 
# 
# 
# ##### extra stuff
# 
# 
# colnames(anovaresults_df) <- c("Estimate", "Stderr", "zvalue", "pval", "driver")
# 
# p1 <- ggplot(data = anovaresults_df[anovaresults_df$Warn == "" & anovaresults_df$Err == "", ], mapping = aes(x = Estimate, y = abs(Zvalue), colour = Padj < .05)) + geom_point(show.legend = F)
# p1 <- p1 + geom_text(mapping = aes(label = gsub(pattern = "^x", replacement = "", x = ifelse(Padj < .05, driver, ""))), nudge_y = 0.15, show.legend = F) + theme_minimal() + scale_colour_manual(values = c('TRUE' = "red", 'FALSE' = "grey"))
# p1 <- p1 + labs(x = "Coefficient estimate", y = "|Z-value|")
# p1
# 
# ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/results/20190130_LMEMPois_no_apo_foci_vs_drivers_agegroup_sex_histol_10hitsmin.pdf",
#        plot = p1, width = 122, height = 91.5, units = "mm")
# 
#  
# # + xCCNE1 + xLINC00290 + xGNA13 + xPARK2 + 
# # x19ptelomere
# # xNF1     xSMARCA4 x19ptelomere       xGNA13
# # x19ptelomere
# 
# lmm3 <- glmer(formula = no_apo_foci ~ xTP53 + xCDKN2A + (1 + xTP53 + xCDKN2A|histology) + (1|donor_age_at_diagnosis) + (1|inferred_sex), data = sinaplotdf2, family = "poisson", control = glmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))
# lmm4 <- glmer(formula = no_apo_foci ~ xTP53 + (1 + xTP53 + xCDKN2A|histology) + (1|donor_age_at_diagnosis) + (1|inferred_sex), data = sinaplotdf2, family = "poisson", control = glmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))
# 
# blmm3 <- bglmer(formula = no_apo_foci ~ xTP53 + xCDKN2A + (1|histology) + (1|donor_age_at_diagnosis) + (1|inferred_sex),
#                 data = sinaplotdf2, family = "poisson", fixef.prior = normal(cov = diag(9,3)))
# 
# blmm4 <- bglmer(formula = no_apo_foci ~ xTP53 + (1|histology) + (1|donor_age_at_diagnosis) + (1|inferred_sex),
#                 data = sinaplotdf2, family = "poisson", fixef.prior = normal(cov = diag(9,2)))
# 
# anova(blmm3,blmm4)
# 
# 
# anova(lmm3,lmm4)
# 
# View(sinaplotdf2[, 1:10])
# 
# lmm1
# summary(lmm1)
# par(mfrow = c(2,2))
# plot(lmm1)
# coef(lmm1)
# 
# lmm2 <- lmer(formula = no_apo_foci ~ TERT + (1+histology) + (1|histology), data = sinaplotdf2[, -c(1,3,4)], REML = F)
# lmm3 <- lmer(formula = no_apo_foci ~ TP53 + SPOP + NF1 + TERT + SETD2 + (1|histology) + (1+histology), data = sinaplotdf2[, -c(1,3,4)], REML = F)
# 
# summary(lmm2)
# 
# anova(lmm2,lmm3)
# 
# checkdf <- sinaplotdf2[, c("sample", "histology", "no_foci", "no_sv_foci", "no_apo_foci", "inferred_sex", "xSPOP", "xTP53", "xTERT", "xNF1", "xFOXA1")]
# checkdf$cohort <- histology_all[match(checkdf$sample, histology_all$samplename), "projectcode"]
# p1 <- ggplot(data = checkdf[checkdf$histology == "Prost-AdenoCA",], mapping = aes(x = histology, y = no_apo_foci, colour = SPOP)) + geom_boxplot()
# p1
# 
# p1 <- ggplot(data = sinaplotdf2[sinaplotdf2$histology == "Prost-AdenoCA",], mapping = aes(x = xERG, y = no_apo_foci, colour = xERG)) + geom_jitter() + scale_y_log10()
# p1
# 
# 
p1 <- ggplot(data = sinaplotdf3, mapping = aes(x = histology, y = no_apo_foci + 1, colour = xBRCA2)) + geom_boxplot() + scale_y_log10() + geom_jitter()
p1 <- p1 + theme(axis.text.x = element_text(angle = 90))
p1
# 
# p1 <- ggplot(data = sinaplotdf2, mapping = aes(x = round(donor_age_at_diagnosis/10), y = no_apo_foci + 1)) + geom_point() + scale_y_log10()
# p1 <- p1 + theme(axis.text.x = element_text(angle = 90))
# p1
# 
# 
# 
# wilcox.test(x = checkdf[checkdf$histology == "Prost-AdenoCA" & checkdf$SPOP, "no_apo_foci"], y = checkdf[checkdf$histology == "Prost-AdenoCA" & !checkdf$SPOP, "no_apo_foci"], paired = F, alternative = "greater")
# wilcox.test(x = checkdf[checkdf$histology == "Skin-Melanoma" & checkdf$TP53, "no_apo_foci"], y = checkdf[checkdf$histology == "Skin-Melanoma" & !checkdf$TP53, "no_apo_foci"], paired = F, alternative = "two.sided")
# wilcox.test(x = checkdf[checkdf$histology == "Skin-Melanoma" & checkdf$TERT, "no_apo_foci"], y = checkdf[checkdf$histology == "Skin-Melanoma" & !checkdf$TERT, "no_apo_foci"], paired = F, alternative = "two.sided")
# wilcox.test(x = checkdf[checkdf$histology == "Bladder-TCC" & checkdf$TERT, "no_apo_foci"], y = checkdf[checkdf$histology == "Bladder-TCC" & !checkdf$TERT, "no_apo_foci"], paired = F, alternative = "two.sided")
# 
# 
# # head(drivers)
# # samplesbydriv <- by(data = drivers, INDICES = drivers$gene, FUN = function(x) unique(x$sample))
# # # emptyvec <- rep("",  max(lengths(samplesbydriv)))
# # 
# samplesbydriv <- do.call(cbind, by(data = drivers, INDICES = drivers$gene, FUN = function(x) {y <- unique(x$sample); c(y, rep("", 900-length(y)))}))
# # write.table(x = samplesbydriv, file = "/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/data/samplesbydriver.gmx", quote = F, sep = "\t", row.names = F)
# 
# 
# 
# 
# ######
# # melalabels <- read.delim(file = "/srv/shared/vanloo/ICGC_annotations/icgc_melanoma_new_label.txt", as.is = T)
# # sinaplotdf2$newhist <- sinaplotdf2$histology
# # sinaplotdf2 <- merge(x = sinaplotdf2, y = melalabels[, c("icgc_aliquot", "subtype")], by.x = "sample", by.y = "icgc_aliquot", all.x = T)
# # 
# # p1 <- ggplot(data = sinaplotdf2[sinaplotdf2$histology == "Skin-Melanoma", ], mapping = aes(x = subtype, y = no_sv_foci_apo))  + geom_boxplot() + geom_jitter()
# # p1
# 
# blmmo3 <- bglmer(formula = no_apo_foci ~ xTP53 + APOBEC3B + (1+xTP53 + APOBEC3B|histology) + (1|donor_age_at_diagnosis_grouped) + (1|inferred_sex),
#                  data = sinaplotdf3, family = "poisson", fixef.prior = normal(cov = diag(9,3)))
# blmmo2 <- bglmer(formula = no_apo_foci ~ xTP53 + (1+xTP53+APOBEC3B|histology) + (1|donor_age_at_diagnosis_grouped) + (1|inferred_sex),
#                  data = sinaplotdf3, family = "poisson", fixef.prior = normal(cov = diag(9,2)))
# 
# blmmo3
# blmmo2
# summary(blmmo3)
# summary(blmmo2)
# anova(blmmo2, blmmo3)
# 
# p1 <- ggplot(data = sinaplotdf3, mapping = aes(x = no_apo_foci + 1, y = APOBEC3G)) + geom_point() + scale_x_log10() + stat_smooth(method = "lm")
# p1 <- p1 + theme(axis.text.x = element_text(angle = 90))
# p1
# 
# 
# temp <- by(data = sinaplotdf3[!sinaplotdf3$histology %in% c("Bone-Benign", "Bone-Epith", "Bone-Osteosarc", "SoftTissue-Leiomyo", "SoftTissue-Liposarc"), ],
#            INDICES = sinaplotdf3[!sinaplotdf3$histology %in% c("Bone-Benign", "Bone-Epith", "Bone-Osteosarc", "SoftTissue-Leiomyo", "SoftTissue-Liposarc"), "histology"], FUN = function(x) glm(formula = no_apo_foci ~ donor_age_at_diagnosis, data = x, family = "poisson"))
# 
# temp <- by(data = sinaplotdf3[!sinaplotdf3$histology %in% c("Bone-Benign", "Bone-Epith", "Bone-Osteosarc", "SoftTissue-Leiomyo", "SoftTissue-Liposarc"), ],
#            INDICES = sinaplotdf3[!sinaplotdf3$histology %in% c("Bone-Benign", "Bone-Epith", "Bone-Osteosarc", "SoftTissue-Leiomyo", "SoftTissue-Liposarc"), "histology"], FUN = function(x) glmRob(formula = no_apo_foci ~ donor_age_at_diagnosis, data = x, family = "poisson"))
# temp2 <- lapply(temp, summary)
# temp2
# 
# 
# glm <- glm(data = sinaplotdf3, formula = no_apo_foci ~ donor_age_at_diagnosis, family = "poisson")
# summary(glm)
# 
# p1 <- ggplot(data = sinaplotdf3, mapping = aes(x = donor_age_at_diagnosis, y = no_apo_foci)) + geom_point(alpha = .5) + stat_smooth(method = "glm", method.args = list(family = "poisson")) + facet_wrap(~histology, scales = "free_y")
# p1
