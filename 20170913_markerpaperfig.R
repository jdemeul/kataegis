## Generating rainfall plot zoom-in cascade

finalsigspersample <- by(data = katres_jmb_repsamples, INDICES = katres_jmb_repsamples$sample, FUN = function(x) unique(x$sigsum_SV))
finalsigspersample[1]

allthreekats <- unlist(lapply(X = finalsigspersample, FUN = function(x) length(intersect(c("APOBEC_SV", "alt-CDA_SV", "C[T>N]T_SV"), x)) == 3))
which(allthreekats)

# save(alldata, file = "075fc96d-6742-4ef3-9369-482592ad3a2f_alldata.Rdata")
# 78100212-65aa-4365-8b64-4b33f77732d5
p1 <- plot_rainfall_local(mutations = alldata[["snv"]], sample = SAMPLE, outdir = OUTDIR, svbreaks = alldata[["sv"]], chr = 1, bsgenome = genome, locus = c(242e6, 244e6), plot_text = F)
p1
ggsave(filename = file.path(paste0(SAMPLE, "_rainfall_chr1-243mb_zoom.pdf")), plot = p1, width = 3, height = 3)

temp <- get_trinuc_context(mutations = alldata[["snv"]], size = 1)

### modified rainfal function code

# plot the intermutation distance across the genome
plot_rainfall_local <- function(mutations, sample = "", outdir = "", svbreaks, bsgenome = genome, chr = NA, locus = NA, plot_text = F) {
  mutations <- mutations[mutations$chr == chr, ]
  mutations <- mutations[order(mutations$pos), ]
  mutations$trinuc <- as.character(get_trinuc_context(mutations = mutations))

  mutations$limd <- log10(do.call(c, by(data = mutations, INDICES = mutations$chr, FUN = function(x) c(diff(x$pos), NA))))
  mutations$color <- ifelse(mutations$ref == "C", ifelse(mutations$alt == "A", "C>A", ifelse(mutations$alt == "G", "C>G", "C>T") ),
                            ifelse(mutations$ref == "G", ifelse(mutations$alt == "T", "C>A", ifelse(mutations$alt == "C", "C>G", "C>T") ),
                                   ifelse(mutations$ref == "T", ifelse(mutations$alt == "A", "T>A", ifelse(mutations$alt == "C", "T>C", "T>G") ),
                                          ifelse(mutations$alt == "T", "T>A", ifelse(mutations$alt == "G", "T>C", "T>G")))))
  
  mutations$color <- factor(mutations$color, levels = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
  mutations$strand <- factor(ifelse(mutations$ref %in% c("C", "T"), "+", "-"), levels = c("+", "-"))
  
  # mutations$cumulpos <- mutations$pos + cumulative_genomic_pos[mutations$chr]
  # mutations$avpos <- do.call(c, by(data = mutations$pos, INDICES = mutations$chr,
                                        # FUN = function(x) filter(x, rep(0.5, 2), sides=2)))
    
    # breaks <- cumsum(c(0, as.numeric(seqlengths(genome)[c(1:22, "X", "Y")])))
    # label_positions <- filter(breaks, rep(.5, 2), sides = 2)[-length(breaks)]
    # annotations <- data.frame(lab_pos = label_positions, label = c(1:22, "X", "Y"))
    
    p1 <- ggplot(data = mutations, mapping = aes(x = pos, y = limd))
    if (plot_text) {
      p1 <- p1 + geom_text(mapping = aes(colour = color, label = trinuc), show.legend = F, alpha = .75, size = 2, na.rm=T)
    } else {
      p1 <- p1 + geom_point(mapping = aes(colour = color, shape = strand), show.legend = F, alpha = .75, size = 1, na.rm=T)
    }
    p1 <- p1 + scale_colour_manual(values = c("C>A" = "#377eb8", "C>G" = "#000000", "C>T" = "#e41a1c", "T>A" = "#984ea3", "T>C" = "#ffff33", "T>G" = "#4daf4a"))
    p1 <- p1 + theme_minimal() + theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5), 
                                       legend.position = "none", axis.ticks = element_blank()) + scale_y_continuous(limits = c(0,7))
    # p1 <- p1 + geom_vline(xintercept = breaks, color = "grey", alpha = 0.5)
    
    if (nrow(svbreaks) > 0) {
      svbreaks <- svbreaks[svbreaks$chrom == chr,]
      p1 <- p1 + geom_point(data = svbreaks, mapping = aes(x = start, y = 7), shape = "|")
    }
    
    # p1 <- p1 + geom_text(data = annotations, mapping = aes(x = lab_pos, y = -0.5, label = label), size = 3, na.rm = T)
    p1 <- p1 + labs(y = "log(intermutation distance)")
    if (all(!is.na(locus)))
      p1 <- p1 + scale_x_continuous(limits = locus)
  # ggsave(filename = file.path(outdir, paste0(sample, "_rainfall.png")), plot = p1, dpi = 300, width = 10, height = 3)
  # ggsave(filename = file.path(outdir, paste0(sample, "_rainfall.pdf")), plot = p1, width = 10, height = 3)
  return(p1)
}


## make plot of the SVs
library(ggforce)
svfile <- "/srv/shared/vanloo/ICGC-structural-variants/78100212-65aa-4365-8b64-4b33f77732d5.pcawg_consensus_1.6.161116.somatic.sv.bedpe.gz"
selectedsampleSVs <- read.delim(file = gzfile(svfile), header = T, sep = "\t", as.is = T)
selectedsampleSVs$cumstart1 <- selectedsampleSVs$start1 + cumsumdf[selectedsampleSVs$chrom1, "cumseqlengths"]
selectedsampleSVs$cumstart2 <- selectedsampleSVs$start2 + cumsumdf[selectedsampleSVs$chrom2, "cumseqlengths"]
vertical_positions <- c(TRA = 3, DUP = 2, DEL = 2, t2tINV = 1, h2hINV = 1)
ctrlpt_offset <- c(TRA = 1, DUP = -1, DEL = 1, t2tINV = 1, h2hINV = -1)*.5

transform_bedpe_to_bezier <- function(bedpesv, vpos, ctrlpt_offset) {
  row1 <- data.frame(x = bedpesv$cumstart1, y = vpos[[bedpesv$svclass]], sv_id = bedpesv$sv_id, svclass = bedpesv$svclass)
  row2 <- data.frame(x = bedpesv$cumstart1, y = vpos[[bedpesv$svclass]] + ctrlpt_offset[[bedpesv$svclass]], sv_id = bedpesv$sv_id, svclass = bedpesv$svclass)
  row3 <- data.frame(x = bedpesv$cumstart2, y = vpos[[bedpesv$svclass]] + ctrlpt_offset[[bedpesv$svclass]], sv_id = bedpesv$sv_id, svclass = bedpesv$svclass)
  row4 <- data.frame(x = bedpesv$cumstart2, y = vpos[[bedpesv$svclass]], sv_id = bedpesv$sv_id, svclass = bedpesv$svclass)
  outdf <- rbind(row1, row2, row3, row4, stringsAsFactors = F, deparse.level = 0)
  # colnames(outdf) <- c("x", "y", "sv_id", "svclass")
  return(outdf)
}

# debug(transform_bedpe_to_bezier)
selectedsampleSVs_plotdf <- do.call(rbind, lapply(X = split(x = selectedsampleSVs, f = 1:nrow(selectedsampleSVs)), FUN = transform_bedpe_to_bezier, vpos = vertical_positions, ctrlpt_offset = ctrlpt_offset))

p2 <- ggplot(data = selectedsampleSVs_plotdf) + geom_bezier(mapping = aes(x = x, y = y, group = sv_id, colour = svclass), alpha =.75, show.legend = F, size = .1)
p2 <- p2 + theme_minimal() + scale_color_manual(values = c(TRA = "#000000", DUP = "#8B4500", DEL = "#483D8B", t2tINV = "#53868B", h2hINV = "#76EE00"))
p2 <- p2 + theme(panel.grid = element_blank(), panel.grid.major.y = element_line(colour = "grey30", linetype = "dotted", size = .25), axis.title = element_blank(), axis.text = element_blank())
p2
ggsave(filename = "20180103_YilongPlot_above_rainfall.pdf", plot = p2, width = 10, height = 2)



## looking for alternative samples
finalsigspersample <- by(data = katres_jmb_repsamples, INDICES = katres_jmb_repsamples$sample, FUN = function(x) unique(x$sigsum_SV))
finalsigspersample[1]

allthreekats <- unlist(lapply(X = finalsigspersample, FUN = function(x) ("APOBEC" %in% x | "APOBEC_SV" %in% x) & ("alt-CDA_SV" %in% x | "alt-CDA" %in% x) & ("C[T>N]T_SV" %in% x | "C[T>N]T" %in% x)))
which(allthreekats)


#### plot CN
SAMPLE <- "fc950c33-faa4-0241-e040-11ac0c486786"

debug(plot_consensus_cn)
p1 <- plot_consensus_cn(tumourname = SAMPLE, chrominfo = cumsumdf)
ggsave(paste0(SAMPLE, "_cnplot_forpanelC.pdf"), plot = p1, width = 10, height = 2)


plot_consensus_cn <- function(tumourname, chrominfo) {
  
  consensus_cn_dir <- "/srv/shared/vanloo/ICGC_consensus_copynumber/20170119_release/"
  consensus_cn <- read.table(file = list.files(path = consensus_cn_dir, pattern = tumourname, full.names = T), header = T, sep = "\t", as.is = T)
  consensus_cn$cumstart <- cumsumdf[as.character(consensus_cn$chromosome), "cumseqlengths"] + consensus_cn$start
  consensus_cn$cumend <- cumsumdf[as.character(consensus_cn$chromosome), "cumseqlengths"] + consensus_cn$end
  
  p3 <- ggplot(data = consensus_cn) + theme_minimal()
  p3 <- p3 + scale_x_continuous(breaks = chrominfo[rownames(chrominfo) %in% unique(c(consensus_cn$chromosome, "Y")), "cumseqlengths"]/1000000)
  p3 <- p3 + scale_y_continuous(limits = c(0,5.1), label = function(x) format(x,nsmall = 0,scientific = FALSE))
  p3 <- p3 + geom_segment(mapping = aes(x = cumstart/1000000, xend = cumend/1000000, y = total_cn, yend = total_cn), colour = "goldenrod2", size = 2)
  p3 <- p3 + geom_segment(mapping = aes(x = cumstart/1000000, xend = cumend/1000000, y = minor_cn, yend = minor_cn), colour = "darkslategrey", size = 2)
  p3 <- p3 + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), axis.text.x = element_blank(), axis.title = element_blank())
}


