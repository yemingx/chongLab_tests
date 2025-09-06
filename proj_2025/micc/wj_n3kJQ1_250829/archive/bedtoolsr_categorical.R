library(ggplot2)
library(dplyr)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap)
library("dendextend") 
library(data.table)
library(IRanges)
library(eulerr)
library(ggsci)
library(UpSetR)
library(ggsignif)
library(patchwork)
library(GenomicRanges)
library(BSgenome) 
library(Biostrings)
library(httpgd)
npg_colors <- pal_npg()(10)
hgd()
setwd('/research/xieyeming1/proj_2025/MICC_paper/peak_hg19/deeptools/tss_divergent')
options(bedtools.path = "/research/xieyeming1/software/Miniconda/envs/fyt_py311/bin/")

# Function to calculate motif density
calculate_motif_density <- function(bed_df,motif, genome="hg19") {
    # Read BED file
    regions <- bed_df
    colnames(regions)[1:3] <- c("chr", "start", "end")
    
    # Create GRanges object
    gr <- makeGRangesFromDataFrame(regions, keep.extra.columns=TRUE)
    
    # Get genome sequence (change according to your genome)
    genome_seq <- getBSgenome(paste0("BSgenome.Hsapiens.UCSC.", genome))
    
    # Get sequences for each region
    sequences <- getSeq(genome_seq, gr)
    
    # Count motif motifs (motif dinucleotides)
    motif_counts <- vcountPattern(motif, sequences)
    
    # Calculate density (motif per base)
    region_lengths <- width(gr)
    motif_density <- motif_counts / region_lengths
    
    # Add results to original data
    regions$motif_count <- motif_counts
    regions$motif_density <- motif_density
    
    return(regions)
}

# Usage example:
# result <- calculate_motif_density(bed,'CG')
# write.table(result, "output_with_motif_density.bed", sep="\t", row.names=FALSE, quote=FALSE)

plot_comparison <- function(overlap_data1, overlap_data2, 
                          plot_title = "Comparison", data1_lab='data1',data2_lab='data2',
                          y_axis = "log2_value",transform_type='log2',
                          feature_col = 8,
                          args_bw = 0.3,
                          fill_colors = c("indianred", "lightblue"),
                          y_limits_fold = c(0.5, 1.5)) {
  
    # Prepare input data
    in_file_1 <- overlap_data1[, c(1, 2, 3, feature_col)]
    colnames(in_file_1) <- c('V1', 'V2', 'V3', 'V4')
    in_file_1$subgroup <- data1_lab
    
    in_file_2 <- overlap_data2[, c(1, 2, 3, feature_col)]
    colnames(in_file_2) <- c('V1', 'V2', 'V3', 'V4')
    in_file_2$subgroup <- data2_lab

    # Combine data
    all_data <- rbind(in_file_1, in_file_2)

    # Transform feature column based on transform_type
    if (transform_type == 'log2') {
        all_data$feature_val <- ifelse(all_data$V4 > 0, log2(all_data$V4), NA)
    } else if (transform_type == 'log10') {
        all_data$feature_val <- ifelse(all_data$V4 > 0, log10(all_data$V4), NA)
    } else if (transform_type == 'NA') {
        all_data$feature_val <- all_data$V4
    }
  
    # Statistical test
    wilcox_result <- wilcox.test(in_file_2$V4, in_file_1$V4, alternative = "two.sided")
    y_limits = c(0, 0)
    y_limit_val_hi <- quantile(all_data$feature_val, 0.999, na.rm = TRUE)
    y_limits[2] <- ifelse(y_limit_val_hi < 0, y_limit_val_hi * y_limits_fold[1], y_limit_val_hi * y_limits_fold[2])
    y_limit_val_low <- quantile(all_data$feature_val, 0.001, na.rm = TRUE)
    y_limits[1] <- ifelse(y_limit_val_low < 0, y_limit_val_low * y_limits_fold[2], y_limit_val_low * y_limits_fold[1])
    y_pos <- ifelse(y_limit_val_hi < 0, y_limit_val_hi * 0.8, y_limit_val_hi * 1.1)

    # Summary stats
    summ <- all_data %>%
        group_by(subgroup) %>%
        dplyr::summarize(
        n = n(), 
        mean = round(mean(feature_val, na.rm = TRUE), 2),
        max_val = round(max(feature_val, na.rm = TRUE), 2),
        sd = sd(feature_val, na.rm = TRUE),
        median = round(median(feature_val, na.rm = TRUE), 2)
        )

    # Create plot
    p <- ggplot(all_data, aes(x = subgroup, y = feature_val, fill = subgroup)) +
        geom_violin(trim = FALSE, bw = args_bw, na.rm = TRUE) +
        geom_boxplot(width = 0.1, outlier.shape = NA, na.rm = TRUE) +
        geom_text(aes(label = paste0('N=', n), y = median), 
                data = summ, size = 4, vjust = 2, hjust = 1.5) +
        geom_text(aes(label = paste0('mean= ', mean), y = median), 
                data = summ, size = 4, vjust = 0, hjust = 1.5) +
        scale_fill_manual(values = fill_colors) + theme_classic() +
        geom_signif(comparisons = list(levels(factor(all_data$subgroup))), 
                    test = "wilcox.test", map_signif_level = TRUE,
                    y_position = y_pos,size=1) +
        labs(title = plot_title,
            subtitle = paste0("Wilcoxon p-value = ", format.pval(wilcox_result$p.value)),
            y = y_axis) +
        theme(axis.text = element_text(face = 'bold'),
            axis.title = element_text(face = "bold"),
            plot.title = element_text(face = "bold"),
            legend.position = "top") +
        ylim(y_limits[1], y_limits[2])
#   pdf(paste0(plot_title, '.pdf'))
#   print(p)
#   dev.off()
  return(p)
}

################### prepare feature table ###################
# rna_exp<-fread('/research/xieyeming1/db/attribute_table/hg19_gene_exp/hg19_rna_seq.bed',sep='\t',header = F)
# head(rna_exp)
# dim(rna_exp)
# rna_exp$transcript<-gsub('\\.\\d+$','',rna_exp$V4)
# rna_exp$tx_length<-rna_exp$V3-rna_exp$V2

# rna_exp_clean<-rna_exp[,c(9,10,6,7)]
# colnames(rna_exp_clean)<-c('transcript','tx_length','rna_exp','gene_id')

MICC_1d<-fread('/research/xieyeming1/proj_2025/MICC_paper/peak_hg19/peak_overlap_1d/hekMiccEcoG_3lanes/hekMiccEcoG_3lanes.narrowpeak',sep='\t',header = F)
pausing_<-fread('pausing.txt',sep='\t',header = T)
head(pausing_)
colnames(pausing_)<-c('transcript','chr','start','end','strand','Length','Copies','annotation','pausing_ratio','Promoter_reads','GeneBody_reads')
# keep the substring after last | in annotation
pausing_$annotation<-gsub('.*\\|','',pausing_$annotation)
table(pausing_$annotation)
dim(pausing_)
promoter_w<-200
pausing<-pausing_
# pausing<-merge(pausing_,rna_exp_clean,by='transcript')
# head(pausing[,c('chr','start','end','transcript')])
dim(pausing)
# add genebody cpg density
genebody_cpg<-calculate_motif_density(pausing[,c('chr','start','end','transcript')],motif='CG')
colnames(genebody_cpg)<-c('chr','start','end','transcript','genebody_cpg_count','genebody_CpG')

divergent<-fread('divergent_peaks_w_header.bed',sep='\t',header = T)
head(divergent)

################### paused/noPaused tss pausing index ###################
percentile_threshold <- quantile(pausing$GeneBody_reads, 0.25)
percentile_threshold
pause_thres <- quantile(pausing$pausing_ratio, 0.25)
pause_thres
p0<-ggplot(pausing, aes(x=log10(GeneBody_reads))) + geom_density() + theme_bw() +
 geom_vline(xintercept = log10(percentile_threshold), color = "red")
p1<-ggplot(pausing, aes(x=log10(pausing_ratio))) + geom_density() + theme_bw() + 
geom_vline(xintercept = log10(pause_thres), color = "red")+xlim(-1,3)

# library(gridExtra)
# grid.arrange(p0, p1, ncol = 2, widths = c(1, 1))
# Print or save the combined plot
combined_plot <- p1 + p0 + 
    plot_layout(ncol = 2, widths = c(1, 1))  # Adjust widths as needed
# Print or save the combined plot
print(combined_plot)

pausing_bed<-pausing[,c('chr','start','end','transcript','pausing_ratio','strand','GeneBody_reads','Promoter_reads','annotation')]
pausing_tss_bed<-pausing_bed
pausing_tss_bed$end<-ifelse(pausing_tss_bed$strand=='+',pausing_tss_bed$start+promoter_w,pausing_tss_bed$end+promoter_w)
pausing_tss_bed$start<-ifelse(pausing_tss_bed$strand=='+',pausing_tss_bed$start-promoter_w,pausing_tss_bed$end-promoter_w)
pausing_tss_bed<-pausing_tss_bed[order(pausing_tss_bed$chr,as.numeric(pausing_tss_bed$start)),]
dim(pausing_tss_bed)
head(pausing_tss_bed)
table(pausing_tss_bed$annotation)
promoter_cpg<-calculate_motif_density(pausing_tss_bed[,c('chr','start','end','transcript')],motif='CG')
colnames(promoter_cpg)<-c('chr','start','end','transcript','promoter_cpg_count','promoter_CpG')
promoter_TATA<-calculate_motif_density(pausing_tss_bed[,c('chr','start','end','transcript')],motif='TATA')
colnames(promoter_TATA)<-c('chr','start','end','transcript','promoter_TATA_count','promoter_TATA')

pausing_tss_bed<-merge(pausing_tss_bed,promoter_cpg[,c(4,5,6)],by='transcript')
pausing_tss_bed<-merge(pausing_tss_bed,promoter_TATA[,c(4,5,6)],by='transcript')
pausing_tss_bed<-merge(pausing_tss_bed,genebody_cpg[,c(4,5,6)],by='transcript')
length(colnames(pausing_tss_bed))
pausing_tss_bed<-pausing_tss_bed[,c(2,3,4,1,5:15)]
head(pausing_tss_bed)
dim(pausing_tss_bed)

paused_tss_bed<-pausing_tss_bed[pausing_tss_bed$pausing_ratio>pause_thres&pausing_tss_bed$GeneBody_reads>percentile_threshold,]
paused_tss_bed<-paused_tss_bed[order(paused_tss_bed$chr,as.numeric(paused_tss_bed$start)),]
dim(paused_tss_bed)
noPaused_tss_bed<-pausing_tss_bed[pausing_tss_bed$pausing_ratio<pause_thres&pausing_tss_bed$GeneBody_reads>percentile_threshold,]
noPaused_tss_bed<-noPaused_tss_bed[order(noPaused_tss_bed$chr,as.numeric(noPaused_tss_bed$start)),]
dim(noPaused_tss_bed)

################### paused/noPaused tss overlapping with Micc1D ###################
filtered_divergent <- divergent %>%
  filter(
 #   Normalized_Tag_Count > 1,          # Minimum read coverage
    # focus_ratio > 0.9,                  # How focused the signal is
   findPeaks_Score > 25,               # Overall peak quality score
    # Fold_Change_vs_Local > 2,           # Enrichment over local background
    `p-value_vs_Local` < 0.05,          # Statistical significance
 #   Clonal_Fold_Change > 1            # Distinguish from PCR artifacts
  )

multi_overlap<-bedtoolsr::bt.multiinter(list(paused_tss_bed,noPaused_tss_bed,filtered_divergent,MICC_1d),cluster=TRUE)
head(multi_overlap)
multi_overlap_mat<-multi_overlap[,6:length(colnames(multi_overlap))]
# multi_overlap_mat<-multi_overlap[multi_overlap$V5!='1,2,3',6:length(colnames(multi_overlap))]
cols<-c('paused_tss_bed','noPaused_tss_bed','filtered_divergent','MICC_1d')
colnames(multi_overlap_mat)<-cols
head(multi_overlap_mat)
combinations <- apply(multi_overlap_mat, 1, function(row) {
    paste(colnames(multi_overlap_mat)[row == 1], collapse = "&")})
counts <- table(combinations)

set_counts <- as.numeric(counts)
names(set_counts) <- names(counts)
fit2 <- euler(set_counts)

p3<-plot(fit2, quantities = TRUE,     fills = npg_colors[c(1,9,2)])
p3

# pdf('euler_hekMiccEcoG_3lanes.pdf')
# print(p3)
# dev.off()

upset_data <- multi_overlap_mat
p2<-upset( data = as.data.frame(upset_data),
  sets = cols,
  order.by = "freq",
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Set Size")
p2

# pdf('upset_tss_pausing.pdf',width=5,height=5, onefile=F)
# print(p2)
# dev.off()

################### paused vs Micc1D numeric feature ###################
paused_tss_MiccOverlap<-bedtoolsr::bt.intersect(a=paused_tss_bed,b=MICC_1d,wa=TRUE,u=TRUE,sorted=TRUE)
paused_tss_NoMiccOverlap<-bedtoolsr::bt.intersect(a=paused_tss_bed,b=MICC_1d,wa=TRUE,sorted=TRUE,v=TRUE)
dim(paused_tss_MiccOverlap)
dim(paused_tss_NoMiccOverlap)
head(paused_tss_MiccOverlap)
head(paused_tss_NoMiccOverlap)

feature_cols<-colnames(paused_tss_bed)
cols_vec<-c(5,7,8,11,13,15)

feature_col<-feature_cols[cols_vec]
feature_col

plot_list <- list()

for (i in seq_along(cols_vec)) {
    col_idx <- cols_vec[i]
    col_name <- feature_cols[col_idx]
    
    p <- plot_comparison(
        data1_lab = 'paused_tss_MiccOverlap',
        data2_lab = 'paused_tss_NoMiccOverlap',
        overlap_data1 = paused_tss_MiccOverlap,
        overlap_data2 = paused_tss_NoMiccOverlap,
        plot_title = paste0("paused_Micc1D_overlap_", col_name),
        y_axis = paste0("log2_", col_name),
        feature_col = col_idx,
        y_limits = c(0.5, 1.3),
        transform_type='log2',
        args_bw = 0.3,
        fill_colors = c("indianred", "lightblue")
    )
    plot_list[[i]] <- p
}

panel_plot <- wrap_plots(plot_list, ncol = 3, guides = "collect") & 
    theme(legend.position = "bottom")
print(panel_plot)

################### paused vs Micc1D categorical feature ###################
category_feature_col=9
data1_lab = 'paused_tss_MiccOverlap'
data2_lab = 'paused_tss_NoMiccOverlap'
plot_title = paste0("paused_Micc1D_overlap_", colnames(paused_tss_bed)[category_feature_col])
y_axis = colnames(paused_tss_bed)[category_feature_col]

paused_tss_MiccOverlap_<-as.data.frame(table(paused_tss_MiccOverlap$V9))
paused_tss_MiccOverlap_$subgroup<-data1_lab
paused_tss_MiccOverlap_$pct<-round(paused_tss_MiccOverlap_$Freq/sum(paused_tss_MiccOverlap_$Freq),3)

paused_tss_NoMiccOverlap_<-as.data.frame(table(paused_tss_NoMiccOverlap$V9))
paused_tss_NoMiccOverlap_$subgroup<-data2_lab
paused_tss_NoMiccOverlap_$pct<-round(paused_tss_NoMiccOverlap_$Freq/sum(paused_tss_NoMiccOverlap_$Freq),3)

all_data<-rbind(paused_tss_MiccOverlap_,paused_tss_NoMiccOverlap_)
all_data
# ... existing code ...

# Create horizontal bar plot
p <- ggplot(all_data, aes(x = Var1, y = pct, fill = subgroup)) +
    geom_col(position = position_dodge(), width = 0.7) +
    geom_text(aes(label = scales::percent(pct, accuracy = 0.1)), 
             position = position_dodge(width = 0.7), 
             hjust = -0.1, size = 3) +
    scale_fill_manual(values = c("indianred", "lightblue")) +
    labs(title = plot_title,
         x = y_axis,
         y = "Percentage") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.text.y = element_text(size = 12), 
          legend.position = "top") +
    coord_flip() +  # Flip to horizontal
    scale_y_continuous(labels = scales::percent, 
                      expand = expansion(mult = c(0, 0.1)))  # Add space for labels

print(p)

################### noPaused vs Micc1D numeric feature ###################
noPaused_tss_MiccOverlap<-bedtoolsr::bt.intersect(a=noPaused_tss_bed,b=MICC_1d,wa=TRUE,u=TRUE,sorted=TRUE)
noPaused_tss_NoMiccOverlap<-bedtoolsr::bt.intersect(a=noPaused_tss_bed,b=MICC_1d,wa=TRUE,sorted=TRUE,v=TRUE)
feature_cols<-colnames(noPaused_tss_bed)
cols_vec<-c(5,7,8,9,11,13,15,17)
feature_col<-feature_cols[cols_vec]
plot_list <- list()

for (i in seq_along(cols_vec)) {
    col_idx <- cols_vec[i]
    col_name <- feature_cols[col_idx]
    
    p <- plot_comparison(
        data1_lab = 'noPaused_tss_MiccOverlap',
        data2_lab = 'noPaused_tss_NoMiccOverlap',
        overlap_data1 = noPaused_tss_MiccOverlap,
        overlap_data2 = noPaused_tss_NoMiccOverlap,
        plot_title = paste0("noPaused_Micc1D_overlap_", col_name),
        y_axis = paste0("log2_", col_name),
        feature_col = col_idx,
        y_limits = c(-5, 10),
        args_bw = 0.5,
        fill_colors = c("indianred", "lightblue")
    )
    plot_list[[i]] <- p
}

panel_plot <- wrap_plots(plot_list, ncol = 4, guides = "collect") & 
    theme(legend.position = "bottom")
print(panel_plot)

################### tss vs Micc1D numeric feature ###################
pausing_tss_bed<-pausing_tss_bed[order(pausing_tss_bed$chr,as.numeric(pausing_tss_bed$start)),]

tss_MiccOverlap<-bedtoolsr::bt.intersect(a=pausing_tss_bed,b=MICC_1d,wa=TRUE,u=TRUE,sorted=TRUE)
tss_NoMiccOverlap<-bedtoolsr::bt.intersect(a=pausing_tss_bed,b=MICC_1d,wa=TRUE,sorted=TRUE,v=TRUE)

dim(tss_MiccOverlap)
dim(tss_NoMiccOverlap)
colnames(pausing_tss_bed)
head(pausing_tss_bed)
feature_cols<-colnames(pausing_tss_bed)
cols_vec<-c(5,7,8,9,11,13,15,17)
feature_col<-feature_cols[cols_vec]
plot_list <- list()

for (i in seq_along(cols_vec)) {
    col_idx <- cols_vec[i]
    col_name <- feature_cols[col_idx]
    
    p <- plot_comparison(
        data1_lab = 'allTss_MiccOverlap',
        data2_lab = 'allTss_NoMiccOverlap',
        overlap_data1 = tss_MiccOverlap,
        overlap_data2 = tss_NoMiccOverlap,
        plot_title = paste0("allTss_Micc1D_overlap_", col_name),
        y_axis = paste0("log2_", col_name),
        feature_col = col_idx,
        y_limits_fold = c(0.5, 1.3),
        args_bw = 0.5,
        fill_colors = c("indianred", "lightblue")
    )
    plot_list[[i]] <- p
}

panel_plot <- wrap_plots(plot_list, ncol = 4, guides = "collect") & 
    theme(legend.position = "bottom")
print(panel_plot)

