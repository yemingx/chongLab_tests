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
library(purrr)
npg_colors <- pal_npg()(10)
hgd()
setwd('/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/archive/wj_n3kJQ1_250829')
options(bedtools.path = "/research/xieyeming1/software/Miniconda/envs/fyt_py311/bin/")

chip_bed<-fread('/research/xieyeming1/db/hek293/histone_chip/hek293_chip_seq.bed.gz')
h3k4me1_bed<-chip_bed[chip_bed$subtype=='h3k4me1',]
h3k4me1_bed<-h3k4me1_bed[order(h3k4me1_bed$chr,as.numeric(h3k4me1_bed$start)),]
head(chip_bed)

Micc1D<-fread('/research/xieyeming1/proj_2025/MICC_paper/peak_hg19/peak_overlap_1d/hekMiccEcoG_3lanes/hekMiccEcoG_3lanes.narrowpeak',sep='\t',header = F)
head(Micc1D)
h3k4me1_bed_Micc1D<-bedtoolsr::bt.intersect(a=h3k4me1_bed,b=Micc1D,wa=TRUE,u=TRUE,sorted=TRUE)
head(h3k4me1_bed_Micc1D)
dim(h3k4me1_bed_Micc1D)

# Create barcode treatment mapping table
id_map <- data.frame(barcode = c(1, 2, 3, 4, 13, 14, 15),
  treatment = c("200nM_2h", "500nM_2h", "200nM_4h", "500nM_4h", "n3k_ctrl", "200nM_8h", "500nM_8h"))
id_map

align_dir='/research/xieyeming1/proj_2025/wj_n3kJQ1_250829/align'

bam_files<-list.files(align_dir,'*.dedup.bam.bed$',recursive = T,full.names = T)
bam_files
bam_files[1]
head(h3k4me1_bed_Micc1D)
dim(h3k4me1_bed_Micc1D)
head(bedtoolsr::bt.coverage(a=h3k4me1_bed_Micc1D,b=bam_files[1],mean=TRUE))
?bedtoolsr::bt.coverage


SE<-read.table('/research/xieyeming1/proj_2025/MICC_paper/misc/super_enhancer/SE/SE_01_0112_SE_ele_hg19.bed',sep='\t',header = T)
head(SE)
colnames(SE)[1:3]<-c('chr','start','end')
peak_list2<-map(bam_files,function(x){
    peak<-bedtoolsr::bt.coverage(a=SE[1:3],b=x,counts=TRUE)
    colnames(peak)<-c('chr','start','end','count_n3k')
    peak$barcode<-gsub(".*/(\\d+)/\\d+\\.dedup\\.bam\\.bed$", "\\1", x)
    peak$treatment<-id_map$treatment[match(peak$barcode,id_map$barcode)]
    # add peak number 
    peak$peak_id<-paste0(peak$chr,':',peak$start,'-',peak$end)
    return_table<-peak[,c('peak_id','treatment','count_n3k')]
    return(return_table)
})
# colnames(peak_compare)[3]<-'count_n3k'
peak_compare<-do.call(rbind,peak_list2)
head(peak_compare)
table(peak_compare$treatment)
# use dplyr for each treatment, normalize the count by the sum of count_n3k in each treatment
# peak_compare<-peak_compare %>%
#   group_by(treatment) %>%
#   mutate(count_n3k=count_n3k/sum(count_n3k))

# replace 0 value in col count_n3k with the min non-zero value
min_nonzero<-min(peak_compare$count_n3k[peak_compare$count_n3k!=0])
peak_compare$count_n3k[peak_compare$count_n3k==0]<-min_nonzero

peak_compare_dt<-as.data.table(peak_compare)
# convert peak_compare to matrix, each row is a peak_id, each col is a treatment
# peak_compare_mat<-dcast(peak_compare_dt,peak_id~treatment,value.var = 'count_n3k')
# peak_compare_mat<-as.data.frame(peak_compare_mat)
# col2 to col7 divided by col8 and remove col8
# peak_compare_mat[,2:7]<-peak_compare_mat[,2:7]/peak_compare_mat[,8]
# peak_compare_mat<-peak_compare_mat[,-8]
# peak_compare_mat<-as.data.table(peak_compare_mat)
# # convert to long format
# peak_compare_long<-melt(peak_compare_mat,id.vars = 'peak_id',variable.name = 'treatment',value.name = 'count_n3k')
# head(peak_compare_long)
# use ggplot to make a violin plot, x is treatment, y is count_n3k, color is treatment
ggplot(peak_compare_dt,aes(x=treatment,y=log10(count_n3k),fill=treatment))+
  geom_violin()+
  geom_boxplot(width=0.1)+
  theme_bw()+ 
  scale_fill_npg()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) # increase y x axis text lab size

#############################################
peak_files<-list.files(align_dir,'*_peaks.narrowPeak',recursive = T,full.names = T)
head(peak_files)

?bedtoolsr::bt.intersect

x<-peak_files[2]
peak<-fread(x,sep='\t',header = F)
colnames(peak)<-c('chr','start','end','name','score','strand','signalValue','pValue','qValue','peak')
peak$barcode<-gsub(".*/(\\d+)/\\d+_peaks\\.narrowPeak", "\\1",x)
head(peak)
peak$treatment<-id_map$treatment[match(peak$barcode,id_map$barcode)]
head(peak)
dim(peak)
peak<-peak[order(peak$chr,as.numeric(peak$start)),]
head(h3k4me1_bed_Micc1D)
dim(h3k4me1_bed_Micc1D)
bedtoolsr::bt.intersect(a=h3k4me1_bed_Micc1D[,1:3],b=peak[,1:3],wa=TRUE,sorted=TRUE)


# use func map to iterate over peak_files
peak_list<-map(peak_files,function(x){
  peak<-fread(x,sep='\t',header = F)
  head(peak)
  overlap<-bedtoolsr::bt.intersect(a=h3k4me1_bed_Micc1D[,1:3],b=peak[,1:3],sorted=TRUE)
#   print(overlap)
  return(overlap)
})

# bind all rows and unique
shared_peak<-do.call(rbind,peak_list)
head(shared_peak)
dim(shared_peak)
shared_peak<-shared_peak[!duplicated(shared_peak),]
dim(shared_peak)
colnames(shared_peak)<-c('chr','start','end')
# merge(a=shared_peak,b=peak,by=c('chr','start','end'))
shared_peak_merge <- merge(shared_peak, as.data.frame(peak), by = c('chr', 'start', 'end'))
return_table<-shared_peak_merge[,c(7,12)]
?merge
peak_list2<-map(peak_files,function(x){
    peak<-fread(x,sep='\t',header = F)
    colnames(peak)<-c('chr','start','end','name','score','strand','signalValue','pValue','qValue','peak')
    peak$barcode<-gsub(".*/(\\d+)/\\d+_peaks\\.narrowPeak", "\\1",x)
    peak$treatment<-id_map$treatment[match(peak$barcode,id_map$barcode)]
    shared_peak_merge <- merge(shared_peak, as.data.frame(peak), by = c('chr', 'start', 'end'),all.x = TRUE)
    head(shared_peak_merge)
    return_table<-shared_peak_merge[,c(7,12)]
    return(return_table)
})

peak_compare<-do.call(rbind,peak_list2)
peak_compare
dim(peak_compare)

################### prepare feature table ###################
# rna_exp<-fread('/research/xieyeming1/db/attribute_table/hg19_gene_exp/hg19_rna_seq.bed',sep='\t',header = F)
# head(rna_exp)
# dim(rna_exp)
# rna_exp$transcript<-gsub('\\.\\d+$','',rna_exp$V4)
# rna_exp$tx_length<-rna_exp$V3-rna_exp$V2

# rna_exp_clean<-rna_exp[,c(9,10,6,7)]
# colnames(rna_exp_clean)<-c('transcript','tx_length','rna_exp','gene_id')

MICC_1d<-fread('/research/xieyeming1/proj_2025/MICC_paper/peak_hg19/shared_peak_1d/hekMiccEcoG_3lanes/hekMiccEcoG_3lanes.narrowpeak',sep='\t',header = F)
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

