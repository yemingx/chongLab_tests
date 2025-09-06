library(ggplot2)
library(dplyr)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap)
library("dendextend") 
library(data.table)
library(IRanges)
library(httpgd)
library(eulerr)
library(ggsci)
library(UpSetR)
library(ggsignif)

npg_colors <- pal_npg()(10)
hgd()
setwd('/research/xieyeming1/proj_2025/MICC_paper/peak_hg19/deeptools/tss_pausing')
options(bedtools.path = "/mnt/software/anaconda3/envs/R4_4/bin/")
MICC_1d<-fread('/research/xieyeming1/proj_2025/MICC_paper/peak_hg19/peak_overlap_1d/hekMiccEcoG_3lanes/hekMiccEcoG_3lanes.narrowpeak',sep='\t',header = F)
pausing<-fread('pausing.txt',sep='\t',header = T)
head(pausing)
colnames(pausing)<-c('transcript','chr','start','end','strand','Length','Copies','annotation','pausing_ratio','Promoter_reads','GeneBody_reads')
promoter_w<-200

################### paused/noPaused tss pausing index ###################
percentile_threshold <- quantile(pausing$GeneBody_reads, 0.25)
pause_thres <- quantile(pausing$pausing_ratio, 0.25)
pause_thres
ggplot(pausing, aes(x=log10(GeneBody_reads))) + geom_density() + theme_bw() +
 geom_vline(xintercept = log10(percentile_threshold), color = "red")
ggplot(pausing, aes(x=log10(pausing_ratio))) + geom_density() + theme_bw() + 
geom_vline(xintercept = log10(pause_thres), color = "red")+xlim(-1,3)

noPaused<-pausing[pausing$pausing_ratio<pause_thres&pausing$GeneBody_reads>percentile_threshold,]
paused<-pausing[pausing$pausing_ratio>pause_thres&pausing$GeneBody_reads>percentile_threshold,]
head(paused)
paused_bed<-paused[,c('chr','start','end','transcript','pausing_ratio','strand','GeneBody_reads','Promoter_reads')]
paused_tss_bed<-paused_bed
paused_tss_bed$end<-ifelse(paused_tss_bed$strand=='+',paused_tss_bed$start+promoter_w,paused_tss_bed$end+promoter_w)
paused_tss_bed$start<-ifelse(paused_tss_bed$strand=='+',paused_tss_bed$start-promoter_w,paused_tss_bed$end-promoter_w)
paused_tss_bed<-paused_tss_bed[order(paused_tss_bed$chr,as.numeric(paused_tss_bed$start)),]

noPaused_bed<-noPaused[,c('chr','start','end','transcript','pausing_ratio','strand','GeneBody_reads','Promoter_reads')]
noPaused_tss_bed<-noPaused_bed
noPaused_tss_bed$end<-ifelse(noPaused_tss_bed$strand=='+',noPaused_tss_bed$start+promoter_w,noPaused_tss_bed$end+promoter_w)
noPaused_tss_bed$start<-ifelse(noPaused_tss_bed$strand=='+',noPaused_tss_bed$start-promoter_w,noPaused_tss_bed$end-promoter_w)
noPaused_tss_bed<-noPaused_tss_bed[order(noPaused_tss_bed$chr,as.numeric(noPaused_tss_bed$start)),]

################### paused/noPaused tss overlapping with Micc1D ###################
multi_overlap<-bedtoolsr::bt.multiinter(list(paused_tss_bed,noPaused_tss_bed,MICC_1d),cluster=TRUE)
head(multi_overlap)
multi_overlap_mat<-multi_overlap[multi_overlap$V5!='1,2,3',6:length(colnames(multi_overlap))]
cols<-c('paused_tss_bed','noPaused_tss_bed','MICC_1d')
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
pdf('euler_hekMiccEcoG_3lanes.pdf')
print(p3)
dev.off()

upset_data <- multi_overlap_mat
p2<-upset( data = as.data.frame(upset_data),
  sets = cols,
  order.by = "freq",
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Set Size")
p2
pdf('upset_tss_pausing.pdf',width=5,height=5, onefile=F)
print(p2)
dev.off()
################### overlapped/no overlapped paused GeneBody_reads ###################
paused_tss_MiccOverlap<-bedtoolsr::bt.intersect(a=paused_tss_bed,b=MICC_1d,wa=TRUE,u=TRUE,sorted=TRUE)
paused_tss_NoMiccOverlap<-bedtoolsr::bt.intersect(a=paused_tss_bed,b=MICC_1d,wa=TRUE,sorted=TRUE,v=TRUE)
dim(paused_tss_MiccOverlap)
dim(paused_tss_NoMiccOverlap)

plot_title='paused_Micc1D_overlap_GeneBody_reads'
y_axis='log2_GeneBody_reads'
in_file1='paused_tss_MiccOverlap'
in_file2='paused_tss_NoMiccOverlap'
in_file_1 <- paused_tss_MiccOverlap[,c(1,2,3,7)]
in_file_1$subgroup<-in_file1
in_file_2 <- paused_tss_NoMiccOverlap[,c(1,2,3,7)]
in_file_2$subgroup<-in_file2

# combine all data
all_data <- rbind(in_file_1, in_file_2)
# all_data$feature_val <- log2(all_data$V7+0.1)
all_data$feature_val <- ifelse(all_data$V7 > 0, log2(all_data$V7), NA)
wilcox_result <- wilcox.test(in_file_2$V7, in_file_1$V7,alternative = "two.sided")
head(all_data)
summ <- all_data %>%
  group_by(subgroup) %>%
  dplyr::summarize(n = n(), mean = round(mean(feature_val, na.rm = TRUE),2),
    max_val = round(max(feature_val, na.rm = TRUE),2),
    sd = sd(feature_val, na.rm = TRUE),median = round(median(feature_val, na.rm = TRUE),2))
summ
levels(factor(all_data$subgroup))

p1<-ggplot(all_data, aes(x=subgroup, y=feature_val, fill=subgroup)) +
  geom_violin(trim=FALSE, bw=0.3, na.rm = TRUE) +
  geom_boxplot(width=0.1, outlier.shape = NA, na.rm = TRUE) +
  geom_text(aes(label = paste0('N=',n), y = max(max_val, na.rm = TRUE)), 
            data = summ, size=4, vjust = 2, hjust = 2) +
  geom_text(aes(label = paste0('mean=',mean), y = max(mean, na.rm = TRUE)), 
            data = summ, size=4, vjust = -2, hjust = 2) +
  scale_fill_manual(values=c("indianred","lightblue")) +
  theme_classic() +  geom_signif(comparisons = list(levels(factor(all_data$subgroup))), 
              test = "wilcox.test", map_signif_level = TRUE,
              y_position = max(all_data$feature_val, na.rm = TRUE) * 1.2) +
  labs(title = plot_title,
    subtitle = paste0("Wilcoxon p-value = ", format.pval(wilcox_result$p.value)),
    y = y_axis) +
  theme(axis.text = element_text(face='bold'),
    axis.title = element_text(face="bold"),plot.title = element_text(face="bold"),
    legend.position="top")+ylim(-5,10)
p1
pdf(paste0(plot_title,'.pdf'))
print(p1)
dev.off()

################### paused Promoter_reads overlapped/no overlapped  ###################
paused_tss_MiccOverlap<-bedtoolsr::bt.intersect(a=paused_tss_bed,b=MICC_1d,wa=TRUE,u=TRUE,sorted=TRUE)
paused_tss_NoMiccOverlap<-bedtoolsr::bt.intersect(a=paused_tss_bed,b=MICC_1d,wa=TRUE,sorted=TRUE,v=TRUE)
dim(paused_tss_MiccOverlap)
dim(paused_tss_NoMiccOverlap)

plot_title='paused_Promoter_reads_Micc1D_overlap'
y_axis='log2_Promoter_reads'
in_file1='paused_tss_MiccOverlap'
in_file2='paused_tss_NoMiccOverlap'
in_file_1 <- paused_tss_MiccOverlap[,c(1,2,3,8)]
colnames(in_file_1)<-c('V1','V2','V3','V4')
in_file_1$subgroup<-in_file1
in_file_2 <- paused_tss_NoMiccOverlap[,c(1,2,3,8)]
colnames(in_file_2)<-c('V1','V2','V3','V4')
in_file_2$subgroup<-in_file2

# combine all data
all_data <- rbind(in_file_1, in_file_2)
all_data$feature_val <- ifelse(all_data$V4 > 0, log2(all_data$V4), NA)
wilcox_result <- wilcox.test(in_file_2$V4, in_file_1$V4,alternative = "two.sided")
head(all_data)
summ <- all_data %>%
  group_by(subgroup) %>%
  dplyr::summarize(n = n(), mean = round(mean(feature_val, na.rm = TRUE),2),
    max_val = round(max(feature_val, na.rm = TRUE),2),
    sd = sd(feature_val, na.rm = TRUE),median = round(median(feature_val, na.rm = TRUE),2))
summ
levels(factor(all_data$subgroup))

p1<-ggplot(all_data, aes(x=subgroup, y=feature_val, fill=subgroup)) +
  geom_violin(trim=FALSE, bw=0.3, na.rm = TRUE) +
  geom_boxplot(width=0.1, outlier.shape = NA, na.rm = TRUE) +
  geom_text(aes(label = paste0('N=',n), y = max(max_val, na.rm = TRUE)), 
            data = summ, size=4, vjust = 2, hjust = 2) +
  geom_text(aes(label = paste0('mean=',mean), y = max(mean, na.rm = TRUE)), 
            data = summ, size=4, vjust = -2, hjust = 2) +
  scale_fill_manual(values=c("indianred","lightblue")) +
  theme_classic() +  geom_signif(comparisons = list(levels(factor(all_data$subgroup))), 
              test = "wilcox.test", map_signif_level = TRUE,
              y_position = max(all_data$feature_val, na.rm = TRUE) * 1.2) +
  labs(title = plot_title,
    subtitle = paste0("Wilcoxon p-value = ", format.pval(wilcox_result$p.value)),
    y = y_axis) +
  theme(axis.text = element_text(face='bold'),
    axis.title = element_text(face="bold"),plot.title = element_text(face="bold"),
    legend.position="top")+ylim(-5,10)
p1
pdf(paste0(plot_title,'.pdf'))
print(p1)
dev.off()
################### overlapped/no overlapped noPaused GeneBody_reads ###################
noPaused_tss_MiccOverlap<-bedtoolsr::bt.intersect(a=noPaused_tss_bed,b=MICC_1d,wa=TRUE,u=TRUE,sorted=TRUE)
noPaused_tss_NoMiccOverlap<-bedtoolsr::bt.intersect(a=noPaused_tss_bed,b=MICC_1d,wa=TRUE,sorted=TRUE,v=TRUE)
dim(noPaused_tss_MiccOverlap)
dim(noPaused_tss_NoMiccOverlap)

plot_title='noPaused_Micc1D_overlap_GeneBody_reads'
y_axis='log2_GeneBody_reads'
in_file2='noPaused_tss_NoMiccOverlap'
in_file1='noPaused_tss_MiccOverlap'
in_file_1 <- noPaused_tss_MiccOverlap[,c(1,2,3,7)]
in_file_1$subgroup<-in_file1
in_file_2 <- noPaused_tss_NoMiccOverlap[,c(1,2,3,7)]
in_file_2$subgroup<-in_file2

# combine all data
all_data <- rbind(in_file_1, in_file_2)
# all_data$feature_val <- log2(all_data$V7+0.1)
all_data$feature_val <- ifelse(all_data$V7 > 0, log2(all_data$V7), NA)
wilcox_result <- wilcox.test(in_file_2$V7, in_file_1$V7,alternative = "two.sided")
head(all_data)
summ <- all_data %>%
  group_by(subgroup) %>%
  dplyr::summarize(n = n(), mean = round(mean(feature_val, na.rm = TRUE),2),
    max_val = round(max(feature_val, na.rm = TRUE),2),
    sd = sd(feature_val, na.rm = TRUE),median = round(median(feature_val, na.rm = TRUE),2))
summ
levels(factor(all_data$subgroup))

p1<-ggplot(all_data, aes(x=subgroup, y=feature_val, fill=subgroup)) +
  geom_violin(trim=FALSE, bw=0.5, na.rm = TRUE) +
  geom_boxplot(width=0.1, outlier.shape = NA, na.rm = TRUE) +
  geom_text(aes(label = paste0('N=',n), y = max(max_val, na.rm = TRUE)), 
            data = summ, size=4, vjust = 2, hjust = 2) +
  geom_text(aes(label = paste0('mean=',mean), y = max(mean, na.rm = TRUE)), 
            data = summ, size=4, vjust = -2, hjust = 2) +
  scale_fill_manual(values=c("indianred","lightblue")) +
  theme_classic() +  geom_signif(comparisons = list(levels(factor(all_data$subgroup))), 
              test = "wilcox.test", map_signif_level = TRUE,
              y_position = max(all_data$feature_val, na.rm = TRUE) * 1.2) +
  labs(title = plot_title,
    subtitle = paste0("Wilcoxon p-value = ", format.pval(wilcox_result$p.value)),
    y = y_axis) +
  theme(axis.text = element_text(face='bold'),
    axis.title = element_text(face="bold"),plot.title = element_text(face="bold"),
    legend.position="top")+ylim(-5,10)
p1
pdf(paste0(plot_title,'.pdf'))
print(p1)
dev.off()

################### noPaused Promoter_reads overlapped/no overlapped  ###################
noPaused_tss_MiccOverlap<-bedtoolsr::bt.intersect(a=noPaused_tss_bed,b=MICC_1d,wa=TRUE,u=TRUE,sorted=TRUE)
noPaused_tss_NoMiccOverlap<-bedtoolsr::bt.intersect(a=noPaused_tss_bed,b=MICC_1d,wa=TRUE,sorted=TRUE,v=TRUE)
dim(noPaused_tss_MiccOverlap)
dim(noPaused_tss_NoMiccOverlap)

plot_title='noPaused_Promoter_reads_Micc1D_overlap'
y_axis='log2_Promoter_reads'
in_file1='noPaused_tss_MiccOverlap'
in_file2='noPaused_tss_NoMiccOverlap'
feature_col=8
args_bw=0.3

in_file_1 <- noPaused_tss_MiccOverlap[,c(1,2,3,feature_col)]
colnames(in_file_1)<-c('V1','V2','V3','V4')
in_file_1$subgroup<-in_file1
in_file_2 <- paused_tss_NoMiccOverlap[,c(1,2,3,feature_col)]
colnames(in_file_2)<-c('V1','V2','V3','V4')
in_file_2$subgroup<-in_file2

# combine all data
all_data <- rbind(in_file_1, in_file_2)
all_data$feature_val <- ifelse(all_data$V4 > 0, log2(all_data$V4), NA)
wilcox_result <- wilcox.test(in_file_2$V4, in_file_1$V4,alternative = "two.sided")
head(all_data)
summ <- all_data %>%
  group_by(subgroup) %>%
  dplyr::summarize(n = n(), mean = round(mean(feature_val, na.rm = TRUE),2),
    max_val = round(max(feature_val, na.rm = TRUE),2),
    sd = sd(feature_val, na.rm = TRUE),median = round(median(feature_val, na.rm = TRUE),2))
summ
levels(factor(all_data$subgroup))

p1<-ggplot(all_data, aes(x=subgroup, y=feature_val, fill=subgroup)) +
  geom_violin(trim=FALSE, bw=args_bw, na.rm = TRUE) +
  geom_boxplot(width=0.1, outlier.shape = NA, na.rm = TRUE) +
  geom_text(aes(label = paste0('N=',n), y = max(mean, na.rm = TRUE)), 
            data = summ, size=4, vjust = 0, hjust = 2) +
  geom_text(aes(label = paste0('mean=',mean), y = max(mean, na.rm = TRUE)), 
            data = summ, size=4, vjust = -2, hjust = 2) +
  scale_fill_manual(values=c("indianred","lightblue")) +
  theme_classic() +  geom_signif(comparisons = list(levels(factor(all_data$subgroup))), 
              test = "wilcox.test", map_signif_level = TRUE,
              y_position = max(all_data$feature_val, na.rm = TRUE) * 1.2) +
  labs(title = plot_title,
    subtitle = paste0("Wilcoxon p-value = ", format.pval(wilcox_result$p.value)),
    y = y_axis) +
  theme(axis.text = element_text(face='bold'),
    axis.title = element_text(face="bold"),plot.title = element_text(face="bold"),
    legend.position="top")+ylim(-5,10)
p1
pdf(paste0(plot_title,'.pdf'))
print(p1)
dev.off()
################### overlapped/no overlapped tss pausing index ###################
head(MICC_1d)
head(pausing)
pausing_bed<-pausing[,c('chr','start','end','transcript','pausing_ratio','strand','GeneBody_reads','Promoter_reads')]
pausing_tss_bed<-pausing_bed
pausing_tss_bed$end<-ifelse(pausing_tss_bed$strand=='+',pausing_tss_bed$start+promoter_w,pausing_tss_bed$end+promoter_w)
pausing_tss_bed$start<-ifelse(pausing_tss_bed$strand=='+',pausing_tss_bed$start-promoter_w,pausing_tss_bed$end-promoter_w)
pausing_tss_bed<-pausing_tss_bed[order(pausing_tss_bed$chr,as.numeric(pausing_tss_bed$start)),]

tss_MiccOverlap<-bedtoolsr::bt.intersect(a=pausing_tss_bed,b=MICC_1d,wa=TRUE,u=TRUE,sorted=TRUE)
tss_NoMiccOverlap<-bedtoolsr::bt.intersect(a=pausing_tss_bed,b=MICC_1d,wa=TRUE,sorted=TRUE,v=TRUE)
dim(tss_MiccOverlap)
dim(tss_NoMiccOverlap)

plot_title='Micc1D_overlap_TSS_pausing_index'
y_axis='log2_pausing_index'
in_file2='tss_NoMiccOverlap'
in_file1='tss_MiccOverlap'
in_file_1 <- tss_MiccOverlap[,c(1,2,3,5)]
in_file_1$subgroup<-in_file1
in_file_2 <- tss_NoMiccOverlap[,c(1,2,3,5)]
in_file_2$subgroup<-in_file2

# combine all data
all_data <- rbind(in_file_1, in_file_2)
all_data$feature_val <- ifelse(all_data$V5 > 0, log2(all_data$V5), NA)

wilcox_result <- wilcox.test(in_file_2$V5, in_file_1$V5, 
                           alternative = "two.sided")
head(all_data)
summ <- all_data %>%
  group_by(subgroup) %>%
  dplyr::summarize(n = n(), mean = round(mean(feature_val, na.rm = TRUE),2),
    max_val = round(max(feature_val, na.rm = TRUE),2),
    sd = sd(feature_val, na.rm = TRUE))
summ
levels(factor(all_data$subgroup))

p1<-ggplot(all_data, aes(x=subgroup, y=feature_val, fill=subgroup)) +
  geom_violin(trim=FALSE, bw=0.3, na.rm = TRUE) +
  geom_boxplot(width=0.1, outlier.shape = NA, na.rm = TRUE) +
  geom_text(aes(label = paste0('N=',n), y = max(max_val, na.rm = TRUE)), 
            data = summ, size=4, vjust = 2, hjust = 2) +
  geom_text(aes(label = paste0('mean=',mean), y = max(mean, na.rm = TRUE)), 
            data = summ, size=4, vjust = -2, hjust = 2) +
  scale_fill_manual(values=c("indianred","lightblue")) +
  theme_classic() +  geom_signif(comparisons = list(levels(factor(all_data$subgroup))), 
              test = "wilcox.test", map_signif_level = TRUE,
              y_position = max(all_data$feature_val, na.rm = TRUE) * 1.2) +
  labs(title = plot_title,
    subtitle = paste0("Wilcoxon p-value = ", format.pval(wilcox_result$p.value)),
    y = y_axis) +
  theme(axis.text = element_text(face='bold'),
    axis.title = element_text(face="bold"),plot.title = element_text(face="bold"),
    legend.position="top")+ylim(0,5)
p1
pdf(paste0(plot_title,'.pdf'))
print(p1)
dev.off()



################### overlapped/no overlapped tss pausing TATA box, cpg motif index at genebody at promoter ###################

