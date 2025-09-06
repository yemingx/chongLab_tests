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
setwd('/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/test_proj/wj_n3kJQ1_250829/scripts')
options(bedtools.path = "/research/xieyeming1/software/Miniconda/envs/fyt_py311/bin/")

# Create barcode treatment mapping table
id_map <- data.frame(barcode = c(1, 2, 3, 4, 13, 14, 15),
  treatment = c("200nM_2h", "500nM_2h", "200nM_4h", "500nM_4h", "n3k_ctrl", "200nM_8h", "500nM_8h"))
id_map

align_dir='/research/xieyeming1/proj_2025/wj_n3kJQ1_250829/align'

bdg_files<-list.files(align_dir,'*.SE.bdg$',recursive = T,full.names = T)
bdg_files
bdg<-fread(bdg_files[1])
head(bdg)
dim(bdg)
peak_list2<-map(bdg_files,function(x){
    peak<-fread(x)
    colnames(peak)<-c('chr','start','end','maxBinCount_n3k')
    peak$barcode<-gsub(".*/(\\d+)/\\d+_.*", "\\1", x)
    peak$treatment<-id_map$treatment[match(peak$barcode,id_map$barcode)]
    peak$peak_id<-paste0(peak$chr,':',peak$start,'-',peak$end)
    return_table<-peak[,c('peak_id','treatment','maxBinCount_n3k')]
    return(return_table)
})

peak_compare<-do.call(rbind,peak_list2)
head(peak_compare)
table(peak_compare$treatment)

peak_compare_dt<-as.data.table(peak_compare)
peak_compare_dt$maxBinCount_n3k<-as.numeric(peak_compare_dt$maxBinCount_n3k)

# convert peak_compare_dt to wide format
wide_data <- dcast(peak_compare_dt, peak_id ~ treatment, value.var = "maxBinCount_n3k")
head(wide_data)
# col2 to col7 subtracted by col8
wide_data<-as.data.frame(wide_data)
wide_data[,2:8]<-(wide_data[,2:8]+0.1)/(wide_data[,8]+0.1)
wide_data<-wide_data[,-8]
wide_data<-as.data.table(wide_data)
# convert wide_data to long format
long_data<-melt(wide_data,id.vars = 'peak_id',variable.name = 'treatment',value.name = 'maxBinCount_n3k')
head(long_data)

# use long_data to plot a density plot panel, each graph is a treatment, for each graph, x is log2(maxBinCount_n3k)
p<-ggplot(long_data,aes(x=log2(maxBinCount_n3k),fill=treatment))+
    geom_density(alpha=0.3)+
    theme_bw()+ 
    scale_fill_npg()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12))+ # increase y x axis text lab size+
    theme_classic()+# add abline of x = 0
    geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
    labs(x='log2(treatment/ctrl)',y='Density')+# add panel title
    ggtitle(paste0('wj_n3kJQ1_250829\nregion: hek, super_enhancer, N=',dim(bdg)[1],'\nsignal: n3k maxBinCount, bin_size=50bp') )+
    facet_wrap(~treatment)+
    theme(strip.text = element_text(size = 14))
png('wj_n3kJQ1_250829_n3k_maxBinCount_super_enhancer_density.png',width = 8,height = 4,units = 'in',res = 100)
print(p)
dev.off()

############################## promoter ##############################
bdg_files<-list.files(align_dir,'*.promoter.bdg$',recursive = T,full.names = T)
bdg_files
bdg<-fread(bdg_files[1])
head(bdg)
dim(bdg)
peak_list2<-map(bdg_files,function(x){
    peak<-fread(x)
    colnames(peak)<-c('chr','start','end','maxBinCount_n3k')
    peak$barcode<-gsub(".*/(\\d+)/\\d+_.*", "\\1", x)
    peak$treatment<-id_map$treatment[match(peak$barcode,id_map$barcode)]
    peak$peak_id<-paste0(peak$chr,':',peak$start,'-',peak$end)
    return_table<-peak[,c('peak_id','treatment','maxBinCount_n3k')]
    return(return_table)
})

peak_compare<-do.call(rbind,peak_list2)
head(peak_compare)
table(peak_compare$treatment)

peak_compare_dt<-as.data.table(peak_compare)
peak_compare_dt$maxBinCount_n3k<-as.numeric(peak_compare_dt$maxBinCount_n3k)

# convert peak_compare_dt to wide format
wide_data <- dcast(peak_compare_dt, peak_id ~ treatment, value.var = "maxBinCount_n3k")
head(wide_data)
# col2 to col7 subtracted by col8
wide_data<-as.data.frame(wide_data)
wide_data[,2:8]<-(wide_data[,2:8]+0.1)/(wide_data[,8]+0.1)
wide_data<-wide_data[,-8]
wide_data<-as.data.table(wide_data)
# convert wide_data to long format
long_data<-melt(wide_data,id.vars = 'peak_id',variable.name = 'treatment',value.name = 'maxBinCount_n3k')
head(long_data)

# use long_data to plot a density plot panel, each graph is a treatment, for each graph, x is log2(maxBinCount_n3k)
p<-ggplot(long_data,aes(x=log2(maxBinCount_n3k),fill=treatment))+
    geom_density(alpha=0.3)+
    theme_bw()+ 
    scale_fill_npg()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12))+ # increase y x axis text lab size+
    theme_classic()+# add abline of x = 0
    geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
    labs(x='log2(treatment/ctrl)',y='Density')+# add panel title
    ggtitle(paste0('wj_n3kJQ1_250829\nregion: hek, promoter, N=',dim(bdg)[1],'\nsignal: n3k maxBinCount, bin_size=50bp') )+
    facet_wrap(~treatment) +
    theme(strip.text = element_text(size = 14))
png('wj_n3kJQ1_250829_n3k_maxBinCount_promoter_density.png',width = 8,height = 4,units = 'in',res = 100)# make image wider
print(p)
dev.off()

############################## end 250905 ##############################
