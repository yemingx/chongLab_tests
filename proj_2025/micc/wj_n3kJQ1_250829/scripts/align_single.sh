fastp='/research/xieyeming1/software/Miniconda/envs/xieyeming1/bin/fastp'
# bwa='/mnt/software/anaconda3/envs/Bioinfo/bin/bwa'
samtools="/research/xieyeming1/software/Miniconda/envs/xieyeming1/bin/samtools"
bedtools="/research/xieyeming1/software/Miniconda/envs/xieyeming1/bin/bedtools"
thread="5"
ref="/research/xieyeming1/db/genome/hg19/Sequence/bt2_index_hg19mt/genome.fa"
macs2="/research/zhangchen/software/anaconda3/envs/BDscRNA_py27/bin/macs2"
# gatk="/mnt/Software/anaconda3/envs/Bioinfo/bin/gatk"
picard="/mnt/software/anaconda3/envs/Bioinfo/bin/picard"
#/research/xieyeming1/software/Miniconda/envs/cyclone/bin/picard
bowtie2="/mnt/software/anaconda3/envs/Bioinfo/bin/bowtie2"
java="/mnt/software/anaconda3/envs/Bioinfo/bin/java"

date
sample="$1"
R1="/research/xieyeming1/proj_2025/wj_n3kJQ1_250829/raw_data/*/V350364939_L02_${sample}_1.fq.gz"
R2="/research/xieyeming1/proj_2025/wj_n3kJQ1_250829/raw_data/*/V350364939_L02_${sample}_2.fq.gz"

mkdir -p ${sample}
cd ${sample}

# ${fastp} -w ${thread} -i ${R1} -o ${sample}_R1.clean.fq -I ${R2} -O ${sample}_R2.clean.fq -h ${sample}.html -j ${sample}.json

# ${bowtie2} --very-sensitive --end-to-end --phred33 -I 10 -X 700 -p ${thread} -x ${ref} -1 ${sample}_R1.clean.fq -2 ${sample}_R2.clean.fq -S ${sample}.clean.sam &> ${sample}.clean.align.txt
# rm ${sample}_R*.clean.fq
# ${samtools} view -@ ${thread} -bS ${sample}.clean.sam > ${sample}.clean.bam
# rm ${sample}.clean.sam
# ${samtools} sort -@ ${thread} ${sample}.clean.bam -o ${sample}.sorted.bam
# ${samtools} index ${sample}.sorted.bam

# ${picard} CreateSequenceDictionary R=${ref}
# ${picard} MarkDuplicates INPUT=${sample}.sorted.bam METRICS_FILE=${sample}.sorted.bam.metrics OUTPUT=${sample}.dedup.bam REMOVE_DUPLICATES=true
# ${samtools} index ${sample}.dedup.bam

# ${macs2} callpeak -f BAM -t ${sample}.dedup.bam --keep-dup all -g hs -n ${sample} > ${sample}.peak

# ${bedtools} bamtobed -i ${sample}.dedup.bam > ${sample}.dedup.bam.bed

# bedGraphToBigWig='/research/zhangchen/software/bedGraphToBigWig'
# # convert dedup bam to bw
# ref_fai=/research/xieyeming1/db/genome/hg19/Sequence/bwa_index_hg19mt/genome.fa.fai
# mkdir -p tmp
# cat ${sample}.dedup.bam.bed |cut -f1-3|shuf|head -18000000|sort -T tmp -k1,1 -k2,2n|awk '{print $0"\t"1}' > ${sample}.dedup.bam.bed.downsampled

# $bedtools genomecov -bg -i ${sample}.dedup.bam.bed.downsampled -g ${ref_fai} > ${sample}.dedup.bam.bedgraph
# ${bedGraphToBigWig} ${sample}.dedup.bam.bedgraph ${ref_fai} ${sample}.dedup.bam.bw

# qc metrics
echo -e "sample\traw_reads\tclean_reads\taligned_reads\tdedup_reads\tpeak_num" > qc_metrics.txt
raw_reads=$(cat ${sample}.json|head|grep total_reads|cut -d ":" -f2|cut -f1 -d ','|awk '{print $1/2}') # divided by 2
clean_reads=$(cat ${sample}.clean.align.txt|awk '{print $1}'|head -1)
aligned_reads=$(cat ${sample}.sorted.bam.metrics|head -8|tail -1| awk '{print $4}')
dup_rate=$(cat ${sample}.sorted.bam.metrics|head -8|tail -1| awk '{print $10}')
dedup_reads=$(cat ${sample}.sorted.bam.metrics|head -8|tail -1| awk '{print $11}')
peak_num=$(cat ${sample}_peaks.narrowPeak|wc -l)

sample_lab=$(cat /research/xieyeming1/proj_2025/wj_n3kJQ1_250829/raw_data/coremail_fetch/tmp/id_mapping.txt|awk -v a=${sample} -F ',' '$1==a'|cut -d ',' -f2)
echo -e "${sample}\t${sample_lab}\t${raw_reads}\t${clean_reads}\t${aligned_reads}\t${dup_rate}\t${dedup_reads}\t${peak_num}"|tr [:blank:] \\t > qc_metrics.txt

bin=100
chrom_size="/research/xieyeming1/db/genome/hg19/Sequence/bwa_index_hg19mt/hg19.chrom.sizes.txt"
# ${bedtools} makewindows -g ${chrom_size} -w ${bin}|sort -T tmp -k1,1 -k2,2n > windows_${bin}.bed
# ${bedtools} coverage -a /research/xieyeming1/proj_2025/wj_n3kJQ1_250829/align/13/windows_${bin}.bed -b ${sample}.dedup.bam.bed.downsampled -sorted -counts|cut -f1-4 > ${sample}_${bin}.bdg
# cat /research/xieyeming1/proj_2025/MICC_paper/misc/super_enhancer/SE/SE_01_0112_SE_ele_hg19.bed|tail -n +2|cut -f1-3|sort -T tmp -k1,1 -k2,2n > SE_01_0112_SE_ele_hg19.bed
# ${bedtools} map -a SE_01_0112_SE_ele_hg19.bed -b ${sample}_${bin}.bdg -c 4 -o max > ${sample}_${bin}.SE.bdg

# zcat /research/xieyeming1/db/attribute_table/hg19_regulatory_element.bed.gz |grep promoter|cut -f1-3|sort -T tmp -k1,1 -k2,2n > promoter.bed
# ${bedtools} map -a promoter.bed -b ${sample}_${bin}.bdg -c 4 -o max > ${sample}_${bin}.promoter.bdg

date
# 18000000


#  /research/zhangchen/software/UCSC_application/bigWigToBedGraph GSM2902630_HEK293_ATAC_medium_depth_bio2_tech1.bw HEK293_ngs_atac_S2.bedgraph
#  cat HEK293_ngs_atac_S2.bedgraph|grep chr7|awk '$4 > 0.003' > hek293_ngs_atac_S2.chr7.0_003.bedgraph


#  ${samtools} view -@ ${thread} -bS ${sample}_S1.clean.sam.unique_map > ${sample}_S1.clean.bam
#  ${samtools} sort -@ ${thread} ${sample}_S1.clean.bam -o ${sample}_S1.sorted.bam
#  ${samtools} index ${sample}_S1.sorted.bam
#
#  ${gatk} MarkDuplicates --INPUT ${sample}_S1.sorted.bam --METRICS_FILE ${sample}_S1.sorted.bam.metrics --OUTPUT ${sample}_S1.dedup.bam  --REMOVE_DUPLICATES true
