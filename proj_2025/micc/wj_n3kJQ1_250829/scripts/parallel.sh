sample_id=($(ls -d /research/xieyeming1/proj_2025/wj_n3kJQ1_250829/raw_data/*/|grep -v fetch|sed 's/\///g'|rev|cut -d "-" -f1|rev))
date
echo ${sample_id[@]}

job_num=10

echo ${sample_id[@]}|sed 's/\ /\n/g'|parallel -j ${job_num} "sh align_single.sh {}"

echo -e "sample_id\tsample_lab\traw_reads\tclean_reads\taligned_reads\tdup_rate\tdedup_reads\tpeak_num"|tr [:blank:] \\t > qc_metrics_sum.txt
cat */qc_metrics.txt >> qc_metrics_sum.txt