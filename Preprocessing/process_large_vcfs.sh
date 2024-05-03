#!/bin/bash
#$ -l h_rt=72:00:00
#$ -t 1-39
#$ -pe smp 4
#$ -binding linear:4
#$ -l h_vmem=8G
#$ -o '/xchip/beroukhimlab/youyun/complexSV/code/Preprocessing/process_vcf_outputs'
#$ -e '/xchip/beroukhimlab/youyun/complexSV/code/Preprocessing/process_vcf_outputs'
#$ -N svaba_vcf2bedpe_annotate

set -e

export PATH="/xchip/beroukhimlab/youyun/miniconda3/bin:$PATH"

# 881 Samples
input_vcf=$(cat /xchip/beroukhimlab/youyun/complexSV/data/large_samples_with_both_sv_and_cna_rerun_20240417.txt | sed -n -e "$SGE_TASK_ID p")
sample_name=$(basename $input_vcf | cut -d '.' -f 1)
output_dir="/xchip/beroukhimlab/youyun/complexSV/data/TCGA/SV/bedpe/"
output_bedpe=$output_dir$sample_name".bedpe"
output_bedpe_clustered=$output_dir$sample_name".clustered"
output_bedpe_clustered_sv_clusters_and_footprints=$output_dir$sample_name".clustered.sv_clusters_and_footprints.tsv"
output_bedpe_clustered_sv_distance_pvals=$output_dir$sample_name".clustered.sv_distance_pvals"
output_bedpe_annotated=$output_dir$sample_name".annotated.bedpe"

# First run the vcf to bedpe which converts the file format and retain the annotation information in various fields
echo "Processing vcf to bedpe"
/xchip/beroukhimlab/youyun/miniconda3/bin/conda run -n rameen Rscript \
    /xchip/beroukhimlab/youyun/complexSV/code/Preprocessing/svaba_vcf2bedpe.r \
    -i $input_vcf -o $output_bedpe 

# Second run the cluster SV 
echo "Clustering SVs"
/xchip/beroukhimlab/youyun/miniconda3/bin/conda run -n rameen Rscript \
    /xchip/beroukhimlab/youyun/complexSV/code/ClusterSV/R/run_cluster_sv.R \
	-chr /xchip/beroukhimlab/youyun/complexSV/code/ClusterSV/references/hg19.chrom_sizes \
	-cen_telo /xchip/beroukhimlab/youyun/complexSV/code/ClusterSV/references/hg19_centromere_and_telomere_coords.txt \
	-bedpe $output_bedpe -n 4 \
    -out $output_bedpe_clustered

# Third run the annotation
echo "Annotating SVs"
/xchip/beroukhimlab/youyun/miniconda3/bin/conda run -n rameen Rscript \
    /xchip/beroukhimlab/youyun/complexSV/code/Preprocessing/annotate_bedpe.r \
    -i $output_bedpe -c $output_bedpe_clustered_sv_clusters_and_footprints \
    -o $output_bedpe_annotated

# Done
echo "Cleaning up"
rm $output_bedpe
rm $output_bedpe_clustered_sv_clusters_and_footprints
rm $output_bedpe_clustered_sv_distance_pvals

echo "Done"