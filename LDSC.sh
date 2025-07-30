#!/bin/bash

#PBS -N ldsc.log
#PBS -o ldsc.txt
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=15:mem=200gb

## Create conda env (DO THIS STEP BEFORE ATTEMPTING TO RUN LDSC) 
#conda env create --file environment.yml

module load anaconda3/personal
source activate ldsc

ANALYSIS_DATE="analysis_ID"


# List your cell types for analysis EG:
cell=('{bed1_file_in_hg19}' '{bed2_file_in_hg19}')

# List of diseases, only keep the disease you are interested in (select from the list below)
disease=('AD_Jan' 'PD_nalls_no23_r_py_2' 'PD_GP2' 'SCZ_PardinasInfo9' 'MS_Andlauer' 'ALS_rheenen')

LDSC="/rds/general/user/$USER/projects/epinott/live/scripts/ldsc/ldsc"
REQUIRED_FILES="/rds/general/user/$USER/projects/epinott/live/scripts/ldsc/required_files"
ANNOT_FILES="/rds/general/user/$USER/projects/epinott/live/user_analysed_data/Janis/ldsc/annot_file_final"
HERITABILITY_OUTPUT="/rds/general/user/$USER/projects/epinott/live/user_analysed_data/Janis/ldsc/heritability_final"
GWAS="/rds/general/user/$USER/projects/epinott/live/scripts/ldsc/required_files/GWAS_munge"
BED_LDSC="/rds/general/user/jlt121/projects/epinott/live/user_analysed_data/Janis/ldsc/bed_files/${ANALYSIS_DATE}"

##############################################
## Step 1: Creating annotation files (Parallel)
##############################################

> annot_commands.txt
for cell_type in "${cell[@]}"; do
  for chrom in {1..22}; do
    echo "python ${LDSC}/make_annot.py \
      --bed-file $BED_LDSC/${cell_type}.sorted.bed \
      --bimfile ${REQUIRED_FILES}/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom}.bim \
      --annot-file ${ANNOT_FILES}/${ANALYSIS_DATE}.${cell_type}.${chrom}.annot.gz" >> annot_commands.txt
  done
done

cat annot_commands.txt | parallel -j 15 --joblog ldsc_annot_joblog.txt

##################################################
## Step 2: Generating ld scores (Parallel)
##################################################

> l2_commands.txt
for cell_type in "${cell[@]}"; do
  for chrom in {1..22}; do
    echo "python ${LDSC}/ldsc.py \
      --l2 \
      --bfile ${REQUIRED_FILES}/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom} \
      --ld-wind-cm 1 \
      --annot ${ANNOT_FILES}/${ANALYSIS_DATE}.${cell_type}.${chrom}.annot.gz \
      --thin-annot \
      --out ${ANNOT_FILES}/${ANALYSIS_DATE}.${cell_type}.${chrom} \
      --print-snps ${REQUIRED_FILES}/list.txt" >> l2_commands.txt
  done
done

cat l2_commands.txt | parallel -j 15 --joblog ldsc_l2_joblog.txt


##################################################
## MODE1: Run in single annot model
##################################################

> ldsc_joint_model_commands.txt
for disease_type in "${disease[@]}"; do
  for cell_type in "${cell[@]}"; do
    echo "python ${LDSC}/ldsc.py \
      --h2 ${GWAS}/${disease_type}_ldsc.sumstats.gz \
      --ref-ld-chr ${REQUIRED_FILES}/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.,${ANNOT_FILES}/${ANALYSIS_DATE}.${cell_type}. \
      --out ${HERITABILITY_OUTPUT}/${disease_type}_${ANALYSIS_DATE}_${cell_type}_SINGLE \
      --overlap-annot \
      --frqfile-chr ${REQUIRED_FILES}/1000G_Phase3_frq/1000G.EUR.QC. \
      --w-ld-chr ${REQUIRED_FILES}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
      --print-coefficients" >> ldsc_single_model_commands.txt
  done
done

cat ldsc_single_model_commands.txt | parallel -j 15 --joblog ldsc_parallel_joblog.txt

exit 0


##################################################
## MODE2: Multiannotation model for comparison of unique contributions
##################################################

for disease_type in "${disease[@]}"
  do
  echo "Running ldsc for: ${disease_type}"

  # Construct the array for all the cell types 
  ref_ld_chr_args=()
  for cell_type in "${cell[@]}"
    do
      ref_ld_chr_args+=("${ANNOT_FILES}/${ANALYSIS_DATE}.${cell_type}.")
    done

  # Join array elements with commas
  ref_ld_chr_args_str=$(IFS=,; echo "${ref_ld_chr_args[*]}")

  # Remove trailing comma, if any
  ref_ld_chr_args_str=${ref_ld_chr_args_str%,}

  python ${LDSC}/ldsc.py \
  --h2 ${GWAS}/${disease_type}_ldsc.sumstats.gz \
  --ref-ld-chr ${REQUIRED_FILES}/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.,${ref_ld_chr_args_str} \
  --out ${HERITABILITY_OUTPUT}/${disease_type}_${ANALYSIS_DATE}_${SUB} \
  --overlap-annot  \
  --frqfile-chr ${REQUIRED_FILES}/1000G_Phase3_frq/1000G.EUR.QC. \
  --w-ld-chr ${REQUIRED_FILES}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
  --print-coefficients
  done

exit 0 

##################################################
## MODE3: Pairwise multiannotation model to decipher overlapping enrichments
##################################################

> ldsc_pairwise_model_commands.txt

for disease_type in "${disease[@]}"; do
for (( i=0; i<${#cell[@]}; i++ )); do
  for (( j=i; j<${#cell[@]}; j++ )); do
      cell_type1=${cell[$i]}
      cell_type2=${cell[$j]}
      echo "python ${LDSC}/ldsc.py \
        --h2 ${GWAS}/${disease_type}_ldsc.sumstats.gz \
        --ref-ld-chr ${REQUIRED_FILES}/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.,${ANNOT_FILES}/${ANALYSIS_DATE}.${cell_type1}.,${ANNOT_FILES}/${ANALYSIS_DATE}.${cell_type2}. \
        --out ${HERITABILITY_OUTPUT}/${disease_type}_${ANALYSIS_DATE}_${cell_type1}_vs_${cell_type2}_PAIRWISE \
        --overlap-annot \
        --frqfile-chr ${REQUIRED_FILES}/1000G_Phase3_frq/1000G.EUR.QC. \
        --w-ld-chr ${REQUIRED_FILES}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
        --print-coefficients" >> ldsc_pairwise_model_commands.txt
    done
  done
done

cat ldsc_pairwise_model_commands.txt | parallel -j 15 --joblog ldsc_pairwise_parallel_joblog.txt


exit 0

