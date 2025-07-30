#!/bin/bash
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=8:mem=8gb
#PBS -j oe
#PBS -N reciprocal_lift
#PBS -o reciprocal_lift.log

# === Activate environment ===
module load anaconda3/personal
source activate ucsc_tools

# === CONFIG ===
CHAIN_HG38_TO_HG19="/liftovers/hg38ToHg19.over.chain.gz"
CHAIN_HG19_TO_HG38="/liftovers/hg19ToHg38.over.chain.gz"
MIN_MATCH=0.5
MAX_LENGTH=25000
OUTPUT_DIR="cuttag_dir/reciprocal_liftover"

# === Cell types to loop over ===
CELLTYPES=("FOXA2" "OLIG2" "NEUN")

# === Create output folder ===
mkdir -p "${OUTPUT_DIR}"

# === Loop over each cell type ===
for CELL in "${CELLTYPES[@]}"; do
  echo "🔄 Processing cell type: $CELL"

  INPUT_BED="cuttag_dir/peak_calling/${CELL}_confidence_peaks.sorted.bed"

  # Output file paths
  LIFTED_HG19="${OUTPUT_DIR}/${CELL}_lifted_hg19.bed"
  UNMAPPED1="${OUTPUT_DIR}/${CELL}_unmapped_to_hg19.bed"
  LIFTED_BACK_HG38="${OUTPUT_DIR}/${CELL}_lifted_back_hg38.bed"
  UNMAPPED2="${OUTPUT_DIR}/${CELL}_unmapped_back_hg38.bed"
  RECIPROCAL="${OUTPUT_DIR}/${CELL}_reciprocal_matches.bed"
  IDS="${OUTPUT_DIR}/${CELL}_reciprocal_ids.txt"
  RECIP_HG19="${OUTPUT_DIR}/${CELL}_human_reciprocal.bed"
  FINAL_OUT="${OUTPUT_DIR}/${CELL}_final_human_peaks.bed"
  FINAL_FILTERED="${OUTPUT_DIR}/${CELL}_final_human_peaks_noXYrandom.bed"

  # === Step 1: hg38 → hg19 ===
  liftOver -minMatch=${MIN_MATCH} ${INPUT_BED} ${CHAIN_HG38_TO_HG19} ${LIFTED_HG19} ${UNMAPPED1}

  # === Step 2: hg19 → hg38 ===
  liftOver -minMatch=${MIN_MATCH} ${LIFTED_HG19} ${CHAIN_HG19_TO_HG38} ${LIFTED_BACK_HG38} ${UNMAPPED2}

  # === Step 3: Intersect to find reciprocal matches ===
  bedtools intersect -f 0.5 -r -wa -wb -a ${INPUT_BED} -b ${LIFTED_BACK_HG38} > ${RECIPROCAL}

  # === Step 4: Extract peak IDs that reciprocally map ===
  cut -f 4 ${RECIPROCAL} | sort | uniq > ${IDS}

  # === Step 5: Filter lifted hg19 BED to retain only reciprocal peaks ===
  grep -Ff ${IDS} ${LIFTED_HG19} > ${RECIP_HG19}

  # === Step 6: Filter for peaks ≤ MAX_LENGTH ===
  awk 'BEGIN{OFS="\t"} {if ($3-$2 <= '${MAX_LENGTH}') print $0}' ${RECIP_HG19} > ${FINAL_OUT}

  # === Step 7: Remove chrX, chrY, and random chromosomes ===
  grep -Ev 'chrX|chrY|random' ${FINAL_OUT} > ${FINAL_FILTERED}

  echo "✅ [${CELL}] Done. Output: ${FINAL_FILTERED}"
done

echo "🎉 All cell types processed. Cleaned BEDs are in: ${OUTPUT_DIR}"

for file in *final_human_peaks_noXYrandom.bed; do
  sorted_file="${file%.bed}.hg38tohg19.sorted.bed"
  sort -k1,1 -k2,2n "$file" > "$sorted_file"
done



