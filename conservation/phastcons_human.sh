#!/bin/bash
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=8:mem=8gb
#PBS -j oe
#PBS -N Phascons
#PBS -o Phascons.log

# === Activate environment ===
module load anaconda3/personal
source activate lecif_env

# === Paths ===
BW_PHASTCONS="/home/PhasCons/conserv_file/hg19.100way.phastCons.bw"
BED_DIR="/cuttag_dir/reciprocal_liftover"
OUT_DIR="${BED_DIR}/phastcons_scores"

mkdir -p "${OUT_DIR}"

# === Loop over all *_final_human_peaks_noXYrandom.bed files ===
for bed_file in ${BED_DIR}/*_final_human_peaks_noXYrandom.bed; do
    # Extract cell type from file name
    base_name=$(basename "${bed_file}")
    cell_type=$(echo "${base_name}" | cut -d'_' -f1)

    # Define output file
    output_file="${OUT_DIR}/${cell_type}_PhastCons_scores.txt"

    echo "🔍 Annotating ${cell_type} peaks with PhastCons..."

    # Run bigWigAverageOverBed
    bigWigAverageOverBed "${BW_PHASTCONS}" "${bed_file}" "${output_file}"
done

echo "✅ All PhastCons scoring complete. Output in: ${OUT_DIR}"



######################################################################################################## below part in R 


```{r}
library(ggplot2)
library(dplyr)

# === Config ===
cell_types <- c("FOXA2", "OLIG2", "NEUN")
input_dir <- "cuttag_dir/reciprocal_liftover"
score_dir <- "cuttag_dir/reciprocal_liftover/phastcons_scores"   # <- set this to where the annotated PhastCons BEDs are
output_dir <- "cuttag_dir/reciprocal_liftover/phastcons_scores" # <- set this to your desired output path
setwd(output_dir)

# === Function to process peaks per cell type ===
process_phastcons <- function(cell_type) {
  bed_file <- file.path(input_dir, paste0(cell_type, "_final_human_peaks_noXYrandom.bed"))
  score_file <- file.path(score_dir, paste0(cell_type, "_PhastCons_scores.txt"))
  
  # Load original BED
  bed <- read.table(bed_file, header=FALSE, stringsAsFactors=FALSE)
  colnames(bed)[1:4] <- c("chr", "start", "end", "Peak_ID")
  
  # Load PhastCons scores
  scores <- read.table(score_file, header=FALSE, stringsAsFactors=FALSE)
  colnames(scores) <- c("Peak_ID", "size", "covered", "sum", "mean0", "PhastCons")
  
  # Merge by Peak_ID
  merged <- merge(bed, scores[, c("Peak_ID", "PhastCons")], by="Peak_ID")
  
  # Remove unwanted chromosomes
  merged <- merged[!grepl("chrX|chrY|random", merged$chr), ]
  
  # Stratify into 3 bins
  merged <- merged %>%
    mutate(Conservation_Level = ntile(PhastCons, 3)) %>%
    mutate(Conservation_Level = factor(Conservation_Level, labels = c("Low", "Mid", "High")),
           Cell_Type = cell_type)

  # Write stratified LDSC-compatible BEDs
  for (level in levels(merged$Conservation_Level)) {
    subset <- merged[merged$Conservation_Level == level, ]
    out_file <- file.path(output_dir, paste0(cell_type, "_", tolower(level), "_conservation_hg19_LDSC_noXY.bed"))
    write.table(subset[, c("chr", "start", "end", "Peak_ID")], out_file,
                sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  }

  return(merged)
}

# === Run for all cell types ===
all_data <- bind_rows(lapply(cell_types, process_phastcons))

# === Plot  ===
pdf(file.path(output_dir, "PhastCons_stratification_summary.pdf"), height = 3, width=8)
ggplot(all_data, aes(x=Conservation_Level, y=PhastCons, fill=Conservation_Level)) +
  geom_boxplot(outlier.shape=NA) +
  facet_wrap(~Cell_Type, scales="free") +
  labs(title="PhastCons Conservation Stratification of Lifted Peaks",
       x="Conservation Level", y="PhastCons Score") +
  scale_fill_manual(values=c("#FFC5A6", "#DC8E90", "#A97882")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
```

######################################################################################################## below part in bash

for file in *_conservation_hg19_LDSC_noXY.bed; do
  sorted_file="${file%.bed}.hg38tohg19.sorted.bed"
  sort -k1,1 -k2,2n "$file" > "$sorted_file"
done
