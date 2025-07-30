#!/bin/bash
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=8:mem=16gb
#PBS -N preprocessing_atac
#PBS -o preprocessing_atac.log
#PBS -j oe

#update - after cutoff analysis a p value of 1e-2 seems more appropriate, 1e-4 is too stringent for ATAC (nuclei and iPSC)


# Load necessary modules
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate deeptools_env

# Updateable variables
ID="nextflow_ID"
PROJECT="project_ID"
CELLS=("FOXA2" "NEUN" "OLIG") #cell types as specified in file name
GENOME="mm" #mm or hs
BLACKLISTED="/rds/general/user/$USER/projects/epinott/live/user_analysed_data/Janis/blacklisted_regions/mm10-blacklist.v2.bed" #ALWAYS UPDATE TO RESPECTIVE GENOME
CHR_SIZES="/rds/general/user/$USER/home/ref_genome/mm10.chrom.sizes" #update if to respective genom, for human: "/rds/general/user/$USER/projects/epinott/live/user_analysed_data/Aydan/chrom_sizes_ncbi_grch38/WholeGenomeFasta/genome.sizes", for mouse "/rds/general/user/$USER/projects/epinott/projects/epinott/live/Ref_genome/mouse/mm10.chrom.sizes" 
EFFECTIVE_GENOME_SIZE="230812534" #effective genome size mouse: 230812534 ; Human: 2701495761 https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html`
PVALUE_MACS2="1e-2"

# Directories
ANALYSIS_DIR="/rds/general/user/$USER/projects/epinott/live/user_analysed_data/Janis/atac_seq/analysed_data/$ID"
FILE_LIST_DIR="/rds/general/user/$USER/projects/epinott/live/user_analysed_data/Janis/atac_seq/file_names/$ID"
BAM_IN_DIR="/rds/general/user/jlt121/projects/epinott/live/Nextflow/atac_seq/analysed_data/$ID/bowtie2/merged_library"
BAM_OUT_DIR="$ANALYSIS_DIR"
MERGED_BAM="$ANALYSIS_DIR/merged_bams"
BEDGRAPHS="$ANALYSIS_DIR/bedgraphs"
PEAK_CALLING="$ANALYSIS_DIR/peak_calling"
FEATURECOUNTS="$ANALYSIS_DIR/featurecounts"
BIGWIG="$ANALYSIS_DIR/bigwigs"

# Create required directories
create_directories() {
    for DIR in "$MERGED_BAM" "$BEDGRAPHS" "$PEAK_CALLING" "$FEATURECOUNTS" "$BIGWIG"; do
        mkdir -p "$DIR"
        if [ ! -d "$DIR" ]; then
            echo "Error: Unable to create directory $DIR. Exiting script."
            exit 1
        fi
    done
}

# Validate the file list
validate_file_list() {
    if [ ! -f "${FILE_LIST_DIR}/file.txt" ]; then
        echo "Error: File list ${FILE_LIST_DIR}/file.txt not found. Exiting script."
        exit 1
    fi
}

# Step 1: Remove low-quality reads, deduplicate, remove blacklisted regions, calling peaks for individual samples
process_individual_samples() {
    for FILE in $(cat "${FILE_LIST_DIR}/file.txt"); do
        if [ ! -f "${BAM_IN_DIR}/${FILE}" ]; then
            echo "Warning: Input BAM file ${FILE} not found. Skipping."
            continue
        fi

        DEDUP_FILE="${FILE%_REP1.mLb.clN.sorted.bam}.dedup.bam"
        NEW_FILE="${FILE%_REP1.mLb.clN.sorted.bam}.dedup.q30.bam"
        FINAL_FILE="${FILE%_REP1.mLb.clN.sorted.bam}.dedup.q30.no_blacklist.bam"

      # Deduplicate using Picard
        picard MarkDuplicates -XX:ParallelGCThreads=8 I="${BAM_IN_DIR}/${FILE}" O="${BAM_OUT_DIR}/${DEDUP_FILE}" M="${BAM_OUT_DIR}/${DEDUP_FILE}.marked_dup_metrics.txt" REMOVE_DUPLICATES=true 

        # Filter BAM by quality
        samtools view -b -F 4 -q 30 "${BAM_OUT_DIR}/${DEDUP_FILE}" > "${BAM_OUT_DIR}/${NEW_FILE}"

  
        # Remove blacklisted regions
        bedtools intersect -v -abam "${BAM_OUT_DIR}/${NEW_FILE}" -b "${BLACKLISTED}" > "${BAM_OUT_DIR}/${FINAL_FILE}"
        samtools index "${BAM_OUT_DIR}/${FINAL_FILE}"

        # Validate final BAM file
        if [ ! -s "${BAM_OUT_DIR}/${FINAL_FILE}" ]; then
            echo "Error: Final BAM file ${FINAL_FILE} is empty. Exiting script."
            exit 1
        fi

        # MACS2 peak calling for individual samples
        macs2 callpeak \
            -t ${BAM_OUT_DIR}/${FINAL_FILE} \
            -q ${PVALUE_MACS2} \
            -f BAMPE \
            --nomodel \
            --shift -37 \
            --extsize 73 \
            --cutoff-analysis \
            -g ${GENOME} \
            --keep-dup all \
            -n ${FINAL_FILE}.narrow.${PVALUE_MACS2}.macs2 \
            --outdir ${PEAK_CALLING}

        # Cleanup temporary files
        rm -f "${BAM_OUT_DIR}/${NEW_FILE}"
    done
}

# Step 2. Merge the BAMs per cell type
merge_bams_per_cell_type() {
    cd "${BAM_OUT_DIR}"

    for CELL in "${CELLS[@]}"; do
        file="${PROJECT}_${CELL}_all"

        echo "Merging files for ${CELL}:"
        ls "${BAM_OUT_DIR}"/*"${CELL}"*.dedup.q30.no_blacklist.bam

        samtools merge "${MERGED_BAM}/${file}.bam" $(ls "${BAM_OUT_DIR}"/*"${CELL}"*.dedup.q30.no_blacklist.bam)
        samtools index "${MERGED_BAM}/${file}.bam"
   
        if [ ! -s "${MERGED_BAM}/${file}.bam" ]; then
            echo "Error: Merging for ${CELL} failed. Exiting script."
            exit 1
        else
            echo "Merging for ${CELL} completed successfully."
        fi

        echo "Peak calling with MACS2 (${PVALUE_MACS2}): ${CELL}"
        macs2 callpeak \
            -t ${MERGED_BAM}/${file}.bam \
            -q ${PVALUE_MACS2} \
            -f BAMPE \
            -g ${GENOME} \
            --nomodel \
            --shift -37 \
            --extsize 73 \
            --cutoff-analysis \
            --keep-dup all \
            -n ${file}.narrow.${PVALUE_MACS2}.macs2 \
            --outdir ${PEAK_CALLING}

        # Filter out random chromosomes
        sed '/^chrUn/ d' <"${PEAK_CALLING}/${PROJECT}_${CELL}_all.narrow.${PVALUE_MACS2}.macs2_peaks.narrowPeak" >"${PEAK_CALLING}/${PROJECT}_${CELL}_all_f1.narrow.${PVALUE_MACS2}.macs2_peaks.narrowPeak"
        sed '/^chrM/ d' <"${PEAK_CALLING}/${PROJECT}_${CELL}_all_f1.narrow.${PVALUE_MACS2}.macs2_peaks.narrowPeak" >"${PEAK_CALLING}/${PROJECT}_${CELL}_all_f2.narrow.${PVALUE_MACS2}.macs2_peaks.narrowPeak"
        sed '/random/ d' <"${PEAK_CALLING}/${PROJECT}_${CELL}_all_f2.narrow.${PVALUE_MACS2}.macs2_peaks.narrowPeak" >"${PEAK_CALLING}/${PROJECT}_${CELL}_all_filtered.narrow.${PVALUE_MACS2}.macs2_peaks.narrowPeak.bed"
        
        #skipping filter by pvalue as we already provide a threshold estimated with cutoff analysis
    done
}

# Step 3. Filtering to only retain peaks that replicate in at least 2 samples
filter_peaks() {
    for CELL in "${CELLS[@]}"; do
        file="${PROJECT}_${CELL}_all"

        MERGED_PEAKS="${PEAK_CALLING}/${PROJECT}_${CELL}_all_filtered.narrow.${PVALUE_MACS2}.macs2_peaks.narrowPeak.bed"
        COMBINED_PEAKS="${PEAK_CALLING}/combined_peaks_${CELL}.bed"

        PEAK_FILES=$(ls "${PEAK_CALLING}"/*"${CELL}"*dedup.q30.no_blacklist*narrow.${PVALUE_MACS2}.macs2*.narrowPeak)
        echo "Peak files for ${CELL}:"
        echo "${PEAK_FILES}"

        for PEAK_FILE in ${PEAK_FILES}; do
            awk -v rep="${PEAK_FILE}" 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, rep}' "${PEAK_FILE}" >> "${COMBINED_PEAKS}"
        done

        #(count[key] >= 2) for at least 2 samples need to have peaks - 1 in at least 1 sample
        bedtools intersect -a "$MERGED_PEAKS" -b "$COMBINED_PEAKS" -wo | \
        awk 'BEGIN {FS="\t"; OFS="\t"} {count[$1 FS $2 FS $3 FS $4]++} END {for (key in count) if (count[key] >= 1) print key}' > "${PEAK_CALLING}/${CELL}_confidence_peaks.bed"

        rm "${COMBINED_PEAKS}"

        echo "Confidence peaks saved to ${PEAK_CALLING}/${CELL}_confidence_peaks.bed"

        awk '{print $1, "macs2", "peak", $2, $3, "0", ".", ".", "gene_id " $4}' "${PEAK_CALLING}/${CELL}_confidence_peaks.bed" > "${PEAK_CALLING}/${CELL}_confidence_peaks.gff"

        sort -k1,1 -k2,2n "${PEAK_CALLING}/${CELL}_confidence_peaks.bed" >> "${PEAK_CALLING}/${CELL}_confidence_peaks.sorted.bed"

        awk '{print $1"\t"$2-1"\t"$3"\t"$4}' FS=" " "${PEAK_CALLING}/${CELL}_confidence_peaks.sorted.bed" > "${PEAK_CALLING}/${CELL}_confidence_peaks.sorted.tmp"

        if [ -s "${PEAK_CALLING}/${CELL}_confidence_peaks.sorted.tmp" ]; then
            echo "Merged peaks file for ${CELL} is not empty."
        else
            echo "Warning: Merged peaks file for ${CELL} is empty or not found."
        fi

        sort -k1,1 -k2,2n "${PEAK_CALLING}/${CELL}_confidence_peaks.sorted.bed" >> "${PEAK_CALLING}/merged_confidence_peaks.bed"
    done

    sort -k1,1 -k2,2n "${PEAK_CALLING}/merged_confidence_peaks.bed" > "${PEAK_CALLING}/merged_confidence_peaks.sorted.bed"
    bedtools merge -i "${PEAK_CALLING}/merged_confidence_peaks.sorted.bed" > "${PEAK_CALLING}/consensus_peaks.bed"

    awk '{count++; print $1"\t"$2"\t"$3"\tNA_peak_"count}' "${PEAK_CALLING}/consensus_peaks.bed" > "${PEAK_CALLING}/consensus_peaks.saf"

    if [ -s "${PEAK_CALLING}/consensus_peaks.saf" ]; then
        echo "Final consensus peaks file is not empty."
    else
        echo "Warning: Final consensus peaks file is empty or not found."
    fi
}

# Step 4. Generate feature counts file
generate_feature_counts() {
    cd "$FEATURECOUNTS"

    awk '{OFS="\t"; print $1, ".", "peak", $2, $3, ".", "+", ".", "gene_id \""$4"\";"}' "${PEAK_CALLING}/consensus_peaks.saf" > "${PEAK_CALLING}/consensus_peaks.gtf"

    featureCounts -p -O -T 8 -a "${PEAK_CALLING}/consensus_peaks.gtf" -o "$FEATURECOUNTS/All_samples_counts.txt" -t "peak" $ANALYSIS_DIR/*.dedup.q30.no_blacklist.bam
}

# Step 5. Calculate FRIP scores
calculate_frip_scores() {
    cd $MERGED_BAM
    FRIP_OUTPUT="${ANALYSIS_DIR}/frip_scores.txt"
    echo -e "Cell\tTotal_Reads_BAM\tTotal_Reads_BED\tTotal_Reads_MACS2\tFRiP_Score_MACS2" > $FRIP_OUTPUT

    for CELL in "${CELLS[@]}"; do
        bedtools bamtobed -i ${PROJECT}_${CELL}_all.bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' > ${PROJECT}_${CELL}_all.PE2SE.tagAlign

        total_reads=$(samtools view -c ${PROJECT}_${CELL}_all.bam)
        echo "Total Reads for ${CELL} in BAM: $total_reads"

        total_reads_bed=$(wc -l < ${PROJECT}_${CELL}_all.PE2SE.tagAlign)
        echo "Total Reads for ${CELL} in BED: $total_reads_bed"

        total_reads_macs2=$(bedtools sort -i ${PEAK_CALLING}/${CELL}_confidence_peaks.sorted.bed | bedtools merge -i stdin | bedtools intersect -u -a ${PROJECT}_${CELL}_all.PE2SE.tagAlign -b stdin | wc -l)
        echo "Total Reads for ${CELL} in MACS2: $total_reads_macs2"

        frip_score_macs2=$(echo "scale=4; $total_reads_macs2 / $total_reads * 100" | bc)
        echo "FRiP Score MACS2: $frip_score_macs2%"

        echo -e "${CELL}\t${total_reads}\t${total_reads_bed}\t${total_reads_macs2}\t${frip_score_macs2}%" >> $FRIP_OUTPUT
    done

    for FILE in $(cat "${FILE_LIST_DIR}/file.txt"); do

        FINAL_FILE="${FILE%_REP1.mLb.clN.sorted.bam}.dedup.q30.no_blacklist.bam"

        bedtools bamtobed -i ${BAM_OUT_DIR}/${FILE%_REP1.mLb.clN.sorted.bam}.dedup.q30.no_blacklist.bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' > ${BAM_OUT_DIR}/${FILE%_REP1.mLb.clN.sorted.bam}.dedup.q30.no_blacklist.bam.PE2SE.tagAlign

        total_reads=$(samtools view -c ${BAM_OUT_DIR}/${FILE%_REP1.mLb.clN.sorted.bam}.dedup.q30.no_blacklist.bam)
        echo "Total Reads for ${FILE} in BAM: $total_reads"

        total_reads_bed=$(wc -l < ${BAM_OUT_DIR}/${FILE%_REP1.mLb.clN.sorted.bam}.dedup.q30.no_blacklist.bam.PE2SE.tagAlign)
        echo "Total Reads for ${FILE} in BED: $total_reads_bed"
        
        total_reads_macs2=$(bedtools sort -i $PEAK_CALLING/${FINAL_FILE}.narrow.${PVALUE_MACS2}.macs2_peaks.narrowPeak | bedtools merge -i stdin | bedtools intersect -u -a ${BAM_OUT_DIR}/${FILE%_REP1.mLb.clN.sorted.bam}.dedup.q30.no_blacklist.bam.PE2SE.tagAlign -b stdin | wc -l)
        echo "Total Reads for ${FILE} in MACS2: $total_reads_macs2"

        frip_score_macs2=$(echo "scale=4; $total_reads_macs2 / $total_reads * 100" | bc)
        echo "FRiP Score MACS2: $frip_score_macs2%"

        echo -e "${FILE}\t${total_reads}\t${total_reads_bed}\t${total_reads_macs2}\t${frip_score_macs2}%" >> $FRIP_OUTPUT
    done
}

# Step 6. Convert to bedgraph
convert_to_bedgraph() {
    for CELL in "${CELLS[@]}"; do
        file="${PROJECT}_${CELL}_all"
        echo "Converting bam into bedgraph (bamCoverage): ${CELL}"
        bamCoverage \
            --bam "${MERGED_BAM}/${file}.bam" \
            -o "${BEDGRAPHS}/${file}.bedgraph" \
            --outFileFormat "bedgraph" \
            --skipNAs \
            --binSize 50 \
            --normalizeUsing CPM \
            --extendReads \
            --effectiveGenomeSize ${EFFECTIVE_GENOME_SIZE}
    done
}

# Step 7. Sort chromosomes
sort_chromosomes() {
    for CELL in "${CELLS[@]}"; do
        file="${PROJECT}_${CELL}_all"
        echo "Sorting chromosome for bedgraph (bedtools sort): ${CELL}"
        bedtools sort -i ${BEDGRAPHS}/${file}.bedgraph > ${BEDGRAPHS}/${file}.sorted.bedgraph
    done
}

# Step 8. Clip bedgraph
clip_bedgraph() {
    for CELL in "${CELLS[@]}"; do
        file="${PROJECT}_${CELL}_all"
        echo "Clipping bedgraph for no overlap beyond chr edge (bedClip): ${CELL}"
        bedClip ${BEDGRAPHS}/${file}.sorted.bedgraph ${CHR_SIZES} ${BEDGRAPHS}/${file}.sorted.clipped.bedgraph
    done

    for CELL in "${CELLS[@]}"; do
        file="${PROJECT}_${CELL}_all"
        if [ ! -s "${BEDGRAPHS}/${file}.sorted.bedgraph" ] || [ ! -s "${BEDGRAPHS}/${file}.sorted.clipped.bedgraph" ]; then
            echo "Error: Sorting and clipping for ${CELL} failed."
            exit 1
        else
            echo "Sorting and clipping for ${CELL} completed successfully."
        fi
    done
}

# Step 9. Generate bigwig files
generate_bigwig_files() {
    mkdir -p "$BIGWIG"
    cd $MERGED_BAM

    for CELL in "${CELLS[@]}"; do
        bamCoverage --bam ${PROJECT}_${CELL}_all.bam -o $BIGWIG/${ID}_${PROJECT}_${CELL}_all.bw --binSize 50 --normalizeUsing CPM --effectiveGenomeSize ${EFFECTIVE_GENOME_SIZE} --minMappingQuality 30 
    done
}

# Main script execution
create_directories
validate_file_list
process_individual_samples
merge_bams_per_cell_type
filter_peaks
generate_feature_counts
calculate_frip_scores
convert_to_bedgraph
sort_chromosomes
clip_bedgraph
generate_bigwig_files

