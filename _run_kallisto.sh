#!/bin/bash

# Define the path to the RNAseq folder
RNAseq_folder="/mnt/x/analysis/RNAseq"

# Define the path to the kallisto index
kallisto_index="/mnt/x/analysis/ref/CriGri-PICRH-1.idx"

# Loop through all trimmed sample folders and run kallisto
WORKDIR=/mnt/x/RNAseq

for sample_dir in "$WORKDIR"/*; do
    if [ -d "$sample_dir" ]; then 
        sample_name=$(basename "$sample_dir")
        echo "Processing $sample_name $sample_dir"
        
        # Define the output directory
        output_dir="/mnt/x/analysis/kallisto_results/$sample_name"
        abundance_file="$output_dir/abundance.tsv"
        
        # Check if abundance.tsv already exists
        if [ -f "$abundance_file" ]; then
            echo "Skipping $sample_name as $abundance_file already exists."
            continue
        else
            echo "Running kallisto quant for $sample_name."
            mkdir -p "$output_dir"
            
            # Run kallisto quant
            # kallisto quant -t 14 -i $kallisto_index -o $output_dir /mnt/x/RNAseq/${sample_name}/trimmed/merged_R1_val_1.fq.gz /mnt/x/RNAseq/${sample_name}/trimmed/merged_R2_val_2.fq.gz
            kallisto quant -t 14 -i $kallisto_index -o $output_dir /mnt/x/RNAseq/${sample_name}/merged_R1.fastq.gz /mnt/x/RNAseq/${sample_name}/merged_R2.fastq.gz
        fi
        
    fi
done
 