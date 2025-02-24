#!/bin/bash

WORKDIR=/mnt/x/RNAseq

# Samples beginning with A
(
  for sample_dir in "$WORKDIR"/A*; do
    if [ -d "$sample_dir" ]; then
      if [ ! -d "$sample_dir/trimmed" ]; then
        if ls "$sample_dir"/*.fastq.gz &> /dev/null; then
          echo "Processing sample: $sample_dir"
          R1=$(ls "$sample_dir"/*_R1*.fastq.gz 2> /dev/null)
          R2=$(ls "$sample_dir"/*_R2*.fastq.gz 2> /dev/null)
          if [ -f "$R1" ] && [ -f "$R2" ]; then
            echo "Running trim_galore on $R1 and $R2"
            trim_galore --paired \
                        --cores 6 \
                        --quality 20 \
                        --phred33 \
                        --stringency 3 \
                        --length 20 \
                        --output_dir "$sample_dir/trimmed" \
                        "$R1" "$R2"
          fi
        fi
      else
        echo "Skipping sample: $sample_dir (already trimmed)"
      fi
    fi
  done
) &

# Samples beginning with B
(
  for sample_dir in "$WORKDIR"/B*; do
    if [ -d "$sample_dir" ]; then
      if [ ! -d "$sample_dir/trimmed" ]; then
        if ls "$sample_dir"/*.fastq.gz &> /dev/null; then
          echo "Processing sample: $sample_dir"
          R1=$(ls "$sample_dir"/*_R1*.fastq.gz 2> /dev/null)
          R2=$(ls "$sample_dir"/*_R2*.fastq.gz 2> /dev/null)
          if [ -f "$R1" ] && [ -f "$R2" ]; then
            echo "Running trim_galore on $R1 and $R2"
            trim_galore --paired \
                        --cores 6 \
                        --quality 20 \
                        --phred33 \
                        --stringency 3 \
                        --length 20 \
                        --output_dir "$sample_dir/trimmed" \
                        "$R1" "$R2"
          fi
        fi
      else
        echo "Skipping sample: $sample_dir (already trimmed)"
      fi
    fi
  done
) &

wait
echo "All parallel runs finished."