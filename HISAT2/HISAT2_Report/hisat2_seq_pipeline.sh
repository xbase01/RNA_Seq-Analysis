#!/bin/bash

# Set timer
SECONDS=0

# pwd = HISAT2

# Build HISAT2 index
hisat2-build /home/Students/M06-310/JudeAneke/HISAT2_indexes/data/GRCm39.genome.fa HISAT2_indexes/GRCm39_idx

# List of input fastq files
FASTQ_FILES=("SRR3414629_1.fastq" "SRR3414630_1.fastq" "SRR3414631_1.fastq" "SRR3414635_1.fastq" "SRR3414636_1.fastq" "SRR3414637_1.fastq")

# Path to Bowtie2 index
HISAT2_INDEX="/home/Students/M06-310/JudeAneke/HISAT2_indexes/data/HISAT2/HISAT2_indexes/GRCm39_idx"

# Loop through each fastq file
for file in "${FASTQ_FILES[@]}"; do
    # Extract the sample ID from the filename
    sample_id=$(echo "$file" | cut -d '_' -f 1)

    # Run HISAT2 alignment for each sample
    hisat2 -p 4 -x "$HISAT2_INDEX" -U "/home/Students/M06-310/data/$file" -S "${sample_id}.sam" 2> "${sample_id}.hisat"
    echo "HISAT2 finished running for ${sample_id}!"

    # Filter the alignment file to include only reads with NH:i:1
    grep -P '^@|NH:i:1$' "${sample_id}.sam" > "${sample_id}.uniq.sam"
    echo "Uniquely mapped reads created successfully for ${sample_id}!"

    # Run samtools to create sorted bam files
    samtools view -b "${sample_id}.uniq.sam" | samtools sort -o "${sample_id}.uniq.sorted.bam"
    echo "Bam file created and sorted successfully for ${sample_id}!"

    # Create Index bam file
    samtools index "${sample_id}.uniq.sorted.bam"
    echo "Bam indexed successfully for ${sample_id}!"

    # Run htseq-count for each sample
    htseq-count --stranded=no "${sample_id}.uniq.sorted.bam" /home/Students/M06-310/JudeAneke/HISAT2_indexes/data/gencode.vM34.primary_assembly.annotation.gtf > "${sample_id}.counts"
    echo "Count files created successfully for ${sample_id}!"

done

# Calculate and display the elapsed time
duration=$SECONDS
time_summary="Alignment, filtering, and counting completed in $(($duration / 60)) minutes and $(($duration % 60)) seconds."

# Echo the time summary to a file
echo "$time_summary" > hisat2_time_summary.txt

# Display the time summary
echo "$time_summary"

# Exit the script
exit



