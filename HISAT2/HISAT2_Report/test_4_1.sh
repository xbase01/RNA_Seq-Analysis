#!/bin/bash

# Set timer
SECONDS=0

# Input fastq file for testing
TEST_FASTQ="SRR3414637_1.fastq"

# Extract the sample ID from the filename
sample_id=$(echo "$TEST_FASTQ" | cut -d '_' -f 1)

# Run HISAT2 alignment for the test sample
hisat2 -p 4 -x /home/Students/M06-310/JudeAneke/HISAT2_indexes/data/HISAT2/HISAT2_indexes/GRCm39_idx -U "/home/Students/M06-310/data/$TEST_FASTQ" -S "${sample_id}.sam" 2> "${sample_id}.hisat"
echo "HISAT2 finished running for ${sample_id}!"

# Filter the alignment file to include only reads with NH:i:1
grep -P '^@|NH:i:1$' "${sample_id}.sam" > "${sample_id}.uniq.sam"
echo "Uniquely mapped reads created successfully for ${sample_id}!"

# Run samtools to create sorted bam file
samtools view -b "${sample_id}.uniq.sam" | samtools sort -o "${sample_id}.uniq.sorted.bam"
echo "Bam file created and sorted successfully for ${sample_id}!"

# Create Index bam file
samtools index "${sample_id}.uniq.sorted.bam"
echo "Bam indexed successfully for ${sample_id}!"

# Run htseq-count for the test sample
htseq-count --stranded=no "${sample_id}.uniq.sorted.bam" /home/Students/M06-310/JudeAneke/HISAT2_indexes/data/gencode.vM34.primary_assembly.annotation.gtf > "${sample_id}.counts"
echo "Count file created successfully for ${sample_id}!"

# Calculate and display the elapsed time
duration=$SECONDS
echo "Alignment, filtering, and counting completed for ${sample_id} in $(($duration / 60)) minutes and $(($duration % 60)) seconds."

# Exit the script
exit

