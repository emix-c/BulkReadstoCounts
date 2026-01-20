#!/bin/bash

#SBATCH --partition=batch 
#SBATCH --account=cedar
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=10gb
#SBATCH --time=5:00:00
#SBATCH --job-name=align_test

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --trimmed_1=*)
            trimmed_1="${1#*=}"
            ;;
        --trimmed_2=*)
            trimmed_2="${1#*=}"
            ;;
        --s_info=*)
            s_info="${1#*=}"
            ;;
        --genome=*)
            genome="${1#*=}"
            ;;
        --threads=*)
            t="${1#*=}"
            ;;
        --outdir=*)
            outdir="${1#*=}"
            ;;
        *)
            echo "Unknown parameter: $1"
            exit 1
            ;;
    esac
    shift
done


if [[ -z "$trimmed_1" || -z "$trimmed_2" || -z "$s_info" || -z "$genome" || -z "$t" || -z "$outdir" ]]; then
    echo "Error: One or more required arguments are missing."
    exit 1
fi


# Check if output directory exists, if not, create it
if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

output_file=${outdir}/${s_info}_aligned.bam


bowtie2 -x $genome \
    -1 $trimmed_1 \
    -2 $trimmed_2 \
    -p $t \
    -S test.sam 

samtools view -bS temp.sam | samtools sort -o $output_file 

rm temp.sam



