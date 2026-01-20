#!/bin/bash
#SBATCH --partition=batch 
#SBATCH --account=cedar
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb
#SBATCH --time=3:00:00
#SBATCH --job-name=get_counts

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --bam=*)
            bam="${1#*=}"
            ;;
        --gtf=*)
            gtf="${1#*=}"
            ;;
        --s_info=*)
            s_info="${1#*=}"
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


if [[ -z "$bam" || -z "$gtf" || -z "$s_info" || -z "$outdir" ]]; then
    echo "Error: One or more required arguments are missing."
    exit 1
fi


# Check if output directory exists, if not, create it
if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

out_counts=${outdir}/${s_info}_feature_counts.txt

featureCounts -p -a $gtf -o $out_counts $bam



