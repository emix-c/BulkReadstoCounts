#!/bin/bash
export _JAVA_OPTIONS="-Xmx4G"

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --r1=*)
            r1="${1#*=}"
            ;;
        --r2=*)
            r2="${1#*=}"
            ;;
        --s_info=*)
            s_info="${1#*=}"
            ;;
        --adapter=*)
            adapter="${1#*=}"
            ;;
        --threads=*)
            t="${1#*=}"
            ;;
        --trailing=*)
            trailing="${1#*=}"
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


if [[ -z "$r1" || -z "$r2" || -z "$s_info" || -z "$adapter" || -z "$t"  || -z "$trailing" || -z "$outdir" ]]; then
    echo "Error: One or more required arguments are missing."
    exit 1
fi


# Check if output directory exists, if not, create it
if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

trimmed_1=${outdir}/${s_info}_1_trimmed.fastq.gz
trimmed_2=${outdir}/${s_info}_2_trimmed.fastq.gz

untrimmed_1=${outdir}/${s_info}_1_untrimmed.fastq.gz
untrimmed_2=${outdir}/${s_info}_2_untrimmed.fastq.gz


trimmomatic PE -threads $t -phred33 $r1 $r2 \
                $trimmed_1 $untrimmed_1 \
                $trimmed_2 $untrimmed_2 \
                ILLUMINACLIP:$adapter  TRAILING:$trailing

fastqc $trimmed_1 -o $outdir
fastqc $trimmed_2 -o $outdir
