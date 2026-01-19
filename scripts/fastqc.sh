#!/bin/bash

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --r1=*)
            r1="${1#*=}"
            ;;
        --r2=*)
            r2="${1#*=}"
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


if [[ -z "$r1" || -z "$r2" || -z "$outdir" ]]; then
    echo "Error: One or more required arguments are missing."
    exit 1
fi


# Check if output directory exists, if not, create it
if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi


fastqc $r1 -o $outdir 
fastqc $r2 -o $outdir 





