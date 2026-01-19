import os 
import logging 
import re 
import pandas as pd 
import glob

main_dir = config['raw_data']
results_dir = config['results']

# raw data files should have both R1 and R2 
R1_NAMES = glob.glob(os.path.join(main_dir, "*_1.fastq.gz"))

# Extract sample names (without _1.fastq.gz)
SAMPLE_NAMES = [os.path.basename(f).replace("_1.fastq.gz", "") for f in R1_NAMES]

# Keep only samples that have both _1 and _2
NAMES = [
    name for name in SAMPLE_NAMES
    if os.path.exists(os.path.join(main_dir, f"{name}_2.fastq.gz"))
]

print("Sample names are:", NAMES)


rule all: 
    fastqc_1 = os.path.join(main_dir, "FASTQC/{sample}_1_fastqc.html"),
    fastqc_2 = os.path.join(main_dir, "FASTQC/{sample}_2_fastqc.html"), 
    trimmed_1 = os.path.join(main_dir, "Trimmed-FASTQC/{sample}_1_trimmed.fastq.gz"), 
    trimmed_2 = os.path.join(main_dir, "Trimmed-FASTQC/{sample}_2_trimmed.fastq.gz") 

rule fastqc: 
    input: 
        raw_1 = os.path.join(main_dir, "{sample}_1.fastq.gz"), 
        raw_2 = os.path.join(main_dir, "{sample}_2.fastq.gz")
    params: 
        fastqc_dir = os.path.join(main_dir, "FASTQC")
    output: 
        fastqc_1 = os.path.join(main_dir, "FASTQC/{sample}_1_fastqc.html"),
        fastqc_2 = os.path.join(main_dir, "FASTQC/{sample}_2_fastqc.html")
    shell: 
        """
        {config[fastqc]} /
            --r1={input.raw_1}
            --r2={input.raw_2}
            --outdir={params.fastqc_dir}

        """

rule trim: 
    input: 
        raw_1 = os.path.join(main_dir, "{sample}_1.fastq.gz"), 
        raw_2 = os.path.join(main_dir, "{sample}_2.fastq.gz")
    params: 
        trim_fastqc_dir = os.path.join(main_dir, "Trimmed-FASTQC")
    output: 
        trimmed_1 = os.path.join(main_dir, "Trimmed-FASTQC/{sample}_1_trimmed.fastq.gz"), 
        trimmed_2 = os.path.join(main_dir, "Trimmed-FASTQC/{sample}_2_trimmed.fastq.gz"), 
        trim_fastqc_1 = os.path.join(main_dir, "Trimmed-FASTQC/{sample}_1_trimmed_fastqc.html"),
        trim_fastqc_2 = os.path.join(main_dir, "Trimmed-FASTQC/{sample}_2_trimmed_fastqc.html")
    shell: 
        """
        {config[trim]} /
            --r1={input.raw_1}
            --r2={input.raw_2}
            --s_info={wildcards.sample}
            --adapter={config[adapter]}
            --threads={config[trim_threads]}
            --trailing={config[trailing_nts]}
            --outdir={params.trim_fastqc_dir}
        """

rule align: 
    input: 
        trimmed_1 = os.path.join(main_dir, "Trimmed-FASTQC/{sample}_1_trimmed.fastq.gz"), 
        trimmed_2 = os.path.join(main_dir, "Trimmed-FASTQC/{sample}_2_trimmed.fastq.gz")
    params: 
        align_dir = os.path.join(main_dir,"Aligned")
    output: 
        



        




