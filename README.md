# Pipeline for Bulk Reads to Counts

## **Introduction**

<details>
<summary><b>Pipeline Overview</b></summary>
<br>
This pipeline consists of the following Snakemake workflow: 

<br>

**Preprocessing Workflow**: 

    Input: PE bulk RNA-seq data
    1. FASTQC 
    2. Trim (Trimmomatic)
    3. Align (bowtie2)
    4. Get Counts (FeatureCounts)

</details>

<details>
<summary><b>Getting Started</b></summary>
<br> 
1. Pull all files from repository: 

```
git clone https://github.com/emix-c/BulkReadstoCounts.git
cd BulkReadstoCounts
``` 

2. Set up the conda environment. 
```
conda env create -f counts_env.yaml
conda activate raw2counts
cd BulkReadstoCounts
```

</details>

## **Running Workflow**

<details>
<summary><b>Steps</b></summary>
<br>

1. Organize input data accordingly. 

    Input data must be paired-end **bulk** sequencing data ending with *_1.fastq.gz and *_2.fastq.gz. 
    
    **An example of input is in the example_input directory.**

<br>

2. Navigate to config/config.yaml. You will need to change the variables to match what you need. 

<br>

3. Navigate to config/cluster/config.v8+.yaml. 
    
    You can change the variables as needed. Most often, you may want to change the slurm account, where the output logs are going to and how the slurm jobs are named. 

    If desired, you can also change the default resources. 
    You are also able to change the resource requests per rule in the smk file.

<br>

4. Perform a Snakemake dry run to confirm that your data will be ran correctly. 
    ```
    cd BulkReadstoCounts
    configfile=[absolute path to config.yaml]

    snakemake -n –-profile config/cluster/ --configfile=$configfile 
    ```

    Pay close attention to the output of this dry run and check that the files Snakemake is expected to generate are correct. 

<br>

5. Now run this workflow using the launch script `scripts/run_pipeline.sh` 
    ```
    cd BulkReadstoCounts
    sbatch scripts/run_pipeline.sh -c $configfile 

    ```

</details>

<details>
<summary><b>Changing Resource Requests</b></summary>
<br>

The default resource requests can be seen by going to the cluster config file. 

I've chosen the default and some rule-specific manual adjustments through some experimenting with small samples. 

I've also manually adjusted the resource requests for some big rules in the Snakemake file. 

**NOTE: If you want to do any resource adjusting, I suggest you adjust it manually first, before trying to change the defaults.** 

To do so: 
1. Navigate to Snakemake. 
2. Look for rule [target_rule] in the code. You can see if there are any parameters for resources by looking for: 
- **threads:** (specifies # of CPUs) 
- **resources:** (specifies time, mem_mb, gres etc. )
3. Feel free to adjust your resource requests by modifying the number in threads and the values in resources. 
    
    If there are no 'threads:' or 'resources:' present under your target rule, you can add them. 


    **Example:**
    ```
    rule target: 
        input: ..
        output: ..
        params: ..
        threads: 1
        resources: 
            time="04:00:00",
            gres="disk:1024", 
            mem_mb=4000  
    ```

Want to know if your resource requests were appropriate? 

Refer to the SlurmJobAssessment tool [here](https://github.com/ohsu-cedar-comp-hub/SlurmStats)



</details>