# Metagenomic Preprocessing and Microbial Community Analyses

# [Linux] Preprocessing metagenomic datasets.
Step 1. Trimgalore
The purpose of TrimGalore is to trim sequences, remove short sequences, and remove sequences with low quality scores.
Before running the analysis, first create a directory for all TrimGalore output files.

    mkdir -p /insert/your/filepath/here/01_Trim_Galore/ # change to your preferred directory path



    #!/bin/bash
    #SBATCH --job-name=trim_galore_G5
    #SBATCH --account=leaph
    #SBATCH --partition=normal_q
    #SBATCH --nodes=1
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=8
    #SBATCH --mem=32G
    #SBATCH --time=24:00:00
    #SBATCH --output=trim_galore.%j.out
    #SBATCH --error=trim_galore.%j.err
    #SBATCH --mail-type=ALL
    #SBATCH --mail-user=hibyeon@vt.edu

    # ============================================================
    # 1. Load Module
    module load Trim_Galore/0.6.10-GCCcore-12.3.0
    
    # ------------------------------
    # 2. Define directories
    RAW_DIR="/projects/leaph/shared/project_data/RockCreekMST/CosmosMicrobiomeData"
    OUT_DIR="/projects/leaph/Haeil/01_Trim_Galore"
    
    # ------------------------------
    # 3. Run Trim Galore for each sample
    cd "$RAW_DIR"

    for R1 in *_R1.fastq.gz
    
    do
        SAMPLE=${R1%_R1.fastq.gz}
        R2=${SAMPLE}_R2.fastq.gz
        echo "Processing sample: $SAMPLE"
        trim_galore --paired --cores 8 --quality 30 --length 60 --fastqc --fastqc_args "--outdir $OUT_DIR" --output_dir "$OUT_DIR" "$RAW_DIR/$R1" "$RAW_DIR/$R2"
    done
    
    # ------------------------------
    # 4. Done
    echo "All samples processed successfully."

Step 2. Superdeduper

# [Python] Downstream analyses of microbial communities
Step 3. Relative abundance bar plots

Step 4. Alpha diversity plots

Step 5. Beta diversity plots



# ----------------------------------------------------------------------------------------------
(1) Useful to yourself in the long run​
 - Include at least 2 notes to yourself about potential mistakes (that you had made during the class), that you don't want to make again.​
(2) Clearly written and easy to follow by a peer​
 - To your future self who is likely to forget what you have done.​
 - Like a "lab notebook", will be the basis for writing a "method section" in your future publications​
(3) Contain both Linux part and (R or Python) visualization ​​
(4) Not the entire pipeline​
 - One or two steps of Linux​
 - One or two figures in R/Python visualization ​
(5) Describe your output and how the visualization supports/illustrates the point (like a "result" section in a paper)
(6) Submit the link to your own GitHub repository
