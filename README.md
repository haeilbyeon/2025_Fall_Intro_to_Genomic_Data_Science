# Metagenomic Preprocessing and Microbial Community Analyses

# [Linux] Preprocessing metagenomic datasets

**Step 1. Trimgalore**

NOTE: The purpose of TrimGalore is to trim sequences, remove short sequences, and remove sequences with low quality scores.

Before running the analysis, first create a directory for all TrimGalore output files.

    mkdir -p /insert/your/filepath/here/01_Trim_Galore/ # change to your preferred directory path

TrimGalore was run using a bash script: 01_trim_galore_job.slurm

**Script values to change for your specific setup:**
 - Input directory path: where all the raw data files are located
 - Output directory path: where all the sequence files coming out of TrimGalore should be saved
 - Minimum quality score value (To replicate the paper, we used 30)
 - Minimum sequence length (To replicate the paper, we used 60)

**Inputs:** Raw paired FASTQ/FASTQ.GZ files for each sample
 - Example: XX_1.fastq and XX_2.fastq

**Outputs:** Trimmed FASTQ files, .html and .zip FastQC output files, and a .txt trimming report file for raw sequence file
 - Example: XX_1_val_1.fq, XX_1_val_1_fastqc.html, XX_1_val_1_fastqc.zip, and XX_1.fastq_trimming_report.txt

**What the script does:**

 1. Loads the TrimGalore module on ARC
 2. Defines the input directory where the raw data files are located and the output directory where the trimmed sequence files will be saved
 3. Changes the working directory to the input directory and loops through all sequence files to run them through TrimGalore

Once all batch script parameters and filepaths are updated in the script, navigate to the directory that your script is located in and use the following line of code to submit the job on ARC.

    sbatch 01_trim_galore_job.slurm

**Script: 01_trim_galore_job.slurm**

    #!/bin/bash
    #SBATCH --job-name=trim_galore
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

**Step 2. Superdeduper**

NOTE: The purpose of SuperDeduper is to remove PCR and optical duplicate reads from paired-end FASTQ files generated after adapter trimming (Trim Galore Step). This step ensures that downstream analyses such as assembly or host mapping are not biased by artificially duplicated reads.

Before running the analysis, first create a directory for all SuperDeduper output files:

    mkdir -p /insert/your/filepath/here/02_after_deduper/

**02a. Preparation of the environment**
To run superdeduper, we first need to build a conda virtual environment that contains the HTStream toolkit, which can provide the hts_SuperDeduper program.

First, we need to load Miniconda

    ## Load Miniconda module
    module load Miniconda3
    
    ## check miniconda version
    conda --version
    
Secondly, we create a virtual environment called htstream12

    ＃ Create and activate the environment
    conda create -y -n htstream12 python=3.10
    conda activate htstream12

Thirdly, Install HTStream (which includes SuperDeduper). HTStream can be installed via Bioconda, which provides precompiled packages for bioinformatics tools.

    # 1. Add the proper Conda channels
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    
    # 2. Install HTStream 
    conda install -y htstream
    
    # 3. Verify the correct installation of SuperDeduper
    # verify installation 
    which hts_SuperDeduper
    # Check version
    hts_SuperDeduper --version

**Run SuperDeduper Analysis**

Once the conda environment was set up, SuperDeduper analysis was run using a bash script: 02_run_superdeduper_cpu.slurm

**Script values to change for your specific setup:**

 - Input directory path: Where all the output files from TrimGalore are located
 - Output directory path: Where all the sequence files coming out of SuperDeduper should be saved

**Inputs:** Paired FASTQ/FASTQ.GZ files from TrimGalore output

 - Example: XX_1_val_1.fq and XX_1_val_.fq

**Outputs:** Deduplicated paired FASTQ files

 - Example: XX_R1.fastq.gz and XX_R2.fastq.gz

**What the script does:**

 1. Loads the conda environment (htstream12) that contains the HTStream toolkit, including SuperDeduper.
 2. Defines both the input and output directory paths.
 3. Creates the output directory if it does not already exist.
 4. Scans through all *_1_val_1.fq files in the input folder with a for loop.
 5. For each detected sample, identifies its mate file (*_2_val_2.fq) and runs SuperDeduper
 6. Checks whether both deduplicated files (_R1.fastq.gz and _R2.fastq.gz) were successfully created.
 7. Prints progress messages such as [INFO], [OK], or [WARN] to help track pipeline execution.

Once all batch script parameters and filepaths are updated in the script, navigate to the directory that your script is located in and use the following line of code to submit the job on ARC.

    sbatch 02_run_superdeduper_cpu.slurm

**Script: 02_run_superdeduper_cpu.slurm**

    #!/bin/bash
    #SBATCH --job-name=superdeduper       # Job name shown in squeue and logs
    #SBATCH --account=leaph          # Project/account to charge compute to
    #SBATCH --partition=normal_q          # Queue/partition to run on
    #SBATCH --time=24:00:00               # Walltime limit (HH:MM:SS)
    #SBATCH --cpus-per-task=8             # CPU threads requested for this job
    #SBATCH --mem=32G                     # Total memory requested
    #SBATCH --output=/projects/leaph/Haeil/02_SuperDeduper/logs/%x_%j.out  # Stdout log file
    #SBATCH --error=/projects/leaph/Haeil/02_SuperDeduper/logs/%x_%j.err   # Stderr log file
    #SBATCH --mail-type=ALL
    #SBATCH --mail-user=hibyeon@vt.edu
    
    shopt -s nullglob                     # If no files match a glob, expand to empty (not the pattern)
    
    # ---- Paths ----
    INDIR="/projects/leaph/Haeil/01_Trim_Galore"
    OUTDIR="/projects/leaph/Haeil/02_SuperDeduper"
    
    # ---- Activate software environment ----
    source /projects/intro2gds/I2GDS2025/tools/miniconda3/etc/profile.d/conda.sh  # Load conda functions
    conda activate /projects/leaph/Haeil/envs/htstream12            # Activate env that provides hts_SuperDeduper
    
    # ---- Collect all R1 files (Trim Galore naming: *_R1_val_1.fq.gz) ----
    R1S=( "$INDIR"/*_R1_val_1.fq.gz )        # Array of R1 fastq files to iterate over
    [[ ${#R1S[@]} -eq 0 ]] && echo "[INFO] No *_R1_val_1.fq.gz under $INDIR; nothing to do."  # Inform when empty set
    
    # ---- Main loop over samples ----
    for R1 in "${R1S[@]}"; do            # Iterate through every R1 file found
      base=$(basename "$R1")             # Get filename only (strip directory)
      sample=${base%_R1_val_1.fq.gz}         # Strip suffix to get sample ID (e.g., SRRxxxx)
    
      R2="$INDIR/${sample}_R2_val_2.fq.gz"   # Expected R2 mate based on Trim Galore naming
      if [[ ! -f "$R2" ]]; then          # If paired R2 is missing,
        echo "[WARN] Missing R2 for $sample; skip."  # warn and skip this sample
        continue
      fi
    
      prefix="$OUTDIR/${sample}"         # Output prefix used by HTStream tools
      echo "[INFO] $sample -> hts_SuperDeduper"   # Progress message
    
      # Run deduplication; -1/-2: paired inputs; -f: output prefix; -F: overwrite if exists
      hts_SuperDeduper -1 "$R1" -2 "$R2" -f "$prefix" -F || { echo "[WARN] hts_SuperDeduper failed for $sample; continue."; continue; }
    
      # Quick sanity check that expected gzipped outputs exist and are non-empty
      if [[ -s "${prefix}_R1.fastq.gz" && -s "${prefix}_R2.fastq.gz" ]]; then
        echo "[OK] $sample done."
      else
        echo "[WARN] Outputs missing/empty for $sample."
      fi
    done
    
    echo "[DONE] Pipeline finished."      # Final message when loop completes

# [Python] Downstream analyses of microbial communities

Step 3. Alpha diversity plots



Step 4. Beta diversity plots



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
