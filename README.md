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

**Step 3. Alpha diversity plots**

NOTE: Alpha diversity analysis measures the diversity within a single sample, describing how many taxa are present (richness) and how evenly they are distributed (evenness).


**0) Import libraries**
 - pandas / numpy: handle tables and simple math.
 - seaborn / matplotlib: make plots.
 - scipy.stats.kruskal: runs the Kruskal–Wallis test (overall difference between groups).
 - statannotations: adds p-value “stars/brackets” on the plot.
 - scikit-bio alpha: calculates alpha diversity indices (Shannon, Simpson, Chao1, ACE).


**1) Tell Python where the data files are (by taxonomic level)**

level_files is a dictionary like:
 - Phylum → path to Phylum count table
 - Class → path to Class count table …and so on.

So the script can repeat the same analysis for each level automatically.


**2) Define which sample belongs to which site**

site_table is a small table that maps:
Sample ID (e.g., 21-2672)
→ Site (e.g., RCR01)

site_order = ["RCR01", "BRB01", "RCR0P"] forces plots to always show sites in this order.


**3) Choose which alpha diversity metrics to compute**

alpha_metrics = ["Shannon", "Simpson", "Chao1", "ACE"]

These are different ways to summarize “how diverse” each sample is.


**4) Compute alpha diversity for each sample**

Function: compute_alpha_indices(count_df)

What it does:
 1. Assumes the first column is the taxon name (like species/genus names).
 2. Converts the remaining columns into numeric counts (fills missing with 0, makes integers).
 3. For each sample column, it calculates:
  - Shannon (diversity with evenness)
  - Simpson (dominance/evenness-related)
  - Chao1 (richness estimator)
  - ACE (another richness estimator)
 4. Returns a new table like:
| Sample | Shannon | Simpson | Chao1 | ACE |


**5) Plot one metric: boxplot + dots + statistics**

Function: plot_alpha_panel(ax, data, level_name, metric)

For one metric (ex: Shannon):
 1. Kruskal–Wallis test across the 3 sites
  - Gives one p-value: “Is any site different overall?”
 2. Draws:
  - Boxplot per site (distribution)
  - Swarmplot dots (each individual sample point)
 3. Runs pairwise Mann–Whitney tests between:
  - RCR01 vs BRB01
  - RCR01 vs RCR0P
  - BRB01 vs RCR0P

Then it annotates the plot with brackets/p-values.

 4. Expands the y-axis a bit so annotations fit.
 5. Writes the Kruskal–Wallis p-value at the top.
 6. Draws red dashed lines showing the mean value for each site.


**6) Loop through every taxonomic level and make a 2×2 panel plot**

For each level (Phylum, Class, Order, Family, Genus, Species):
 1. Read the CSV (pd.read_csv(path)).
 2. Compute alpha metrics (compute_alpha_indices).
 3. Merge with site info (pd.merge(..., site_table ...)) so each sample has a Site label.
 4. Create a multi-panel figure (2 columns × enough rows).
 5. Plot all 4 metrics (Shannon/Simpson/Chao1/ACE), one per panel.
 6. Show the figure.

**Script**

    import numpy as np
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy.stats import kruskal
    from statannotations.Annotator import Annotator
    from skbio.diversity import alpha as skbio_alpha
    
    # ============================================
    # 0. File paths for each taxonomic level
    # ============================================
    
    level_files = {
        "Phylum":  r"C:\Users\dygks\Downloads\Research_MST\Original_data\MG_Bracken_using_DB_LEAPH_output\Phylum\bracken_Phylum_new_est_reads_matrix.csv",
        "Class":   r"C:\Users\dygks\Downloads\Research_MST\Original_data\MG_Bracken_using_DB_LEAPH_output\Class\bracken_Class_new_est_reads_matrix.csv",
        "Order":   r"C:\Users\dygks\Downloads\Research_MST\Original_data\MG_Bracken_using_DB_LEAPH_output\Order\bracken_Order_new_est_reads_matrix.csv",
        "Family":  r"C:\Users\dygks\Downloads\Research_MST\Original_data\MG_Bracken_using_DB_LEAPH_output\Family\bracken_Family_new_est_reads_matrix.csv",
        "Genus":   r"C:\Users\dygks\Downloads\Research_MST\Original_data\MG_Bracken_using_DB_LEAPH_output\Genus\bracken_Genus_new_est_reads_matrix.csv",
        "Species": r"C:\Users\dygks\Downloads\Research_MST\Original_data\MG_Bracken_using_DB_LEAPH_output\Species\bracken_Species_new_est_reads_matrix.csv",
    }
    
    # ============================================
    # 1. Sample ↔ Site mapping and site order
    # ============================================
    
    site_order = ["RCR01", "BRB01", "RCR0P"]
    
    site_table = pd.DataFrame({
        "Sample": [
            "21-2672", "21-2688", "21-3862",
            "21-2673", "21-2689", "21-3863",
            "21-2674", "21-2690", "21-3864",
        ],
        "Site": [
            "RCR01", "RCR01", "RCR01",
            "BRB01", "BRB01", "BRB01",
            "RCR0P", "RCR0P", "RCR0P",
        ],
    })
    site_table["Site"] = pd.Categorical(
        site_table["Site"],
        categories=site_order,
        ordered=True
    )
    
    # Alpha metrics to calculate
    alpha_metrics = ["Shannon", "Simpson", "Chao1", "ACE"]
    
    
    # ============================================
    # 2. Compute alpha indices
    # ============================================
    
    def compute_alpha_indices(count_df: pd.DataFrame) -> pd.DataFrame:
        """
        Compute Shannon, Simpson, Chao1, ACE for each sample using scikit-bio.
        """
        taxon_col = count_df.columns[0]
    
        # Clean & convert counts to integers (required for ACE)
        counts_only = (
            count_df
            .set_index(taxon_col)
            .apply(pd.to_numeric, errors="coerce")
            .fillna(0)
            .astype(int)
        )
    
        records = []
        for sample in counts_only.columns:
            counts = counts_only[sample].values
    
            shannon = skbio_alpha.shannon(counts)
            simpson = skbio_alpha.simpson(counts)
            chao1 = skbio_alpha.chao1(counts)
            ace = skbio_alpha.ace(counts)
    
            records.append({
                "Sample": sample,
                "Shannon": shannon,
                "Simpson": simpson,
                "Chao1": chao1,
                "ACE": ace
            })
    
        return pd.DataFrame(records)
    
    
    # ============================================
    # 3. Plotting function (ticks visible, no rotation)
    # ============================================
    
    def plot_alpha_panel(ax, data: pd.DataFrame, level_name: str, metric: str):
        """
        Creates a boxplot + swarmplot for a single alpha metric,
        includes KW test, pairwise Mann-Whitney, annotation,
        and visible ticks (no rotation on x-axis labels).
        """
        data = data.copy()
        data["Site"] = pd.Categorical(data["Site"], categories=site_order, ordered=True)
    
        # ---------- Kruskal–Wallis ----------
        grouped = [group[metric].values for _, group in data.groupby("Site")]
        kw_stat, kw_p = kruskal(*grouped)
    
        # ---------- Pairwise comparisons ----------
        pairs = [("RCR01", "BRB01"),
                 ("RCR01", "RCR0P"),
                 ("BRB01", "RCR0P")]
    
        # ---------- Plot ----------
        sns.boxplot(x="Site", y=metric, data=data, ax=ax,
                    palette="Set2", showfliers=False)
        sns.swarmplot(x="Site", y=metric, data=data,
                      ax=ax, color="grey", s=6, alpha=0.8, edgecolor="black")
    
        # ---------- Stat annotation ----------
        annot = Annotator(ax, pairs, data=data, x="Site", y=metric)
        annot.configure(test='Mann-Whitney', comparisons_correction=None, verbose=0)
        annot.apply_test()
        annot.annotate()
    
        # ---------- Extend Y limits ----------
        y_min, y_max = ax.get_ylim()
        ax.set_ylim(y_min - 0.10 * (y_max - y_min), y_max + 0.10 * (y_max - y_min))
    
        # ---------- KW P-value ----------
        ax.text(0.5, 0.94, f"KW P = {kw_p:.3e}",
                fontsize=12, ha="center", transform=ax.transAxes)
    
        # ---------- Group means (red dashed line) ----------
        site_means = data.groupby("Site")[metric].mean()
        for i, site in enumerate(site_order):
            ax.hlines(y=site_means[site], xmin=i - 0.4, xmax=i + 0.4,
                      colors="red", linestyles="--", linewidth=1.2)
    
        # ---------- Force ticks visible ----------
        ax.tick_params(axis="x", which="both", bottom=True, labelbottom=True)
        ax.tick_params(axis="y", which="both", left=True, labelleft=True)
    
        # ---------- Labels ----------
        ax.set_title(f"{level_name} - {metric} by Site", fontsize=13)
        ax.set_xlabel("Site")
        ax.set_ylabel(metric)
    
    
    # ============================================
    # 4. Loop over levels: compute + plot
    # ============================================
    
    for level, path in level_files.items():
        print("\n==============================================")
        print(f"[INFO] Processing level: {level}")
        print("==============================================")
    
        df_counts = pd.read_csv(path)
    
        # Compute alpha diversity
        alpha_df = compute_alpha_indices(df_counts)
    
        # Attach site metadata
        alpha_site = pd.merge(alpha_df, site_table, on="Sample", how="inner")
    
        print("[INFO] Alpha diversity with site mapping:")
        print(alpha_site)
    
        # Create multi-panel figure (2 × 2)
        n_metrics = len(alpha_metrics)
        ncols = 2
        nrows = int(np.ceil(n_metrics / ncols))
    
        fig, axes = plt.subplots(nrows, ncols, figsize=(3.5 * ncols, 4.2 * nrows))
        axes = np.array(axes).reshape(-1)
    
        # Plot each metric
        for i, metric in enumerate(alpha_metrics):
            plot_alpha_panel(axes[i], alpha_site, level, metric)
    
        # Remove unused axes
        for j in range(n_metrics, len(axes)):
            fig.delaxes(axes[j])
    
        plt.tight_layout()
        plt.show()

![Image Alt](https://github.com/haeilbyeon/2025_Fall_Intro_to_Genomic_Data_Science/blob/989cb43d1c7a3c66c2b303831118c728475406b3/Alpha%20diversity%20plot.png)

Step 4. Beta diversity plots

![Image Alt](https://github.com/haeilbyeon/2025_Fall_Intro_to_Genomic_Data_Science/blob/989cb43d1c7a3c66c2b303831118c728475406b3/Beta%20diversity%20plot.png)













