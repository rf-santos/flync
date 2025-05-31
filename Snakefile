import yaml

# Load configuration
with open("test/config.yaml", "r") as f:
    config = yaml.safe_load(f)

# Define application directory and working directory from config or use defaults
appdir = config.get("appdir", "workflow/scripts")
workdir = config.get("workdir", "results")

# Get SRA accessions from config
# Assuming sra_list_file in config points to a file with one accession per line
with open(config["sra_list_file"], "r") as f:
    SRA_ACCESSIONS = [line.strip() for line in f if line.strip()]

# Get sample names from config for local FASTQ files if provided
# This part needs to be more robust based on how samples are defined.
# For now, let's assume a simple list of sample names if not using SRA.
# SAMPLES = config.get("samples", []) # This will be replaced by keys from config["samples_config"]

# Define what samples to run based on the samples_config in config.yaml
# This allows mixed SRA and local FASTQ samples
ALL_SAMPLES = list(config.get("samples_config", {}).keys())

# Parse tracksFile.tsv for remote tracks
def get_tracks_from_file(file_path):
    tracks = {}
    try:
        with open(file_path) as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    try:
                        name, path = line.strip().split('\t', 1) # Split only on the first tab
                        tracks[name] = path
                    except ValueError:
                        print(f"Warning: Skipping malformed line in {file_path}: {line.strip()}")
    except FileNotFoundError:
        print(f"Warning: Tracks file {file_path} not found. Remote tracks will be empty.")
    return tracks

TRACKS_REMOTE_FILE = config.get("tracks_file", "static/tracksFile.tsv")
TRACKS_REMOTE = get_tracks_from_file(TRACKS_REMOTE_FILE)
ALL_TRACK_NAMES = list(TRACKS_REMOTE.keys())
# Placeholder for local tracks logic - for now, we only use remote tracks
# LOCAL_TRACK_DIR = config.get("local_track_dir", os.path.join(appdir, "local_tracks"))
# LOCAL_TRACK_FILES = ... scan LOCAL_TRACK_DIR ...
# ALL_TRACK_NAMES.extend(list(LOCAL_TRACK_FILES.keys()))

# Define basename for non-coding BED file, used in several rules
NONCODING_BED_BASENAME = config.get("noncoding_bed_basename", "new-non-coding.chr")


rule all:
    input:
        # Genome preparation outputs
        f"{appdir}/genome/genome.fa",
        f"{appdir}/genome/genome.gtf",
        # Actual index files from build_genome_index
        expand(f"{appdir}/genome/index/genome_index.{{i}}.ht2", i=range(1,9)),
        f"{appdir}/genome/index/genome.ss", # Splice sites file
        # SRA info outputs (only for samples specified as SRA in config)
        expand(f"{workdir}/sra_info/{{sample}}.info", sample=[s for s in ALL_SAMPLES if "sra" in config["samples_config"][s]]),
        # Mapping outputs
        expand(f"{workdir}/mapped/{{sample}}.bam", sample=ALL_SAMPLES),
        # Assembly outputs
        f"{workdir}/assembly/merged/merged_transcripts.gtf",
        f"{workdir}/cov/ballgown_paths.txt",
        expand(f"{workdir}/cov/{{sample}}/ballgown_input.done", sample=ALL_SAMPLES),
        # Coding potential outputs
        f"{workdir}/results/new-coding.gtf",
        f"{workdir}/results/new-non-coding.gtf",
        # DGE outputs
        f"{workdir}/dge_analysis_complete.txt",
        # Feature extraction outputs
        f"{workdir}/results/non-coding/features/paths.csv",
        # Prediction outputs
        f"{workdir}/results/non-coding_{NONCODING_BED_BASENAME}_NONCODE_PREDICTED.bed"


rule prepare_genome:
    output:
        fasta=f"{appdir}/genome/genome.fa",
        gtf=f"{appdir}/genome/genome.gtf",
    conda:
        "env/infoMod.yml"
    shell:
        # Ensure the output directory for the script is the parent of one of its outputs
        # output.fasta is f"{appdir}/genome/genome.fa", so output.fasta.parent is f"{appdir}/genome"
        """
        python scripts/get_genome.py --output_dir {output.fasta.parent}
        """

rule get_sra_info:
    input:
        # sra_list is not directly used by the script, but ensures SRA_ACCESSIONS is populated
        # The script takes individual sra_accession via wildcards
        # No explicit input file needed here if SRA_ACCESSIONS is used by expand in `all` rule
        # This rule is now targeted by specific SRA samples defined in config
    output:
        sra_info=f"{workdir}/sra_info/{{sample}}.info" # Wildcard is now 'sample'
    conda:
        "env/infoMod.yml"
    shell:
        # The script expects sra_accession, which is derived from the sample's config
        """
        python scripts/get_sra_info.py --workdir {workdir} --sra_accession {config[samples_config][wildcards.sample][sra]}
        """

rule build_genome_index:
    input:
        genome_fasta=f"{appdir}/genome/genome.fa", # Use ancient if this rule is run independently often
        genome_gtf=f"{appdir}/genome/genome.gtf"    # GTF needed for splice sites
    output:
        # HISAT2 creates 8 index files (prefix.1.ht2 to prefix.8.ht2)
        # Touch a marker file, or list all files. Listing is more explicit.
        index_files=expand(f"{appdir}/genome/index/genome_index.{{i}}.ht2", i=range(1,9)),
        splice_site_file=f"{appdir}/genome/index/genome.ss"
    params:
        index_prefix=f"{appdir}/genome/index/genome_index",
        threads=config.get("threads", 1)
    conda:
        "env/mapMod.yml"
    shell:
        """
        python scripts/mapping.py build_index \
            --genome_fasta {input.genome_fasta} \
            --gtf_file {input.genome_gtf} \
            --index_prefix {params.index_prefix} \
            --threads {params.threads}
        """

def get_map_reads_inputs(wildcards):
    """Determine inputs for map_reads based on sample config."""
    sample_conf = config["samples_config"][wildcards.sample]
    inputs = {
        "index_files": expand(f"{appdir}/genome/index/genome_index.{{i}}.ht2", i=range(1,9)),
        "splice_site_file": f"{appdir}/genome/index/genome.ss" # ensure splice site file is an input
    }
    if "sra" in sample_conf:
        # Depends on SRA info if we need to fetch layout or other details from it later
        # For now, assume layout is in config. sra_info rule produces {workdir}/sra_info/{sample}.info
        inputs["sra_info_file"] = f"{workdir}/sra_info/{wildcards.sample}.info"
    # No explicit FASTQ path input here, they are passed as params from config to the script
    return inputs

rule map_reads:
    input:
        unpack(get_map_reads_inputs)
    output:
        bam=f"{workdir}/mapped/{{sample}}.bam"
    params:
        index_prefix=f"{appdir}/genome/index/genome_index",
        threads=config.get("threads", 1),
        sample_config=lambda wildcards: config["samples_config"][wildcards.sample],
        workdir_sra_data=f"{workdir}/sra_downloads" # Centralized SRA download location
    conda:
        "env/mapMod.yml"
    shell:
        """
        cmd_args="map_reads \
            --index_prefix {params.index_prefix} \
            --output_bam {output.bam} \
            --threads {params.threads} \
            --sample_name {wildcards.sample} \
            --layout {params.sample_config[layout]}"

        if [[ -n "{params.sample_config.get(sra, '')}" ]]; then
            cmd_args="$cmd_args --sra_accession {params.sample_config[sra]} \
                                --workdir_sra_data {params.workdir_sra_data}"
        elif [[ -n "{params.sample_config.get(fastq1, '')}" ]]; then
            cmd_args="$cmd_args --fastq1 {params.sample_config[fastq1]}"
            if [[ -n "{params.sample_config.get(fastq2, '')}" ]]; then
                cmd_args="$cmd_args --fastq2 {params.sample_config[fastq2]}"
            fi
        else
            echo "Error: Sample {wildcards.sample} has no sra or fastq1 defined in config."
            exit 1
        fi

        python scripts/mapping.py $cmd_args
        """

rule assemble_transcripts:
    input:
        bam=f"{workdir}/mapped/{{sample}}.bam",
        ref_gtf=f"{appdir}/genome/genome.gtf"
    output:
        gtf=f"{workdir}/assembly/{{sample}}/transcripts.gtf",
        cov=f"{workdir}/assembly/{{sample}}/coverage.cov"
    params:
        threads=config.get("threads", 1)
    conda:
        "env/assembleMod.yml"
    shell:
        """
        python scripts/assembly.py assemble \
            --bam_file {input.bam} \
            --ref_gtf {input.ref_gtf} \
            --gtf_output {output.gtf} \
            --cov_output {output.cov} \
            --threads {params.threads}
        """

rule merge_transcripts:
    input:
        # Gather all per-sample assembled GTFs
        sample_gtfs=expand(f"{workdir}/assembly/{{sample}}/transcripts.gtf", sample=ALL_SAMPLES),
        ref_gtf=f"{appdir}/genome/genome.gtf"
    output:
        merged_gtf=f"{workdir}/assembly/merged/merged_transcripts.gtf",
        ballgown_paths=f"{workdir}/cov/ballgown_paths.txt"
    params:
        threads=config.get("threads", 1),
        workdir=workdir, # Pass the main workdir
        assembly_list_file=f"{workdir}/assembly/assembly_list.txt" # Temporary file
    conda:
        "env/assembleMod.yml"
    shell:
        """
        # Create the list of GTF files for StringTie --merge
        (echo '{input.sample_gtfs}' | tr ' ' '\\n' > {params.assembly_list_file}) && \
        python scripts/assembly.py merge \
            --assembly_list_file {params.assembly_list_file} \
            --ref_gtf {input.ref_gtf} \
            --merged_gtf_output {output.merged_gtf} \
            --threads {params.threads} \
            --workdir {params.workdir}
        """

rule prepare_ballgown_inputs:
    input:
        bam=f"{workdir}/mapped/{{sample}}.bam",
        merged_gtf=f"{workdir}/assembly/merged/merged_transcripts.gtf"
    output:
        # StringTie -B creates multiple files in the directory.
        # We'll use a done file to mark completion for this sample.
        # The directory itself can also be an output.
        # directory(f"{workdir}/cov/{{sample}}"), # This tells Snakemake to expect a directory
        done_file=f"{workdir}/cov/{{sample}}/ballgown_input.done"
    params:
        threads=config.get("threads", 1),
        ballgown_dir=f"{workdir}/cov/{{sample}}"
    conda:
        "env/assembleMod.yml"
    shell:
        """
        python scripts/assembly.py count \
            --bam_file {input.bam} \
            --merged_gtf {input.merged_gtf} \
            --ballgown_dir_output {params.ballgown_dir} \
            --threads {params.threads} && \
        touch {output.done_file}
        """

rule predict_coding_potential:
    input:
        merged_gtf=f"{workdir}/assembly/merged/merged_transcripts.gtf",
        genome_fasta=f"{appdir}/genome/genome.fa",
        # Assuming CPAT static files are in appdir/static/ as per original script
        # These should ideally be specified in config.yaml if paths can vary
        cpat_hexamer_table=f"{appdir}/static/fly_Hexamer.tsv",
        cpat_logit_model=f"{appdir}/static/Fly_logitModel.RData"
    output:
        # Key CPAT output files used by the classification step
        cpat_best_orf_prob=f"{workdir}/cpat/cpat.ORF_prob.best.tsv",
        cpat_no_orf=f"{workdir}/cpat/cpat.no_ORF.txt",
        # Other CPAT files that are created, can be listed if needed by other rules
        # cpat_rRNA=f"{workdir}/cpat/cpat.rRNA.txt",
        # cpat_CDS=f"{workdir}/cpat/cpat.CDS.txt",
    params:
        cpat_output_dir=f"{workdir}/cpat",
        cpat_output_prefix=f"{workdir}/cpat/cpat", # Prefix for cpat.py -o
        threads=config.get("threads", 1) # CPAT itself is single-threaded, for future use
    conda:
        "env/codMod.yml"
    shell:
        """
        python scripts/coding_probability.py predict_coding_potential \
            --merged_gtf {input.merged_gtf} \
            --genome_fasta {input.genome_fasta} \
            --cpat_hexamer_table {input.cpat_hexamer_table} \
            --cpat_logit_model {input.cpat_logit_model} \
            --cpat_output_dir {params.cpat_output_dir} \
            --cpat_output_prefix {params.cpat_output_prefix} \
            --threads {params.threads}
        """

rule classify_transcripts:
    input:
        cpat_best_orf_prob=f"{workdir}/cpat/cpat.ORF_prob.best.tsv",
        cpat_no_orf=f"{workdir}/cpat/cpat.no_ORF.txt",
        merged_gtf=f"{workdir}/assembly/merged/merged_transcripts.gtf"
    output:
        coding_gtf=f"{workdir}/results/new-coding.gtf",
        noncoding_gtf=f"{workdir}/results/new-non-coding.gtf"
    params:
        # cpat_prob_cutoff can be added to config if needed
        cpat_prob_cutoff=config.get("cpat_prob_cutoff", 0.39)
    conda:
        "env/codMod.yml" # Assuming pandas is in codMod or base env
    shell:
        """
        python scripts/coding_probability.py classify_transcripts \
            --cpat_best_orf_prob_file {input.cpat_best_orf_prob} \
            --cpat_no_orf_file {input.cpat_no_orf} \
            --merged_gtf_file {input.merged_gtf} \
            --output_coding_gtf {output.coding_gtf} \
            --output_noncoding_gtf {output.noncoding_gtf} \
            --cpat_prob_cutoff {params.cpat_prob_cutoff}
        """

rule run_ballgown_dge:
    input:
        # Depends on all samples having their Ballgown tables ready
        ballgown_sample_data=expand(f"{workdir}/cov/{{sample}}/ballgown_input.done", sample=ALL_SAMPLES),
        ballgown_paths_file=f"{workdir}/cov/ballgown_paths.txt", # List of sample dirs for Ballgown
        # Metadata file from config, could be None or empty if DGE is not run
        metadata_file=config.get("metadata_file", "")
    output:
        # Key output tables from the R script (adjust paths if R script changes)
        # The R script currently writes to workdir/results/
        transcript_results=f"{workdir}/results/dge_transcripts.csv",
        gene_results=f"{workdir}/results/dge_exons.csv", # Original R script calls this dge_exons.csv
        # Marker file for completion
        dge_done=f"{workdir}/dge_analysis_complete.txt"
    params:
        workdir=workdir,
        r_script_path=f"{appdir}/scripts/ballgown.R"
    log:
        f"{workdir}/logs/ballgown_dge.log"
    conda:
        "env/dgeMod.yml"
    shell:
        """
        # Check if metadata_file is provided and exists
        if [ -z "{input.metadata_file}" ] || [ ! -f "{input.metadata_file}" ]; then
            echo "Metadata file not provided or not found ({input.metadata_file}). Skipping DGE analysis." > {log}
            # Create dummy output files and the done marker to satisfy Snakemake
            touch {output.transcript_results} {output.gene_results} {output.dge_done}
        else
            # Ensure the results directory exists (R script might not create it if it expects workdir/results)
            mkdir -p {workdir}/results
            Rscript {params.r_script_path} {params.workdir} {input.metadata_file} > {log} 2>&1 && \
            touch {output.dge_done}
        fi
        """

rule gtf_to_bed:
    input:
        noncoding_gtf=f"{workdir}/results/new-non-coding.gtf"
    output:
        # Outputting to workdir/results as per plan, original script was less specific
        noncoding_bed=f"{workdir}/results/new-non-coding.chr.bed"
    conda:
        "env/infoMod.yml" # bedops (convert2bed, sort-bed) assumed here
    shell:
        """
        python scripts/feature_extraction.py gtf_to_bed \
            --input_gtf {input.noncoding_gtf} \
            --output_bed {output.noncoding_bed}
        """

rule extract_single_feature:
    input:
        bed_file=f"{workdir}/results/new-non-coding.chr.bed",
        # Dynamically get track_path based on wildcard track_name
        track_path=lambda wildcards: TRACKS_REMOTE[wildcards.track_name]
    output:
        feature_tsv=f"{workdir}/results/non-coding/features/{{track_name}}.tsv"
    params:
        track_name="{track_name}", # Pass wildcard to script
        output_dir=f"{workdir}/results/non-coding/features"
    conda:
        "env/featureMod.yml" # UCSC utils (bigWigAverageOverBed etc.)
    shell:
        """
        python scripts/feature_extraction.py extract_single_feature \
            --track_name {params.track_name} \
            --track_path '{input.track_path}' \
            --input_bed {input.bed_file} \
            --output_dir {params.output_dir}
        """

rule create_feature_paths_csv:
    input:
        # Depends on all feature TSVs being created
        feature_tsvs=expand(f"{workdir}/results/non-coding/features/{{track_name}}.tsv", track_name=ALL_TRACK_NAMES)
    output:
        paths_csv=f"{workdir}/results/non-coding/features/paths.csv"
    params:
        feature_dir=f"{workdir}/results/non-coding/features"
    conda:
        "env/infoMod.yml" # pandas for df.to_csv
    shell:
        """
        python scripts/feature_extraction.py create_paths_csv \
            --feature_dir {params.feature_dir} \
            --output_csv {output.paths_csv}
        """

rule create_feature_table:
    input:
        paths_csv=f"{workdir}/results/non-coding/features/paths.csv",
        noncoding_bed=f"{workdir}/results/{NONCODING_BED_BASENAME}.bed"
    output:
        feature_table=f"{workdir}/results/non-coding/final_table_transcripts_NONCODE_FEATURES.tsv"
    params:
        appdir=appdir,
        workdir=workdir,
        # The script feature-table.py expects the full path to the bed file as the third argument
        bed_file_path=lambda wildcards, input: input.noncoding_bed
    conda:
        "env/predictMod.yml"
    shell:
        "python {params.appdir}/scripts/feature-table.py {params.appdir} {params.workdir} {params.bed_file_path}"

rule run_prediction:
    input:
        feature_table=f"{workdir}/results/non-coding/final_table_transcripts_NONCODE_FEATURES.tsv",
        # Model path should be configurable, but using appdir path as per subtask description
        model=f"{appdir}/model/rf_dm6_lncrna_classifier.model"
    output:
        predictions_tsv=f"{workdir}/results/non-coding/rf_predictions_transcripts_NONCODE_FEATURES.tsv"
    params:
        appdir=appdir,
        workdir=workdir,
        # The script predict.py expects the basename of the BED file (without .bed) as the third argument
        bed_basename=NONCODING_BED_BASENAME
    conda:
        "env/predictMod.yml"
    shell:
        "python {params.appdir}/scripts/predict.py {params.appdir} {params.workdir} {params.bed_basename}"

rule create_final_output_table:
    input:
        predictions_tsv=f"{workdir}/results/non-coding/rf_predictions_transcripts_NONCODE_FEATURES.tsv",
        noncoding_bed=f"{workdir}/results/{NONCODING_BED_BASENAME}.bed"
    output:
        final_bed=f"{workdir}/results/non-coding_{NONCODING_BED_BASENAME}_NONCODE_PREDICTED.bed"
    params:
        appdir=appdir,
        workdir=workdir,
        # The script final-table.py expects the full path to the bed file as the third argument
        bed_file_path=lambda wildcards, input: input.noncoding_bed
    conda:
        "env/predictMod.yml"
    shell:
        "python {params.appdir}/scripts/final-table.py {params.appdir} {params.workdir} {params.bed_file_path}"
