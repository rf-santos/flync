import argparse
import os
import subprocess
import sys
import shutil

def run_command(command, cwd=None):
    """Helper function to run a shell command and print output/error."""
    print(f"Running command: {' '.join(command)}", flush=True)
    try:
        process = subprocess.run(command, check=True, text=True, capture_output=True, cwd=cwd)
        if process.stdout:
            print("STDOUT:\n", process.stdout, flush=True)
        if process.stderr:
            # HISAT2 often writes informational messages to stderr
            print("STDERR:\n", process.stderr, flush=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {' '.join(e.cmd)}", flush=True)
        print(f"Return code: {e.returncode}", flush=True)
        if e.stdout:
            print(f"STDOUT:\n{e.stdout}", flush=True)
        if e.stderr:
            print(f"STDERR:\n{e.stderr}", flush=True)
        return False
    except Exception as e:
        print(f"An unexpected error occurred: {e}", flush=True)
        return False

def build_index(genome_fasta, gtf_file, index_prefix, threads):
    """Builds HISAT2 genome index and extracts splice sites."""
    print(f"Building HISAT2 index with prefix: {index_prefix}")

    index_dir = os.path.dirname(index_prefix)
    if index_dir: # Ensure index directory exists
        os.makedirs(index_dir, exist_ok=True)

    # 1. Build HISAT2 index
    cmd_build = [
        "hisat2-build",
        "-p", str(threads),
        genome_fasta,
        index_prefix
    ]
    if not run_command(cmd_build):
        print("HISAT2 index building failed.")
        return False

    # 2. Extract splice sites
    splice_sites_file = os.path.join(index_dir, "genome.ss") # Store it with the index
    print(f"Extracting splice sites to: {splice_sites_file}")
    cmd_extract_ss = [
        "hisat2_extract_splice_sites.py",
        gtf_file
    ]
    # Need to capture stdout to a file
    try:
        print(f"Running command: {' '.join(cmd_extract_ss)} > {splice_sites_file}", flush=True)
        with open(splice_sites_file, "w") as f_out:
            process = subprocess.run(cmd_extract_ss, check=True, text=True, stdout=f_out, stderr=subprocess.PIPE)
            if process.stderr:
                 print("STDERR (hisat2_extract_splice_sites.py):\n", process.stderr, flush=True)
        print("Splice site extraction complete.")
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {' '.join(e.cmd)}", flush=True)
        print(f"Return code: {e.returncode}", flush=True)
        if e.stdout: # Should be in file, but just in case
            print(f"STDOUT:\n{e.stdout}", flush=True)
        if e.stderr:
            print(f"STDERR:\n{e.stderr}", flush=True)
        return False

    print("Genome indexing and splice site extraction complete.")
    return True

def map_reads(index_prefix, output_bam, threads, layout,
              sra_accession=None, fastq1=None, fastq2=None,
              sample_name=None, workdir_sra_data=None):
    """Maps reads using HISAT2 and processes to sorted BAM."""

    if not sample_name:
        sample_name = sra_accession if sra_accession else os.path.basename(fastq1).split('_')[0]

    os.makedirs(os.path.dirname(output_bam), exist_ok=True)

    # Determine paths for temporary files
    sam_output_path = output_bam + ".sam"
    unsorted_bam_path = output_bam + ".unsorted.bam"

    # Path to the splice sites file (assuming it's with the index)
    splice_sites_file = os.path.join(os.path.dirname(index_prefix), "genome.ss")
    if not os.path.exists(splice_sites_file):
        print(f"Error: Splice sites file {splice_sites_file} not found. Run build_index first.")
        return False

    hisat2_cmd = [
        "hisat2",
        "-p", str(threads),
        "-x", index_prefix,
        "--dta", "--dta-cufflinks", # Standard options from original scripts
        "--known-splicesite-infile", splice_sites_file,
        "-S", sam_output_path
    ]

    local_fastq_files_to_remove = []

    if sra_accession:
        if not workdir_sra_data:
            print("Error: workdir_sra_data must be provided for SRA processing.")
            return False

        sra_sample_dir = os.path.join(workdir_sra_data, sample_name)
        os.makedirs(sra_sample_dir, exist_ok=True)

        print(f"Downloading SRA {sra_accession} to {sra_sample_dir}...")
        # Using fasterq-dump, preferred over fastq-dump
        # Ensure fasterq-dump is in PATH
        # We need to handle splitting for paired-end data
        fastq_dump_cmd = [
            "fasterq-dump",
            "--split-files", # Creates _1.fastq and _2.fastq if paired
            "-O", sra_sample_dir,
            "-e", str(threads),
            "-f", # force overwrite
            "-3", # output R1/R2/single reads in new FASTQ standard format
            sra_accession
        ]
        if not run_command(fastq_dump_cmd, cwd=sra_sample_dir): # Run in sample dir to avoid clutter
            print(f"fasterq-dump failed for {sra_accession}.")
            return False

        # Define expected fastq paths after download
        # fasterq-dump names files {accession}_1.fastq, {accession}_2.fastq, or {accession}.fastq
        # It does not automatically gzip them.
        sra_fastq1_raw = os.path.join(sra_sample_dir, f"{sra_accession}_1.fastq")
        sra_fastq2_raw = os.path.join(sra_sample_dir, f"{sra_accession}_2.fastq")
        sra_fastq_single_raw = os.path.join(sra_sample_dir, f"{sra_accession}.fastq")

        # Gzip the outputs of fasterq-dump
        print("Gzipping downloaded FASTQ files...")
        sra_fastq1_gz = ""
        sra_fastq2_gz = ""

        if layout.upper() == "PAIRED":
            if os.path.exists(sra_fastq1_raw) and os.path.exists(sra_fastq2_raw):
                sra_fastq1_gz = sra_fastq1_raw + ".gz"
                sra_fastq2_gz = sra_fastq2_raw + ".gz"
                if not run_command(["gzip", "-f", sra_fastq1_raw]): return False
                if not run_command(["gzip", "-f", sra_fastq2_raw]): return False
                hisat2_cmd.extend(["-1", sra_fastq1_gz, "-2", sra_fastq2_gz])
                local_fastq_files_to_remove.extend([sra_fastq1_gz, sra_fastq2_gz])
            else:
                print(f"Error: Paired-end SRA download did not produce expected _1 and _2 files for {sra_accession}")
                return False
        else: # SINGLE
            if os.path.exists(sra_fastq_single_raw):
                sra_fastq1_gz = sra_fastq_single_raw + ".gz"
                if not run_command(["gzip", "-f", sra_fastq_single_raw]): return False
                hisat2_cmd.extend(["-U", sra_fastq1_gz])
                local_fastq_files_to_remove.append(sra_fastq1_gz)
            else:
                print(f"Error: Single-end SRA download did not produce expected .fastq file for {sra_accession}")
                return False

    elif fastq1: # Local FASTQ files
        if layout.upper() == "PAIRED":
            if not fastq2:
                print("Error: Paired-end layout specified but fastq2 is missing.")
                return False
            hisat2_cmd.extend(["-1", fastq1, "-2", fastq2])
        else: # SINGLE
            hisat2_cmd.extend(["-U", fastq1])
    else:
        print("Error: No input specified (neither SRA accession nor FASTQ files).")
        return False

    # Run HISAT2
    print("Starting HISAT2 mapping...")
    if not run_command(hisat2_cmd):
        print("HISAT2 mapping failed.")
        return False

    # Samtools processing: view SAM -> BAM, sort BAM, index BAM
    print("Converting SAM to BAM...")
    samtools_view_cmd = ["samtools", "view", "-@_threads_placeholder_", "-bS", "-o", unsorted_bam_path, sam_output_path]
    samtools_view_cmd = [s.replace("_threads_placeholder_", str(max(1, threads // 2))) for s in samtools_view_cmd]

    if not run_command(samtools_view_cmd):
        print("samtools view failed.")
        # Clean up SAM if view fails to prevent issues on retry
        if os.path.exists(sam_output_path): os.remove(sam_output_path)
        return False

    # Remove SAM file after successful conversion
    if os.path.exists(sam_output_path):
        os.remove(sam_output_path)

    print("Sorting BAM...")
    samtools_sort_cmd = ["samtools", "sort", "-@_threads_placeholder_", "-o", output_bam, unsorted_bam_path]
    samtools_sort_cmd = [s.replace("_threads_placeholder_", str(max(1, threads // 2))) for s in samtools_sort_cmd]
    if not run_command(samtools_sort_cmd):
        print("samtools sort failed.")
        # Clean up unsorted BAM if sort fails
        if os.path.exists(unsorted_bam_path): os.remove(unsorted_bam_path)
        return False

    # Remove unsorted BAM after successful sort
    if os.path.exists(unsorted_bam_path):
        os.remove(unsorted_bam_path)

    # Index the sorted BAM
    # print("Indexing BAM...")
    # samtools_index_cmd = ["samtools", "index", output_bam]
    # if not run_command(samtools_index_cmd):
    #     print("samtools index failed.")
    #     return False
    # Indexing is typically a separate step in Snakemake if the .bai file is an explicit output.
    # For now, the mapping rule just produces the BAM.

    # Clean up downloaded SRA FASTQ files
    for f_path in local_fastq_files_to_remove:
        if os.path.exists(f_path):
            print(f"Removing temporary FASTQ: {f_path}")
            os.remove(f_path)
    # Clean up SRA sample directory if it's empty (or contains only logs etc.)
    if sra_accession and workdir_sra_data:
        sra_sample_dir = os.path.join(workdir_sra_data, sample_name)
        try:
            if not os.listdir(sra_sample_dir): # Check if empty
                os.rmdir(sra_sample_dir)
            # else: # Optionally remove all contents of sra_sample_dir if it's truly temporary
            #    shutil.rmtree(sra_sample_dir)
        except OSError as e:
            print(f"Could not remove or check SRA sample directory {sra_sample_dir}: {e}", flush=True)


    print(f"Mapping for {sample_name} complete. Output BAM: {output_bam}")
    return True


def main():
    parser = argparse.ArgumentParser(description="HISAT2 indexing and mapping script.")
    subparsers = parser.add_subparsers(dest="action", required=True)

    # Subparser for build_index
    parser_idx = subparsers.add_parser("build_index", help="Build HISAT2 genome index.")
    parser_idx.add_argument("--genome_fasta", required=True, help="Path to genome FASTA file.")
    parser_idx.add_argument("--gtf_file", required=True, help="Path to genome GTF file (for splice sites).")
    parser_idx.add_argument("--index_prefix", required=True, help="Prefix for HISAT2 index files (e.g., path/to/genome_index).")
    parser_idx.add_argument("--threads", type=int, default=1, help="Number of threads to use.")

    # Subparser for map_reads
    parser_map = subparsers.add_parser("map_reads", help="Map reads with HISAT2 and process to BAM.")
    parser_map.add_argument("--index_prefix", required=True, help="Prefix of HISAT2 index files.")
    parser_map.add_argument("--output_bam", required=True, help="Path for the final sorted BAM file.")
    parser_map.add_argument("--threads", type=int, default=1, help="Number of threads for mapping and samtools.")
    parser_map.add_argument("--layout", required=True, choices=["SINGLE", "PAIRED"], help="Read layout (SINGLE or PAIRED).")

    # Input options for map_reads
    group_input = parser_map.add_mutually_exclusive_group(required=True)
    group_input.add_argument("--sra_accession", help="SRA accession number. Will use fasterq-dump.")
    group_input.add_argument("--fastq1", help="Path to FASTQ file (R1 for paired-end, or single-end).")

    parser_map.add_argument("--fastq2", help="Path to R2 FASTQ file (for paired-end).")
    parser_map.add_argument("--sample_name", help="Sample name for output structuring (defaults to SRA or fastq1 prefix).")
    parser_map.add_argument("--workdir_sra_data", help="Directory for storing intermediate SRA downloads (e.g. workdir/sra_downloads). Required if --sra_accession is used.")


    args = parser.parse_args()

    success = False
    if args.action == "build_index":
        success = build_index(args.genome_fasta, args.gtf_file, args.index_prefix, args.threads)
    elif args.action == "map_reads":
        if args.sra_accession and not args.workdir_sra_data:
            parser_map.error("--workdir_sra_data is required when using --sra_accession.")
        if args.fastq1 and args.layout.upper() == "PAIRED" and not args.fastq2:
            parser_map.error("--fastq2 is required for PAIRED layout with local FASTQ files.")
        success = map_reads(
            args.index_prefix, args.output_bam, args.threads, args.layout,
            args.sra_accession, args.fastq1, args.fastq2,
            args.sample_name, args.workdir_sra_data
        )

    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()
