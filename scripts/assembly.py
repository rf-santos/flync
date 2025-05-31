import argparse
import os
import subprocess
import sys

def run_command(command, cwd=None):
    """Helper function to run a shell command and print output/error."""
    print(f"Running command: {' '.join(command)}", flush=True)
    try:
        process = subprocess.run(command, check=True, text=True, capture_output=True, cwd=cwd)
        if process.stdout:
            print("STDOUT:\n", process.stdout, flush=True)
        if process.stderr:
            print("STDERR:\n", process.stderr, flush=True) # StringTie often writes info to stderr
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

def assemble(bam_file, ref_gtf, gtf_output, cov_output, threads):
    """Assembles transcripts using StringTie."""
    print(f"Assembling transcripts for {bam_file}...")
    os.makedirs(os.path.dirname(gtf_output), exist_ok=True)
    if cov_output: # Ensure directory for coverage output also exists
        os.makedirs(os.path.dirname(cov_output), exist_ok=True)

    cmd = [
        "stringtie", bam_file,
        "-G", ref_gtf,
        "-o", gtf_output,
        "-p", str(threads),
        "-l", os.path.splitext(os.path.basename(gtf_output))[0] # Use sample name for label
    ]
    if cov_output:
        cmd.extend(["-c", cov_output]) # Output coverage file

    if not run_command(cmd):
        print(f"StringTie assembly failed for {bam_file}.")
        return False
    print("StringTie assembly complete.")
    return True

def merge(assembly_list_file, ref_gtf, merged_gtf_output, threads, workdir):
    """Merges multiple transcript assemblies using StringTie."""
    print(f"Merging transcript assemblies listed in {assembly_list_file}...")
    os.makedirs(os.path.dirname(merged_gtf_output), exist_ok=True)

    cmd_merge = [
        "stringtie", "--merge",
        assembly_list_file,
        "-G", ref_gtf,
        "-o", merged_gtf_output,
        "-p", str(threads)
    ]
    if not run_command(cmd_merge):
        print("StringTie merge failed.")
        return False
    print("StringTie merge complete.")

    # Create ballgown_paths.txt
    # This file should list the subdirectories in workdir/cov/
    # These subdirectories are created by the 'count' step.
    ballgown_base_dir = os.path.join(workdir, "cov")
    ballgown_paths_file = os.path.join(ballgown_base_dir, "ballgown_paths.txt")

    print(f"Creating {ballgown_paths_file}...")
    if not os.path.exists(ballgown_base_dir):
        print(f"Warning: Directory {ballgown_base_dir} does not exist. Cannot create ballgown_paths.txt yet.")
        # It's okay if it doesn't exist if no count rules have run yet.
        # The file will be empty or not created, to be populated as samples are processed by 'count'.
        # However, for the rule structure, we create it, possibly empty.
        os.makedirs(ballgown_base_dir, exist_ok=True) # Ensure cov dir exists

    sample_dirs_for_ballgown = []
    if os.path.exists(ballgown_base_dir):
        for item in os.listdir(ballgown_base_dir):
            item_path = os.path.join(ballgown_base_dir, item)
            if os.path.isdir(item_path):
                # Check if it's a directory that looks like a sample output (e.g., contains .ctab files)
                # For simplicity now, just list all directories.
                # A more robust check might be needed if other dirs are created in workdir/cov.
                sample_dirs_for_ballgown.append(item_path) # Store full path, or relative to workdir/cov

    try:
        with open(ballgown_paths_file, "w") as f:
            for sample_dir_path in sorted(sample_dirs_for_ballgown):
                 # Ballgown typically expects paths relative to where analysis is run,
                 # or full paths. Let's write the directory name.
                f.write(os.path.basename(sample_dir_path) + "\n")
        print(f"Successfully created/updated {ballgown_paths_file}.")
    except IOError as e:
        print(f"Error writing {ballgown_paths_file}: {e}")
        return False

    return True

def count(bam_file, merged_gtf, ballgown_dir_output, threads):
    """Generates Ballgown input files using StringTie."""
    print(f"Generating Ballgown inputs for {bam_file} in {ballgown_dir_output}...")
    os.makedirs(ballgown_dir_output, exist_ok=True)

    # StringTie with -B creates files in the output directory.
    # The -o path is for the *covered transcripts GTF reference* for this sample.
    # e.g. workdir/cov/sampleA/sampleA.gtf
    output_gtf_in_ballgown_dir = os.path.join(ballgown_dir_output, os.path.basename(ballgown_dir_output) + ".gtf")

    cmd = [
        "stringtie", bam_file,
        "-G", merged_gtf,
        "-o", output_gtf_in_ballgown_dir,
        "-eB", # -e: only estimate expression of transcripts in merged_gtf; -B: create Ballgown tables
        "-p", str(threads)
    ]
    if not run_command(cmd):
        print(f"StringTie count failed for {bam_file}.")
        return False
    print("StringTie count (Ballgown input generation) complete.")
    return True

def main():
    parser = argparse.ArgumentParser(description="StringTie operations: assemble, merge, count.")
    subparsers = parser.add_subparsers(dest="action", required=True)

    # Subparser for assemble
    parser_asm = subparsers.add_parser("assemble", help="Assemble transcripts with StringTie.")
    parser_asm.add_argument("--bam_file", required=True, help="Path to input sorted BAM file.")
    parser_asm.add_argument("--ref_gtf", required=True, help="Path to reference GTF file.")
    parser_asm.add_argument("--gtf_output", required=True, help="Path for output assembled GTF file.")
    parser_asm.add_argument("--cov_output", help="Path for output coverage file (optional).")
    parser_asm.add_argument("--threads", type=int, default=1, help="Number of threads.")

    # Subparser for merge
    parser_mrg = subparsers.add_parser("merge", help="Merge transcript assemblies with StringTie.")
    parser_mrg.add_argument("--assembly_list_file", required=True, help="File containing a list of GTF paths to merge.")
    parser_mrg.add_argument("--ref_gtf", required=True, help="Path to reference GTF file.")
    parser_mrg.add_argument("--merged_gtf_output", required=True, help="Path for the merged output GTF file.")
    parser_mrg.add_argument("--threads", type=int, default=1, help="Number of threads.")
    parser_mrg.add_argument("--workdir", required=True, help="Main working directory (to find workdir/cov).")


    # Subparser for count
    parser_cnt = subparsers.add_parser("count", help="Generate Ballgown inputs with StringTie.")
    parser_cnt.add_argument("--bam_file", required=True, help="Path to input sorted BAM file.")
    parser_cnt.add_argument("--merged_gtf", required=True, help="Path to merged GTF file (from stringtie --merge).")
    parser_cnt.add_argument("--ballgown_dir_output", required=True, help="Directory for Ballgown output files (e.g., workdir/cov/sample_name).")
    parser_cnt.add_argument("--threads", type=int, default=1, help="Number of threads.")

    args = parser.parse_args()

    success = False
    if args.action == "assemble":
        success = assemble(args.bam_file, args.ref_gtf, args.gtf_output, args.cov_output, args.threads)
    elif args.action == "merge":
        success = merge(args.assembly_list_file, args.ref_gtf, args.merged_gtf_output, args.threads, args.workdir)
    elif args.action == "count":
        success = count(args.bam_file, args.merged_gtf, args.ballgown_dir_output, args.threads)

    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()
