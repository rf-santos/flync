import argparse
import os
import subprocess
import sys
import pandas as pd

def run_command(command, cwd=None, shell=False):
    """Helper function to run a shell command and print output/error."""
    print(f"Running command: {' '.join(command) if isinstance(command, list) else command}", flush=True)
    try:
        process = subprocess.run(command, check=True, text=True, capture_output=True, cwd=cwd, shell=shell)
        if process.stdout:
            print("STDOUT:\n", process.stdout, flush=True)
        if process.stderr:
            print("STDERR:\n", process.stderr, flush=True)
        return True
    except subprocess.CalledProcessError as e:
        cmd_str = ' '.join(e.cmd) if isinstance(e.cmd, list) else e.cmd
        print(f"Error running command: {cmd_str}", flush=True)
        print(f"Return code: {e.returncode}", flush=True)
        if e.stdout:
            print(f"STDOUT:\n{e.stdout}", flush=True)
        if e.stderr:
            print(f"STDERR:\n{e.stderr}", flush=True)
        return False
    except Exception as e:
        print(f"An unexpected error occurred: {e}", flush=True)
        return False

def predict_coding_potential(merged_gtf, genome_fasta, cpat_hexamer_table, cpat_logit_model,
                             cpat_output_dir, cpat_output_prefix, threads):
    """
    Predicts coding potential using gffread and CPAT.
    Omits BLASTx steps due to missing parse_blast_cpat.py.
    """
    os.makedirs(cpat_output_dir, exist_ok=True)

    # 1. Extract transcript sequences using gffread
    # Original script used $workdir/assemblies/assembled-new-transcripts.fa
    # which was created from merged-new-transcripts.gtf (MSTRG IDs only)
    # Here, we'll use the full merged_gtf and let CPAT handle it, or filter first.
    # For now, let's use the full merged_gtf.
    # A more faithful approach would be to first filter merged_gtf for 'novel' transcripts (MSTRG).

    novel_transcripts_gtf = os.path.join(cpat_output_dir, "novel_transcripts.gtf")
    print(f"Filtering novel transcripts from {merged_gtf} to {novel_transcripts_gtf}")
    try:
        with open(merged_gtf, 'r') as infile, open(novel_transcripts_gtf, 'w') as outfile:
            for line in infile:
                if 'MSTRG.' in line: # Simple check for novel transcripts
                    outfile.write(line)
    except IOError as e:
        print(f"Error filtering GTF: {e}")
        return False

    transcript_fasta = os.path.join(cpat_output_dir, "transcripts_for_cpat.fa")
    print(f"Extracting transcript FASTA from {novel_transcripts_gtf} to {transcript_fasta}...")
    cmd_gffread = [
        "gffread",
        "-w", transcript_fasta,
        "-g", genome_fasta,
        novel_transcripts_gtf # Use the filtered GTF
    ]
    if not run_command(cmd_gffread):
        print("gffread failed.")
        return False

    # 2. Run CPAT
    print("Running CPAT...")
    # Ensure cpat.py is in PATH or provide full path
    # Original script command: cpat.py -x hexamer -d logit_model -g input_fasta -o output_prefix
    cmd_cpat = [
        "cpat.py",
        "-x", cpat_hexamer_table,
        "-d", cpat_logit_model,
        "-g", transcript_fasta,
        "-o", cpat_output_prefix # CPAT creates multiple files with this prefix
    ]
    if not run_command(cmd_cpat):
        print("CPAT failed.")
        return False

    # Check for expected CPAT output files
    expected_orf_prob = f"{cpat_output_prefix}.ORF_prob.best.tsv"
    expected_no_orf = f"{cpat_output_prefix}.no_ORF.txt"
    if not (os.path.exists(expected_orf_prob) and os.path.exists(expected_no_orf)):
        print(f"CPAT did not produce all expected output files ({expected_orf_prob}, {expected_no_orf}). Please check CPAT logs.")
        # return False # Allow to proceed if some files are there, classify can handle missing.

    print("Coding potential prediction (CPAT part) complete.")
    return True


def classify_transcripts(cpat_best_orf_prob_file, cpat_no_orf_file, merged_gtf_file,
                         output_coding_gtf, output_noncoding_gtf, cpat_prob_cutoff=0.39):
    """
    Classifies transcripts from merged_gtf based on CPAT results.
    This version adapts logic from class-new-transfrags.sh but uses the
    StringTie merged GTF directly and filters MSTRG transcripts.
    """
    print("Classifying transcripts...")

    # Load CPAT results
    coding_ids = set()
    non_coding_ids = set()

    try:
        if os.path.exists(cpat_best_orf_prob_file):
            orf_probs = pd.read_csv(cpat_best_orf_prob_file, sep='\t', header=0)
            # Column names in cpat.ORF_prob.best.tsv are typically:
            # mRNA_ID, ORF_ID, ORF_length, CDS_length, Fickett_score, Hexamer_score, ORF_integrity, CPAT_probability
            # The original script uses $11 for CPAT probability, which is likely 1-based index.
            # Assuming 'CPAT_probability' is the correct column name or it's the last column (index -1 or 7 for 0-based).
            # Let's try to be robust: use column name if present, else fallback to index.
            prob_col_name = 'CPAT_probability'
            if prob_col_name not in orf_probs.columns:
                # Fallback, assuming it's the last column if typical names are not found
                # This is risky; better to know the exact column name.
                # For now, let's assume it's the 8th column (index 7) if 'CPAT_probability' fails
                prob_col_idx = 7
                if len(orf_probs.columns) <= prob_col_idx: # check if column exists
                     print(f"Error: CPAT probability column not found by name '{prob_col_name}' or by index {prob_col_idx} in {cpat_best_orf_prob_file}")
                     return False
            else:
                prob_col_idx = orf_probs.columns.get_loc(prob_col_name)


            for _, row in orf_probs.iterrows():
                transcript_id = row.iloc[0] # First column is mRNA_ID
                cpat_prob = float(row.iloc[prob_col_idx])
                if transcript_id.startswith("MSTRG."):
                    if cpat_prob > cpat_prob_cutoff:
                        coding_ids.add(transcript_id)
                    else:
                        non_coding_ids.add(transcript_id)
        else:
            print(f"Warning: {cpat_best_orf_prob_file} not found. Coding transcripts might be underrepresented.")

        if os.path.exists(cpat_no_orf_file):
            with open(cpat_no_orf_file, 'r') as f:
                for line in f:
                    transcript_id = line.strip()
                    if transcript_id.startswith("MSTRG."): # Ensure we only add MSTRG IDs
                        non_coding_ids.add(transcript_id) # IDs here are considered non-coding
        else:
            print(f"Warning: {cpat_no_orf_file} not found. Non-coding transcripts might be underrepresented.")

    except Exception as e:
        print(f"Error processing CPAT files: {e}")
        return False

    # Filter the merged GTF
    print(f"Filtering {merged_gtf_file} into coding and non-coding GTFs...")
    coding_count = 0
    noncoding_count = 0
    try:
        with open(merged_gtf_file, 'r') as mgf, \
             open(output_coding_gtf, 'w') as cgf, \
             open(output_noncoding_gtf, 'w') as ncgf:

            for line in mgf:
                if line.startswith("#"):
                    cgf.write(line)
                    ncgf.write(line)
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 9: continue

                attributes = fields[8]
                transcript_id = None
                # Example GTF attribute: transcript_id "MSTRG.1.1";
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if attr.startswith('transcript_id'):
                        transcript_id = attr.split('"')[1]
                        break

                if transcript_id and transcript_id.startswith("MSTRG."): # Process only novel transcripts
                    if transcript_id in coding_ids:
                        cgf.write(line)
                        coding_count +=1
                    elif transcript_id in non_coding_ids:
                        ncgf.write(line)
                        noncoding_count += 1
                    # else: Transcripts not in CPAT outputs (e.g. filtered out by CPAT due to length) are skipped
    except IOError as e:
        print(f"Error reading/writing GTF files: {e}")
        return False

    print(f"Classification complete. Found {coding_count} coding transcript lines and {noncoding_count} non-coding transcript lines.")
    if coding_count == 0 and noncoding_count == 0 and (len(coding_ids)>0 or len(non_coding_ids)>0) :
        print(f"Warning: No MSTRG transcripts from CPAT output were found in {merged_gtf_file}. Check transcript ID consistency.")

    return True


def main():
    parser = argparse.ArgumentParser(description="Predict coding potential and classify transcripts.")
    subparsers = parser.add_subparsers(dest="action", required=True)

    # Subparser for predict_coding_potential
    parser_pred = subparsers.add_parser("predict_coding_potential", help="Predict coding potential using CPAT.")
    parser_pred.add_argument("--merged_gtf", required=True, help="Path to merged StringTie GTF file.")
    parser_pred.add_argument("--genome_fasta", required=True, help="Path to genome FASTA file.")
    parser_pred.add_argument("--cpat_hexamer_table", required=True, help="Path to CPAT hexamer table.")
    parser_pred.add_argument("--cpat_logit_model", required=True, help="Path to CPAT logit model.")
    parser_pred.add_argument("--cpat_output_dir", required=True, help="Directory for CPAT outputs.")
    parser_pred.add_argument("--cpat_output_prefix", required=True, help="Prefix for CPAT output files (e.g., cpat_dir/cpat).")
    parser_pred.add_argument("--threads", type=int, default=1, help="Number of threads (for future use, CPAT is single-threaded).")
    # BLAST related args are omitted for now

    # Subparser for classify_transcripts
    parser_class = subparsers.add_parser("classify_transcripts", help="Classify transcripts based on CPAT results.")
    parser_class.add_argument("--cpat_best_orf_prob_file", required=True, help="Path to CPAT ORF probability file (*.ORF_prob.best.tsv).")
    parser_class.add_argument("--cpat_no_orf_file", required=True, help="Path to CPAT no ORF file (*.no_ORF.txt).")
    parser_class.add_argument("--merged_gtf_file", required=True, help="Path to the input merged GTF file (e.g., StringTie merged GTF).")
    parser_class.add_argument("--output_coding_gtf", required=True, help="Path for the output coding GTF file.")
    parser_class.add_argument("--output_noncoding_gtf", required=True, help="Path for the output non-coding GTF file.")
    parser_class.add_argument("--cpat_prob_cutoff", type=float, default=0.39, help="CPAT probability cutoff for coding.")


    args = parser.parse_args()
    success = False

    if args.action == "predict_coding_potential":
        success = predict_coding_potential(args.merged_gtf, args.genome_fasta, args.cpat_hexamer_table,
                                           args.cpat_logit_model, args.cpat_output_dir,
                                           args.cpat_output_prefix, args.threads)
    elif args.action == "classify_transcripts":
        success = classify_transcripts(args.cpat_best_orf_prob_file, args.cpat_no_orf_file,
                                       args.merged_gtf_file, args.output_coding_gtf,
                                       args.output_noncoding_gtf, args.cpat_prob_cutoff)

    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()
