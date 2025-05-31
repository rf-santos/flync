import argparse
import os
import subprocess
import sys
import pandas as pd # For gtf_to_bed if needed, and create_paths_csv

def run_command(command_list, cwd=None, shell=False, check=True):
    """Helper function to run a shell command."""
    command_str = ' '.join(command_list) if isinstance(command_list, list) else command_list
    print(f"Running command: {command_str}", flush=True)
    try:
        process = subprocess.run(command_list, check=check, text=True, capture_output=True, cwd=cwd, shell=shell)
        if process.stdout and process.stdout.strip():
            print("STDOUT:\n", process.stdout.strip(), flush=True)
        if process.stderr and process.stderr.strip():
            print("STDERR:\n", process.stderr.strip(), flush=True)
        return process
    except subprocess.CalledProcessError as e:
        cmd_str_err = ' '.join(e.cmd) if isinstance(e.cmd, list) else e.cmd
        print(f"Error running command: {cmd_str_err}", flush=True)
        print(f"Return code: {e.returncode}", flush=True)
        if e.stdout and e.stdout.strip():
            print(f"STDOUT:\n{e.stdout.strip()}", flush=True)
        if e.stderr and e.stderr.strip():
            print(f"STDERR:\n{e.stderr.strip()}", flush=True)
        if check:
            raise
        return None # Should not be reached if check=True
    except Exception as e:
        print(f"An unexpected error occurred: {e}", flush=True)
        if check:
            raise
        return None


def gtf_to_bed(input_gtf, output_bed):
    """
    Converts a GTF file to BED6 format, replicating the logic from gtf-to-bed.sh.
    Outputs a '.chr.bed' file.
    """
    print(f"Converting GTF {input_gtf} to BED {output_bed}...")

    # Intermediate steps will write to stdout and pipe to next command
    # Base command: convert2bed -i gtf
    # awk '$5="1000"'
    # sed 's/\ gene_name\ [^~]*"//g'
    # awk '$7=$13$15'
    # cut -f1-7 -d' '
    # sed 's/;$//g' (twice)
    # sed 's/;/./g'
    # sed 's/"//g'
    # awk '$4=$7'
    # cut -f1-6 -d' '
    # sed 's/\ /\t/g'
    # sort-bed -

    # Using shell=True for complex pipes. Ensure inputs are safe.
    # This version targets the .chr.bed output directly by modifying chromosome names early.
    # awk '$5="1000" && $1="chr"$1'
    cmd_chain = (
        f"convert2bed -i gtf < {input_gtf} | "
        f"awk 'BEGIN{{OFS=\"\\t\"}} {{ if (!/^chr/.test($1)) $1=\"chr\"$1; $5=1000; print $0 }}' | " # Add chr if missing, set score
        # Attempt to replicate the complex name creation. Original: sed 's/\ gene_name\ [^~]*"//g' | awk '$7=$13$15'
        # This part is tricky. The original script seems to rely on specific GTF attribute order or presence.
        # A robust Python GTF parser would be better here, but for direct refactoring:
        # Let's assume transcript_id and exon_number are key for the unique name.
        # A simplified name for now: transcript_id.exon_number or just transcript_id if exon_number is not standard.
        # The original script creates a unique name in field 7 ($7=$13$15 after some sed), then copies to field 4 ($4=$7)
        # For now, let's parse attributes to get transcript_id for the name field.
        # This is a simplified version of the original's complex awk/sed chain for name generation.
        f"awk 'BEGIN{{OFS=\"\\t\"}} {{"
        f"    attr_str = $9; gsub(/; /, \";\", attr_str);" # Consolidate attribute separator
        f"    split(attr_str, attrs, \";\"); transcript_id=\"\"; exon_num=\"\";"
        f"    for (i in attrs) {{ "
        f"        if (attrs[i] ~ /transcript_id\\s*\".*\"/) {{ gsub(/transcript_id\\s*\"|\"/, \"\", attrs[i]); transcript_id=attrs[i]; }} "
        f"        if (attrs[i] ~ /exon_number\\s*\".*\"/) {{ gsub(/exon_number\\s*\"|\"/, \"\", attrs[i]); exon_num=attrs[i]; }} "
        f"    }} "
        f"    if (exon_num != \"\") {{ name = transcript_id \".\" exon_num; }} else {{ name = transcript_id; }} "
        f"    $4 = (name != \"\" ? name : \"NoName\"); print $1, $2, $3, $4, $5, $6;" # Output BED6
        f"}}' | "
        f"sort-bed - > {output_bed}"
    )

    process = run_command(cmd_chain, shell=True, check=False) # sort-bed can return non-zero on empty input
    if process and process.returncode != 0 :
        # sort-bed (and other BEDOPS tools) can exit with specific codes for empty inputs, etc.
        # Check if output file is empty or if it's a real error
        if os.path.exists(output_bed) and os.path.getsize(output_bed) == 0:
            print(f"Warning: {output_bed} is empty. This might be expected if input GTF was empty or had no matching entries.")
        elif not os.path.exists(output_bed) : # No output created means likely error
            print(f"gtf_to_bed conversion failed or produced no output. Exit code: {process.returncode}")
            return False

    print(f"GTF to BED conversion complete: {output_bed}")
    return True

def extract_single_feature(track_name, track_path, input_bed, output_dir):
    """
    Extracts features using bigWigAverageOverBed or bigBedAverageOverBed.
    Simplified: Does not implement per-line summaries or special CAGE TSS logic.
    """
    os.makedirs(output_dir, exist_ok=True)
    output_tsv = os.path.join(output_dir, f"{track_name}.tsv")

    print(f"Extracting features for {track_name} from {track_path} using {input_bed}...")

    if not os.path.exists(input_bed) or os.path.getsize(input_bed) == 0:
        print(f"Input BED file {input_bed} is missing or empty. Skipping feature extraction for {track_name}.")
        # Create an empty output file to satisfy Snakemake
        with open(output_tsv, 'w') as f_out:
            # Potentially write a header if a standard empty format is expected
            # For bigWigAverageOverBed, it's usually: name, size, covered, sum, mean0, mean
            # For bigBedAverageOverBed, it's: name, size, covered, sum, mean0, mean
             f_out.write("name\tsize\tcovered\tsum\tmean0\tmean\n")
        print(f"Created empty output file: {output_tsv}")
        return True

    cmd = []
    track_path_lower = track_path.lower()

    if track_path_lower.endswith((".bw", ".bigwig")):
        cmd = [
            "bigWigAverageOverBed",
            track_path,
            input_bed,
            output_tsv,
            "-minMax" # As in the original get_avg_over_bed for .bw
        ]
    elif track_path_lower.endswith((".bb", ".bigbed")):
        # Original script used bigBedSummary per line.
        # Using bigBedAverageOverBed for consistency with bigWigAverageOverBed.
        # This will produce a single file with averages, not per-line summaries.
        cmd = [
            "bigBedAverageOverBed", # This is a change from bigBedSummary
            track_path,
            input_bed,
            output_tsv
        ]
    else:
        print(f"Unsupported file type for track {track_name}: {track_path}. Skipping.")
        # Create an empty file to satisfy Snakemake if this rule is targeted directly
        with open(output_tsv, 'w') as f_out:
             f_out.write("name\tsize\tcovered\tsum\tmean0\tmean\n") # Dummy header
        return True # Allow pipeline to continue for other features

    if not run_command(cmd):
        print(f"Feature extraction failed for {track_name}.")
        # If the command fails, it might leave an empty/incomplete tsv.
        # Snakemake might require the file to exist.
        if not os.path.exists(output_tsv):
            with open(output_tsv, 'w') as f_out: # Create empty if it failed and didn't create one
                 f_out.write("name\tsize\tcovered\tsum\tmean0\tmean\n")
        return False

    print(f"Feature extraction complete for {track_name}: {output_tsv}")
    return True

def create_paths_csv(feature_dir, output_csv):
    """
    Creates a CSV file listing feature names and paths to their TSV files.
    """
    print(f"Creating paths CSV {output_csv} for features in {feature_dir}...")
    features = []
    if os.path.isdir(feature_dir):
        for f_name in os.listdir(feature_dir):
            if f_name.endswith(".tsv"):
                feature_name = os.path.splitext(f_name)[0]
                full_path = os.path.abspath(os.path.join(feature_dir, f_name))
                features.append({"feature_name": feature_name, "path": full_path})

    if not features:
        print(f"Warning: No .tsv files found in {feature_dir}. Output CSV will be empty or have only headers.")

    df = pd.DataFrame(features)
    try:
        df.to_csv(output_csv, index=False, columns=["feature_name", "path"]) # Specify column order
    except Exception as e:
        print(f"Error writing CSV {output_csv}: {e}")
        # Create empty CSV with headers if write fails, to satisfy Snakemake
        with open(output_csv, 'w') as f_out:
            f_out.write("feature_name,path\n")
        return False

    print(f"Paths CSV created: {output_csv}")
    return True

def main():
    parser = argparse.ArgumentParser(description="Feature extraction and GTF to BED conversion.")
    subparsers = parser.add_subparsers(dest="action", required=True)

    # Subparser for gtf_to_bed
    parser_gtf = subparsers.add_parser("gtf_to_bed", help="Convert GTF to sorted BED6 format.")
    parser_gtf.add_argument("--input_gtf", required=True, help="Input GTF file.")
    parser_gtf.add_argument("--output_bed", required=True, help="Output BED file (will be sorted, with .chr prefix if not already).")

    # Subparser for extract_single_feature
    parser_feat = subparsers.add_parser("extract_single_feature", help="Extract features for a single track over a BED file.")
    parser_feat.add_argument("--track_name", required=True, help="Name of the track (for output file naming).")
    parser_feat.add_argument("--track_path", required=True, help="Path or URL to the track file (BigWig/BigBed).")
    parser_feat.add_argument("--input_bed", required=True, help="Input BED file.")
    parser_feat.add_argument("--output_dir", required=True, help="Directory to save the feature TSV file.")

    # Subparser for create_paths_csv
    parser_csv = subparsers.add_parser("create_paths_csv", help="Create a CSV listing feature names and their TSV paths.")
    parser_csv.add_argument("--feature_dir", required=True, help="Directory containing feature TSV files.")
    parser_csv.add_argument("--output_csv", required=True, help="Output CSV file path.")

    args = parser.parse_args()
    success = False

    if args.action == "gtf_to_bed":
        success = gtf_to_bed(args.input_gtf, args.output_bed)
    elif args.action == "extract_single_feature":
        success = extract_single_feature(args.track_name, args.track_path, args.input_bed, args.output_dir)
    elif args.action == "create_paths_csv":
        success = create_paths_csv(args.feature_dir, args.output_csv)

    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()
