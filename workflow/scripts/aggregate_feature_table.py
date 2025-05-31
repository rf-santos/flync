#!/usr/bin/env python

# Original Author: Prof. Folkers
# Modified for Snakemake integration

import pandas as pd
import os
import sys

def load_bed_data(bed_file_path):
    """
    Loads relevant data from a BED file.
    Extracts transcript name (column 4) and calculates length (end - start).
    Assumes BED-like format where column 0 is chrom, 1 is start, 2 is end, 3 is name.
    """
    try:
        bed_df = pd.read_csv(
            bed_file_path,
            sep='\t',
            header=None,
            usecols=[0, 1, 2, 3],
            names=['chrom', 'start', 'end', 'transcript_id']
        )
        bed_df['length'] = bed_df['end'] - bed_df['start']
        # Keep only transcript_id and length for the feature table
        return bed_df[['transcript_id', 'length']].set_index('transcript_id')
    except FileNotFoundError:
        sys.exit(f"Error: BED file not found at {bed_file_path}")
    except Exception as e:
        sys.exit(f"Error processing BED file {bed_file_path}: {e}")


def load_feature_data(feature_file_path):
    """
    Loads data from a feature TSV file.
    Assumes a two-column TSV: transcript_id, feature_value.
    The feature name is derived from the file name (e.g., CAGE_pos.feature.tsv -> CAGE_pos).
    """
    try:
        feature_df = pd.read_csv(feature_file_path, sep='\t', header=None, names=['transcript_id', 'value'])
        feature_name = os.path.basename(feature_file_path).replace(".feature.tsv", "")
        feature_df.rename(columns={'value': feature_name}, inplace=True)
        return feature_df.set_index('transcript_id')
    except FileNotFoundError:
        print(f"Warning: Feature file not found: {feature_file_path}. Skipping.")
        return None
    except pd.errors.EmptyDataError:
        print(f"Warning: Feature file is empty: {feature_file_path}. Skipping.")
        return None
    except Exception as e:
        print(f"Warning: Error processing feature file {feature_file_path}: {e}. Skipping.")
        return None

def main(input_feature_tsvs, input_bed_file, output_csv_path):
    """
    Main function to aggregate features.
    """
    print(f"Loading base transcript data from BED file: {input_bed_file}")
    final_df = load_bed_data(input_bed_file)

    if final_df.empty:
        print("BED file was empty or contained no transcripts. Resulting feature table will be empty.")
        # Create an empty DataFrame with expected columns if final_df is empty
        # This depends on what feature names would have been.
        # For now, just write an empty CSV.
        pd.DataFrame().to_csv(output_csv_path)
        print(f"Empty feature table (or header only) written to {output_csv_path}")
        return

    print(f"Found {len(final_df)} transcripts in BED file.")

    num_features_loaded = 0
    for feature_file in input_feature_tsvs:
        print(f"Processing feature file: {feature_file}")
        feature_df = load_feature_data(feature_file)
        if feature_df is not None:
            if not feature_df.empty:
                # Use outer join to keep all transcripts from the base BED
                # and all transcripts from the feature file.
                # Transcripts not in the feature file will get NaN for that feature.
                # Transcripts in feature file but not BED (should not happen if BED is primary list)
                # would also be kept with NaN for length.
                final_df = final_df.join(feature_df, how='left')
                num_features_loaded +=1
            else:
                 print(f"Feature file {feature_file} was empty after loading. Skipping join.")
        else:
            print(f"Skipped feature file {feature_file} due to loading error/file not found.")

    if num_features_loaded == 0 :
        print("Warning: No feature files were successfully loaded or they were all empty. The output table will only contain transcript IDs and lengths.")

    # Ensure the index is named 'transcript_id' for the final CSV
    final_df.index.name = 'transcript_id'

    # Fill NaN values with 0 or a specific placeholder if appropriate,
    # for example, if a feature being NaN means it was 0 (e.g., no CAGE signal)
    # final_df.fillna(0, inplace=True) # Example: fill NaNs with 0 - use with caution

    print(f"Writing final feature table to: {output_csv_path}")
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
    final_df.to_csv(output_csv_path, index=True) # index=True to write the transcript_id

    print("Feature table aggregation finished.")

if __name__ == "__main__":
    try:
        # Snakemake integration
        main(input_feature_tsvs=snakemake.input.feature_tsvs,
             input_bed_file=snakemake.input.bed_file,
             output_csv_path=snakemake.output.final_table
            )
    except NameError:
        # Fallback for testing
        print("Running in standalone mode for testing. Requires dummy files.")
        # Create dummy files for testing:
        # Dummy BED: results/classification/new-non-coding.bed
        # Dummy feature TSVs: results/feature_extraction/F1.feature.tsv, results/feature_extraction/F2.feature.tsv

        # Example usage:
        # base_dir = "results"
        # os.makedirs(f"{base_dir}/classification", exist_ok=True)
        # os.makedirs(f"{base_dir}/feature_extraction", exist_ok=True)

        # dummy_bed_path = f"{base_dir}/classification/new-non-coding.bed"
        # with open(dummy_bed_path, "w") as f:
        #     f.write("chr1\t10\t100\tT1\n")
        #     f.write("chr1\t200\t300\tT2\n")

        # dummy_feat1_path = f"{base_dir}/feature_extraction/FeatA.feature.tsv"
        # with open(dummy_feat1_path, "w") as f:
        #     f.write("T1\t0.5\n")
        #     f.write("T2\t0.8\n")

        # dummy_feat2_path = f"{base_dir}/feature_extraction/FeatB.feature.tsv"
        # with open(dummy_feat2_path, "w") as f:
        #     f.write("T1\t100\n")
        #     # T2 missing for FeatB - will result in NaN, then fillna if used.

        # main(input_feature_tsvs=[dummy_feat1_path, dummy_feat2_path],
        #      input_bed_file=dummy_bed_path,
        #      output_csv_path=f"{base_dir}/ml_feature_table_test.csv")
        pass
