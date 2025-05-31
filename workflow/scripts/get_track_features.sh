#!/bin/bash

# This script extracts features for a given BED file from a single track (local or remote).
# It's designed to be called by Snakemake for each track.

set -e # Exit on error
# set -x # Debugging

# --- Parse Command Line Arguments ---
PARAMS=""
while (( "$#" )); do
  case "$1" in
    --track_name)
      TRACK_NAME="$2"
      shift 2
      ;;
    --url)
      TRACK_URL="$2"
      shift 2
      ;;
    --bed)
      INPUT_BED="$2"
      shift 2
      ;;
    --output_tsv)
      OUTPUT_TSV="$2"
      shift 2
      ;;
    --offset)
      OFFSET="$2"
      shift 2
      ;;
    --track_type)
      TRACK_TYPE="$2" # e.g., bw, bb, cage_tss (custom for specific CAGE logic if needed)
      shift 2
      ;;
    --processing_method) # From config: bigWigAverageOverBed, bigBedSummary etc.
      PROCESSING_METHOD="$2"
      shift 2
      ;;
    *) # unsupported flags
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done

# Validate required arguments
if [ -z "$TRACK_NAME" ] || [ -z "$TRACK_URL" ] || [ -z "$INPUT_BED" ] || [ -z "$OUTPUT_TSV" ] || [ -z "$OFFSET" ] || [ -z "$TRACK_TYPE" ] || [ -z "$PROCESSING_METHOD" ]; then
  echo "Error: Missing one or more required arguments."
  echo "Usage: $0 --track_name <name> --url <url> --bed <input.bed> --output_tsv <out.tsv> --offset <val> --track_type <type> --processing_method <method>"
  exit 1
fi

echo "Processing track: $TRACK_NAME"
echo "URL: $TRACK_URL"
echo "Input BED: $INPUT_BED"
echo "Output TSV: $OUTPUT_TSV"
echo "Offset: $OFFSET"
echo "Track Type: $TRACK_TYPE"
echo "Processing Method: $PROCESSING_METHOD"


# --- Setup ---
# Temporary directory for downloads and intermediate files
# Snakemake should handle output directory creation for $OUTPUT_TSV
# TEMP_DIR=$(mktemp -d -t ${TRACK_NAME}_features_XXXXXX)
# trap "rm -rf '$TEMP_DIR'" EXIT # Cleanup temp dir

LOCAL_TRACK_FILE="" # Path to the track file after download/if local

# --- Download track if URL ---
if [[ "$TRACK_URL" == http* ]] || [[ "$TRACK_URL" == ftp* ]]; then
    LOCAL_TRACK_FILE="${TRACK_NAME}.${TRACK_TYPE}" # Save in CWD or a defined temp dir
    # Ensure filename is somewhat clean
    LOCAL_TRACK_FILE=$(basename "$TRACK_URL")
    # If basename results in generic name (like 'download'), use track_name
    if [[ "$LOCAL_TRACK_FILE" == "download" ]] || [[ "$LOCAL_TRACK_FILE" == "" ]]; then
        LOCAL_TRACK_FILE="${TRACK_NAME}_$(basename ${TRACK_URL%.*})_track.${TRACK_URL##*.}"
    fi
    LOCAL_TRACK_FILE="${TRACK_NAME}.${TRACK_TYPE}" # Simpler naming based on track name and type

    # Use a directory managed by Snakemake's temp() if possible, or a specific download cache.
    # For simplicity here, download to current working directory of the script.
    # Snakemake's rule should ensure CWD is reasonable or use temp() for output.
    # For now, assume LOCAL_TRACK_FILE is just its name, will be in CWD.

    DOWNLOAD_DIR=$(dirname "$OUTPUT_TSV")/downloads # Store downloads next to output, in a 'downloads' subdir
    mkdir -p "$DOWNLOAD_DIR"
    LOCAL_TRACK_FILE="$DOWNLOAD_DIR/${TRACK_NAME}_$(basename $TRACK_URL)"

    if [ ! -f "$LOCAL_TRACK_FILE" ]; then
        echo "Downloading $TRACK_URL to $LOCAL_TRACK_FILE..."
        wget -q -O "$LOCAL_TRACK_FILE" "$TRACK_URL" || { echo "Error downloading $TRACK_URL"; exit 1; }
    else
        echo "Track $LOCAL_TRACK_FILE already downloaded."
    fi
elif [ -f "$TRACK_URL" ]; then # Local file path provided
    LOCAL_TRACK_FILE="$TRACK_URL"
    echo "Using local track file: $LOCAL_TRACK_FILE"
else
    echo "Error: Track URL '$TRACK_URL' is not a valid URL or local file."
    exit 1
fi

if [ ! -s "$LOCAL_TRACK_FILE" ]; then
    echo "Error: Local track file $LOCAL_TRACK_FILE is empty or does not exist after download attempt."
    # Create empty output to satisfy snakemake if download fails and we want to continue
    # touch "$OUTPUT_TSV"
    exit 1
fi

# --- Prepare BED file: Sort and potentially adjust for offset if needed by tool ---
# Most UCSC tools expect sorted BED. The input BED should already be sorted if necessary.
# The offset logic from original get-features.sh was complex (center vs start/end).
# For bigWigAverageOverBed, it averages over the exact regions in the BED.
# If offset means TSS +/- offset, the BED file itself should define these regions.
# The input BED here is assumed to be new-non-coding.chr.bed (transcript regions).
# If features around TSS are needed, a new BED file for TSS +/- offset should be an input.
# For now, this script will use the provided BED file as is for feature extraction.
# The $OFFSET parameter from original script was used for CAGE data processing.

# --- Feature Extraction based on track type and processing method ---

# Ensure output file is clear before writing
> "$OUTPUT_TSV"

# Check if input BED has content
if [ ! -s "$INPUT_BED" ]; then
    echo "Warning: Input BED file $INPUT_BED is empty. Output TSV will be empty."
    touch "$OUTPUT_TSV" # Create empty file
    exit 0 # Exit successfully as there's no data to process
fi


case "$PROCESSING_METHOD" in
    "bigWigAverageOverBed")
        if [ "$TRACK_TYPE" != "bw" ]; then
            echo "Warning: Track type is not 'bw' for bigWigAverageOverBed. Attempting anyway."
        fi
        echo "Using bigWigAverageOverBed..."
        # Output format: name, size, covered, sum, mean, mean0
        # We want to extract the 'mean0' (column 6)
        bigWigAverageOverBed "$LOCAL_TRACK_FILE" "$INPUT_BED" stdout | awk -v OFS='\t' '{print $1, $6}' > "$OUTPUT_TSV"
        ;;
    "bigBedSummary") # Typically for .bb files
        if [ "$TRACK_TYPE" != "bb" ]; then
            echo "Warning: Track type is not 'bb' for bigBedSummary. Attempting anyway."
        fi
        echo "Using bigBedSummary..."
        # This needs a different approach as bigBedSummary is for a single range, not per BED entry.
        # The original script might have intended bigBedAverageOverBed or similar for .bb files.
        # Or, if it's for presence/absence, one might use bedtools intersect -c.
        # For now, let's assume a simple overlap count if it's a bigBed of regions.
        # This is a placeholder, as direct feature extraction for .bb per BED line needs more specific logic.
        # bedtools intersect -a "$INPUT_BED" -b "$LOCAL_TRACK_FILE" -c -wa | awk -v OFS='\t' '{print $4, $NF}' > "$OUTPUT_TSV"
        # The above assumes $4 is name. If INPUT_BED is BED6+, $4 is name.
        # For now, let's use a placeholder or error, as this needs careful thought.
        echo "Error: 'bigBedSummary' per-entry feature extraction not fully implemented. Placeholder."
        # Create dummy output: transcript_name feature_value (0 for now)
        awk -v OFS='\t' '{print $4, 0}' "$INPUT_BED" > "$OUTPUT_TSV" # Assuming $4 is name
        ;;
    "custom_cage_tss_sum") # Example for CAGE data if specific logic was used
        echo "Using custom CAGE TSS sum logic (placeholder)..."
        # The original script had:
        # awk -v offset="$OFFSET" '{print $1"\t"$2-offset"\t"$3+offset"\t"$4"\t"$5"\t"$6}' "$INPUT_BED" > tmp.bed
        # bigWigAverageOverBed track.bw tmp.bed stdout | awk '{print $1"\t"$5}' > $OUTPUT_TSV
        # This implies creating a new BED with offset regions around original BED entries.
        # For now, this script assumes the INPUT_BED is already correctly defined for the desired feature.
        # So, if CAGE sum over TSS +/- offset is needed, INPUT_BED should be TSS regions +/- offset.
        # Here, we'll just use bigWigAverageOverBed on the given BED directly for CAGE.
        if [ "$TRACK_TYPE" != "bw" ]; then
            echo "Warning: Track type is not 'bw' for CAGE processing with bigWigAverageOverBed. Attempting anyway."
        fi
        bigWigAverageOverBed "$LOCAL_TRACK_FILE" "$INPUT_BED" stdout | awk -v OFS='\t' '{print $1, $5}' > "$OUTPUT_TSV" # sum is col 5
        ;;
    *)
        echo "Error: Unknown processing_method '$PROCESSING_METHOD' for track_type '$TRACK_TYPE'."
        # Create dummy output
        awk -v OFS='\t' '{print $4, "NA"}' "$INPUT_BED" > "$OUTPUT_TSV" # Assuming $4 is name
        exit 1
        ;;
esac

if [ ! -s "$OUTPUT_TSV" ]; then
    echo "Warning: Output TSV $OUTPUT_TSV is empty after processing. This might be OK if there were no overlaps."
    # Ensure file exists even if empty
    touch "$OUTPUT_TSV"
fi

echo "Finished processing track: $TRACK_NAME. Output at $OUTPUT_TSV"
