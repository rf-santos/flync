import pandas as pd
import pybedtools
import sys

def get_features(bed_file, tracks_file, output_dir):
    """
    Extracts features from bigWig/bigBed files for a given set of genomic intervals in a BED file.
    """
    bed = pybedtools.BedTool(bed_file)
    tracks = pd.read_csv(tracks_file, sep='\t', header=None, names=['track', 'url'])

    for _, row in tracks.iterrows():
        track_name = row['track']
        track_url = row['url']
        output_file = f"{output_dir}/{track_name}.tsv"

        if track_url.endswith('.bw'):
            # Use map function of pybedtools to get the average value
            # The score column (column 5) will have the average value
            bed.map(b=track_url, c=5, o='mean').saveas(output_file)
        elif track_url.endswith('.bb'):
            # For bigBed files, we can use the coverage tool
            bed.coverage(b=track_url).saveas(output_file)

if __name__ == '__main__':
    bed_file = sys.argv[1]
    tracks_file = sys.argv[2]
    output_dir = sys.argv[3]
    get_features(bed_file, tracks_file, output_dir)
