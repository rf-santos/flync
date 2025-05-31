import argparse
import os
import requests
import gzip
import shutil
import subprocess

def download_file(url, output_path):
    """Downloads a file from a URL to a given path."""
    print(f"Downloading {url} to {output_path}...")
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()  # Raise an exception for bad status codes
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print("Download complete.")
    except requests.exceptions.RequestException as e:
        print(f"Error downloading {url}: {e}")
        return False
    return True

def decompress_gz(gz_path, output_path):
    """Decompresses a .gz file."""
    print(f"Decompressing {gz_path} to {output_path}...")
    try:
        with gzip.open(gz_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(gz_path)  # Remove the .gz file after decompression
        print("Decompression complete.")
    except Exception as e:
        print(f"Error decompressing {gz_path}: {e}")
        return False
    return True

def main():
    parser = argparse.ArgumentParser(description="Download and prepare genome files.")
    parser.add_argument("--output_dir", required=True, help="Directory to save genome files.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Read URLs from static/required_links.txt
    required_links_path = os.path.join(os.path.dirname(__file__), "../static/required_links.txt")
    if not os.path.exists(required_links_path):
        print(f"Error: {required_links_path} not found.")
        return

    genome_fasta_urls = []
    genome_gtf_url = None

    with open(required_links_path, "r") as f:
        for line in f:
            line = line.strip()
            if "Drosophila_melanogaster.BDGP6.32.dna.primary_assembly" in line and line.endswith(".fa.gz"):
                genome_fasta_urls.append(line)
            elif "Drosophila_melanogaster.BDGP6.32.106.chr.gtf.gz" in line: # Specific GTF link from original script
                genome_gtf_url = line

    if not genome_fasta_urls:
        print("Error: No genome FASTA URLs found in required_links.txt.")
        return
    if not genome_gtf_url:
        print("Error: Genome GTF URL not found in required_links.txt.")
        return

    # Download and process genome FASTA files
    downloaded_fasta_parts = []
    for url in genome_fasta_urls:
        filename_gz = os.path.join(args.output_dir, os.path.basename(url))
        filename_fa = filename_gz[:-3]  # Remove .gz
        if download_file(url, filename_gz):
            if decompress_gz(filename_gz, filename_fa):
                downloaded_fasta_parts.append(filename_fa)
            else:
                print(f"Failed to decompress {filename_gz}")
                return # Stop if decompression fails
        else:
            print(f"Failed to download {url}")
            return # Stop if download fails

    # Concatenate FASTA parts
    final_fasta_path = os.path.join(args.output_dir, "genome.fa")
    print(f"Concatenating FASTA parts to {final_fasta_path}...")
    with open(final_fasta_path, 'wb') as outfile:
        for fa_part_path in downloaded_fasta_parts:
            with open(fa_part_path, 'rb') as infile:
                shutil.copyfileobj(infile, outfile)
            os.remove(fa_part_path) # Clean up individual parts
    print("FASTA concatenation complete.")

    # Download and process genome GTF file
    gtf_filename_gz = os.path.join(args.output_dir, os.path.basename(genome_gtf_url))
    final_gtf_path = os.path.join(args.output_dir, "genome.gtf")
    if download_file(genome_gtf_url, gtf_filename_gz):
        if not decompress_gz(gtf_filename_gz, final_gtf_path):
            print(f"Failed to decompress {gtf_filename_gz}")
            return
    else:
        print(f"Failed to download {genome_gtf_url}")
        return

    print(f"Genome preparation complete. Files are in {args.output_dir}")

if __name__ == "__main__":
    main()
