import os
import sys
import gzip
import shutil
import requests # Using requests for downloading

def download_file(url, target_path, chunk_size=8192):
    """Downloads a file from a URL to a target path with progress."""
    try:
        response = requests.get(url, stream=True, timeout=60) # Added timeout
        response.raise_for_status()  # Raise an exception for bad status codes
        file_size = int(response.headers.get('content-length', 0))

        if os.path.exists(target_path) and os.path.getsize(target_path) == file_size and file_size != 0 :
            print(f"File {target_path} already exists and is complete. Skipping download.")
            return True

        print(f"Downloading {url} to {target_path}...")
        with open(target_path, 'wb') as f:
            downloaded_size = 0
            for chunk in response.iter_content(chunk_size=chunk_size):
                f.write(chunk)
                downloaded_size += len(chunk)
                progress = int(50 * downloaded_size / file_size) if file_size else 0
                sys.stdout.write(f"\r[{'=' * progress}{' ' * (50 - progress)}] {downloaded_size}/{file_size} bytes")
                sys.stdout.flush()
        sys.stdout.write("\nDownload complete.\n")
        return True
    except requests.exceptions.RequestException as e:
        print(f"Error downloading {url}: {e}")
        if os.path.exists(target_path): # Clean up partial download
            os.remove(target_path)
        return False
    except IOError as e:
        print(f"Error writing file {target_path}: {e}")
        if os.path.exists(target_path): # Clean up partial download
            os.remove(target_path)
        return False


def decompress_gz(gz_path, output_path):
    """Decompresses a .gz file."""
    if os.path.exists(output_path):
        print(f"File {output_path} already exists. Skipping decompression.")
        return True

    print(f"Decompressing {gz_path} to {output_path}...")
    try:
        with gzip.open(gz_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print("Decompression complete.")
        return True
    except (gzip.BadGzipFile, IOError) as e:
        print(f"Error decompressing {gz_path}: {e}")
        if os.path.exists(output_path): # Clean up partial file
             os.remove(output_path)
        return False

def main(genome_dir, ensembl_base_url, ensembl_gtf_url, genome_fasta_name, chromosomes, output_fasta_path, output_gtf_path):
    os.makedirs(genome_dir, exist_ok=True)

    # Download and concatenate chromosome FASTA files
    concatenated_fasta_path_tmp = os.path.join(genome_dir, f"{genome_fasta_name}.all_chromo.fa.tmp")
    final_fasta_path = output_fasta_path

    if os.path.exists(final_fasta_path):
        print(f"Final genome FASTA {final_fasta_path} already exists. Skipping chromosome download and concatenation.")
    else:
        print("Processing chromosome FASTA files...")
        with open(concatenated_fasta_path_tmp, 'wb') as final_fasta_file:
            for chrom in chromosomes:
                fasta_file_name = f"{genome_fasta_name}.chromosome.{chrom}.fa.gz"
                fasta_url = f"{ensembl_base_url}{fasta_file_name}"
                local_fasta_gz_path = os.path.join(genome_dir, fasta_file_name)
                local_fasta_path = os.path.join(genome_dir, f"{genome_fasta_name}.chromosome.{chrom}.fa")

                if not download_file(fasta_url, local_fasta_gz_path):
                    print(f"Failed to download chromosome {chrom}. Aborting genome preparation.")
                    if os.path.exists(concatenated_fasta_path_tmp): os.remove(concatenated_fasta_path_tmp)
                    return False

                if not decompress_gz(local_fasta_gz_path, local_fasta_path):
                    print(f"Failed to decompress chromosome {chrom}. Aborting genome preparation.")
                    if os.path.exists(concatenated_fasta_path_tmp): os.remove(concatenated_fasta_path_tmp)
                    return False

                with open(local_fasta_path, 'rb') as chrom_fasta:
                    shutil.copyfileobj(chrom_fasta, final_fasta_file)
                print(f"Appended {chrom} to {concatenated_fasta_path_tmp}")
                os.remove(local_fasta_path) # Clean up individual chromosome file
                os.remove(local_fasta_gz_path) # Clean up downloaded gz file


        shutil.move(concatenated_fasta_path_tmp, final_fasta_path)
        print(f"Concatenated genome FASTA created at {final_fasta_path}")

    # Download and decompress GTF file
    gtf_gz_file_name = os.path.basename(ensembl_gtf_url) # e.g., Drosophila_melanogaster.BDGP6.32.106.chr.gtf.gz
    local_gtf_gz_path = os.path.join(genome_dir, gtf_gz_file_name)
    final_gtf_path = output_gtf_path # e.g. genome_dir/genome.gtf

    if os.path.exists(final_gtf_path):
        print(f"Final GTF file {final_gtf_path} already exists. Skipping download and decompression.")
    else:
        if not download_file(ensembl_gtf_url, local_gtf_gz_path):
            print("Failed to download GTF file. Aborting.")
            return False
        if not decompress_gz(local_gtf_gz_path, final_gtf_path):
            print("Failed to decompress GTF file. Aborting.")
            if os.path.exists(local_gtf_gz_path): os.remove(local_gtf_gz_path) # Clean up .gz if decompression fails
            return False
        print(f"GTF file processed and saved to {final_gtf_path}")
        # os.remove(local_gtf_gz_path) # Clean up downloaded .gz file after successful decompression

    print("Genome preparation script finished.")
    return True

if __name__ == "__main__":
    # This part is for Snakemake integration
    # Snakemake will provide 'snakemake' object
    try:
        genome_dir_param = snakemake.params.genome_dir
        ensembl_base_url_param = snakemake.params.ensembl_base_url
        ensembl_gtf_url_param = snakemake.params.ensembl_gtf_url
        genome_fasta_name_param = snakemake.params.genome_fasta_name
        chromosomes_param = snakemake.params.chromosomes

        output_fasta_param = snakemake.output.genome_fasta
        output_gtf_param = snakemake.output.genome_gtf

        if not main(genome_dir_param, ensembl_base_url_param, ensembl_gtf_url_param, genome_fasta_name_param, chromosomes_param, output_fasta_param, output_gtf_param):
            sys.exit(1) # Indicate failure to Snakemake

    except NameError:
        # Fallback for testing outside Snakemake (optional)
        print("Running in standalone mode (not via Snakemake). Using placeholder parameters.")
        # Define some defaults for testing if you want to run it directly
        # For example:
        # main(
        #     genome_dir="workflow/genome_data_test",
        #     ensembl_base_url="https://ftp.ensembl.org/pub/release-106/fasta/drosophila_melanogaster/dna/",
        #     ensembl_gtf_url="https://ftp.ensembl.org/pub/release-106/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.106.chr.gtf.gz",
        #     genome_fasta_name="Drosophila_melanogaster.BDGP6.32.dna.primary_assembly",
        #     chromosomes=["4", "mitochondrion_genome"], # Smaller set for testing
        #     output_fasta_path="workflow/genome_data_test/genome.fa",
        #     output_gtf_path="workflow/genome_data_test/genome.gtf"
        # )
        pass
