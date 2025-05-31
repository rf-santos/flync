import csv
import subprocess
import os
import sys

def fetch_sra_info(sra_accession):
    """Fetches SRA information for a single accession using efetch."""
    print(f"Fetching info for {sra_accession}...")
    try:
        command = [
            "efetch",
            "-db", "sra",
            "-id", sra_accession,
            "-format", "runinfo"
        ]
        # The output of efetch -format runinfo is already CSV-like with a header
        # We will capture stdout and stderr
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate(timeout=60) # Added timeout

        if process.returncode != 0:
            print(f"Error fetching SRA info for {sra_accession}: {stderr.strip()}")
            return None, None

        if not stdout:
            print(f"No output received from efetch for {sra_accession}.")
            return None, None

        # The first line is the header, subsequent lines are data
        lines = stdout.strip().splitlines()
        if len(lines) < 2: # Should have at least header and one data line
            print(f"Unexpected output format for {sra_accession}: {stdout}")
            return None, None

        header = lines[0].split(',')
        data_line = lines[1].split(',') # Assuming one run per accession for simplicity here
                                        # More complex SRA entries might yield multiple lines

        return header, data_line

    except subprocess.TimeoutExpired:
        print(f"Timeout fetching SRA info for {sra_accession}.")
        return None, None
    except Exception as e:
        print(f"An unexpected error occurred while fetching info for {sra_accession}: {e}")
        return None, None

def main(sra_accession_file, output_csv_path):
    """
    Reads SRA accessions from a file, fetches their info,
    and writes it to a CSV file.
    """
    sra_ids = []
    try:
        with open(sra_accession_file, 'r') as f:
            for line in f:
                sra_ids.append(line.strip())
    except FileNotFoundError:
        print(f"Error: SRA accession file {sra_accession_file} not found.")
        sys.exit(1)

    if not sra_ids:
        print("No SRA accessions found in the input file.")
        # Create an empty CSV with default header if no SRA IDs are provided or if needed
        # For now, we let Snakemake handle empty file creation if that's desired behavior.
        # with open(output_csv_path, 'w', newline='') as csvfile:
        #    writer = csv.writer(csvfile)
        #    writer.writerow(["Run","ReleaseDate","LoadDate",...]) # Example default header
        # print(f"Empty CSV created at {output_csv_path} with headers.")
        return


    all_data = []
    processed_header = None

    for sra_id in sra_ids:
        if not sra_id: continue # Skip empty lines
        header, data_line = fetch_sra_info(sra_id)
        if data_line:
            if processed_header is None: # Store header from the first successful fetch
                processed_header = header
            all_data.append(data_line)
        else:
            print(f"Skipping {sra_id} due to previous errors.")

    if not all_data:
        print("No data fetched for any SRA accessions. Output CSV will not be created or will be empty.")
        # Optionally create an empty file or a file with only headers
        # For now, if all_data is empty, the script won't write anything,
        # and Snakemake might complain if the output file is expected but not created.
        # To ensure an empty file with headers is created:
        if processed_header is None: # If no successful fetch, use a default header
            # This default header might not match actual efetch output if all fetches fail.
            # It's a fallback. A better approach might be to error out if no data can be fetched.
            # For now, we'll let it be, or Snakemake creates an empty file.
            # Consider adding a specific efetch call for a known valid SRA ID to get a header if needed.
            print("No header could be determined as all SRA fetches failed or returned no data.")
            print(f"The output file {output_csv_path} might not be created or will be empty.")
            # To create an empty file:
            # open(output_csv_path, 'w').close()
            return


    try:
        file_exists = os.path.exists(output_csv_path)
        # We will overwrite the file each time with the new set of SRA info.
        # Snakemake handles timestamp checking for re-runs.
        with open(output_csv_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            if processed_header:
                writer.writerow(processed_header) # Write header
            writer.writerows(all_data) # Write all data rows
        print(f"SRA info successfully written to {output_csv_path}")
    except IOError as e:
        print(f"Error writing CSV file {output_csv_path}: {e}")
        sys.exit(1)


if __name__ == "__main__":
    try:
        sra_accession_file_param = snakemake.input.sra_list
        output_csv_param = snakemake.output.sra_csv
        main(sra_accession_file_param, output_csv_param)
    except NameError:
        print("Running in standalone mode (not via Snakemake). Using placeholder parameters.")
        # Example for standalone testing:
        # Create a dummy sra_list.txt with a few SRR accessions
        # with open("sra_list_test.txt", "w") as f:
        #     f.write("SRR123456\n")
        #     f.write("ERR654321\n")
        # main("sra_list_test.txt", "results/runinfo_test.csv")
        pass
