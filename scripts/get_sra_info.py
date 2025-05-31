import argparse
import os
import subprocess

def get_sra_info(workdir, sra_accession):
    """
    Fetches SRA information using esearch and efetch.
    Saves selected fields to a .info file.
    """
    output_dir = os.path.join(workdir, "sra_info")
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{sra_accession}.info")

    print(f"Fetching SRA info for {sra_accession}...")

    try:
        # Construct the efetch command to get runinfo, then select specific columns
        # The original script used: cut -f1,7,13,14,15,16,22 -d,
        # These fields correspond to: Run,AvgLength,Experiment,LibraryStrategy,LibrarySelection,LibraryLayout,BioProject
        # We will fetch the full runinfo and then process it, or directly ask for CSV.

        efetch_command = [
            "efetch",
            "-db", "sra",
            "-format", "runinfo", # Fetches a CSV format
            "-id", sra_accession
        ]

        print(f"Running command: {' '.join(efetch_command)}")
        result = subprocess.run(efetch_command, capture_output=True, text=True, check=True)

        # The output of efetch -format runinfo is a CSV.
        # We need to parse it and select the desired columns.
        # Header is usually the first line.
        lines = result.stdout.strip().split('\n')
        if len(lines) < 2:
            print(f"Error: No data returned by efetch for {sra_accession}")
            return False

        header = lines[0].split(',')
        data_line = lines[1].split(',') # Assuming one data line per accession

        # Define the columns we want based on the original script's cut command
        # Run,AvgSpotLen,Experiment,LibraryStrategy,LibrarySelection,LibraryLayout,BioProject (example header names)
        # Original cut fields: 1,7,13,14,15,16,22
        # We need to map these to actual column names in the runinfo CSV if possible,
        # or rely on their positions if the format is stable.
        # For simplicity, let's try to get them by typical names or fallback to indices.

        # A more robust way would be to find indices by header name:
        try:
            run_idx = header.index("Run")
            avglen_idx = header.index("AvgSpotLen") # Or "AvgLength"
            experiment_idx = header.index("Experiment")
            strategy_idx = header.index("LibraryStrategy")
            selection_idx = header.index("LibrarySelection")
            layout_idx = header.index("LibraryLayout")
            bioproject_idx = header.index("BioProject")

            selected_info = [
                f"Run: {data_line[run_idx]}",
                f"AvgLength: {data_line[avglen_idx]}",
                f"Experiment: {data_line[experiment_idx]}",
                f"LibraryStrategy: {data_line[strategy_idx]}",
                f"LibrarySelection: {data_line[selection_idx]}",
                f"LibraryLayout: {data_line[layout_idx]}",
                f"BioProject: {data_line[bioproject_idx]}"
            ]
        except ValueError as e:
            # Fallback to original indices if headers are not as expected (less robust)
            print(f"Warning: Could not find all headers, falling back to indices for {sra_accession}. Error: {e}")
            # Indices from `cut -f1,7,13,14,15,16,22` (0-based for list access)
            indices = [0, 6, 12, 13, 14, 15, 21]
            header_names_fallback = ["Run", "AvgLength", "Experiment", "LibraryStrategy", "LibrarySelection", "LibraryLayout", "BioProject"]
            selected_info = []
            for i, header_name in zip(indices, header_names_fallback):
                if i < len(data_line):
                    selected_info.append(f"{header_name}: {data_line[i]}")
                else:
                    selected_info.append(f"{header_name}: <Not Available>")


        with open(output_file, "w") as f:
            f.write("\n".join(selected_info))
            f.write("\n") # Add a newline at the end

        print(f"SRA info saved to {output_file}")
        return True

    except subprocess.CalledProcessError as e:
        print(f"Error running E-utilities for {sra_accession}: {e}")
        print(f"Stderr: {e.stderr}")
        print(f"Stdout: {e.stdout}")
        return False
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Get SRA information.")
    parser.add_argument("--workdir", required=True, help="Main working directory.")
    parser.add_argument("--sra_accession", required=True, help="SRA accession number.")
    args = parser.parse_args()

    if not get_sra_info(args.workdir, args.sra_accession):
        exit(1)

if __name__ == "__main__":
    main()
