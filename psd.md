Product Specification Document: FLYNC v2

Version: 1.0

Date: October 26, 2023

Author: Bard (AI Assistant)

1. Introduction

1.1. Project Overview

FLYNNC (FLY Long Non-Coding RNA Classifier) is a bioinformatics pipeline designed to identify and classify novel long non-coding RNAs (lncRNAs) from transcriptomic data. FLYNC v1, implemented primarily in shell scripts, leverages the Tuxedo2 protocol and CPAT for bioinformatics analysis, and a pre-trained machine learning (ML) model for final lncRNA classification. FLYNC integrates bioinformatics and artificial intelligence stages to achieve this. FLYNC v1 needs modernization to improve its robustness, performance, usability, and maintainability. FLYNC v2 aims to address these needs by adopting a Python-first approach and incorporating several enhancements.

1.2. Purpose of this Document

This Product Specification Document (PSD) outlines the requirements and specifications for FLYNC v2, a modernized and improved version of the application. The primary goal of FLYNC v2 is to enhance usability, performance, maintainability, and extensibility by adopting a Python-first approach, improving the underlying architecture, and incorporating best practices in software development. This document serves as a guide for developers to implement FLYNC v2. It provides a detailed roadmap for developers to create a significantly improved application ready for peer review and wider adoption.

2. Goals

The primary goals for FLYNC v2 are to:

**Modernize the codebase:**
    * Refactor the existing shell script-based pipeline into a Python-based application for improved maintainability, readability, and extensibility.
    * Refactor shell scripts into Python scripts for improved maintainability and extensibility.

**Enhance User Experience:**
    * Improve user feedback, logging, configuration, and output formats for better usability and interpretability.
    * Provide comprehensive application logging for debugging and monitoring.
    * Offer detailed console feedback to users during execution.
    * Create user-friendly, human-readable output files and formats.
    * Enable easier configuration and job submission through an abstract configuration file.

**Improve Performance and Efficiency:**
    * Optimize performance through parallelism, efficient data handling, and faster data download strategies.
    * Implement parallelism and multithreading to accelerate computationally intensive steps.
    * Optimize data download speeds from external resources.
    * Ensure reliable access to external data sources.

**Increase Robustness and Reliability:**
    * Implement checkpointing, error handling, and ensure external resource availability for a more robust pipeline.
    * Incorporate checkpoint logic to handle potential step failures and allow for restarts.
    * Thoroughly test the application to ensure smooth operation.

**Increase Flexibility and Customization:**
    * Facilitate Parameter Tuning and Customization: Provide users with more control over key parameters in the bioinformatics pipeline.
    * Provide mechanisms for users to tweak key parameters of read mapping and transcript assembly.
    * Improve the management and modularity of the bioinformatics pipeline.

**Prepare for Future Expansion:**
    * Design a modular and extensible architecture to accommodate future feature additions, including model retraining and integration of new data sources.
    * Facilitate model retraining with new features and datasets.

3. Target Audience

The target audience for FLYNC v2 remains the same as v1:

* Bioinformaticians and Computational Biologists: Researchers who need a robust and user-friendly pipeline for lncRNA discovery and classification from RNA-Seq data.
* Researchers in genomics and transcriptomics interested in lncRNA discovery.
* Molecular Biologists and Geneticists: Researchers interested in identifying and characterizing novel lncRNAs in their biological systems.
* Data scientists and machine learning practitioners working in the field of non-coding RNA biology.
* Researchers with limited command-line experience: While command-line usage will still be necessary, FLYNC v2 should strive for improved clarity and ease of use.

4. Functional Requirements

FLYNC v2 will maintain the core two-stage architecture (Bioinformatics and AI) but with significant improvements in each stage.

4.1. Core Pipeline Refactoring (Python-first Approach)

This stage will be implemented entirely in Python, replacing existing shell scripts.

FR-4.1.1. Python Scripting: All shell scripts from FLYNC v1 must be refactored into Python scripts. This includes:

* Workflow orchestration and pipeline management.
* Data processing steps (mapping, assembly, abundance estimation, etc.).
* Interaction with external tools (e.g., read mappers, assemblers, CPAT).
* Feature extraction from UCSC Genome Browser.
* ML model prediction.

FR-4.1.2. Modular Design: The pipeline should be designed with a modular architecture, separating distinct steps into independent Python modules or functions. This will improve code organization, maintainability, and reusability.

FR-4.1.3. Workflow Management System Integration (Recommended): Explore integrating a workflow management system like Snakemake or Nextflow to manage the pipeline execution, dependencies, and parallelization. This is highly recommended for long-term maintainability and scalability and is considered a potential future enhancement. For v2, focus on clear Python code structure and modularity.

4.2. Improved Logging

FR-4.2.1. Comprehensive Logging: Implement a robust logging system using Python's logging module.

FR-4.2.2. Log Levels: Support different log levels (DEBUG, INFO, WARNING, ERROR, CRITICAL) to control the verbosity of logging output. Users should be able to configure the log level.

FR-4.2.3. Log Output: Log messages should include timestamps, log levels, module names, and informative descriptions of events.

FR-4.2.4. Log Destinations: Support logging to both console (stdout/stderr) and log files. Users should be able to specify the log file path. Log file should be configurable in config.yaml (path, filename).

FR-4.2.5. Error and Warning Reporting: Clearly log errors and warnings, providing context and potential solutions where possible. Log important events: pipeline start/end, step start/end, errors, warnings, resource downloads, parameter settings, etc.

4.3. Enhanced User Feedback

FR-4.3.1. Console Progress Updates: Provide informative and real-time progress updates to the user through console prints during long-running steps. This should include:

* Current step being executed.
* Percentage completion for long tasks (e.g., data download, mapping, assembly).
* Estimated time remaining (where feasible).
* Use progress bars for long-running steps (e.g., downloads, read mapping, assembly).

FR-4.3.2. Clear Status Messages: Display clear status messages indicating success, failure, or warnings for each step in the pipeline.

FR-4.3.3. User Prompts and Instructions: Provide clear prompts and instructions to the user when input is required or when choices need to be made (e.g., selector steps).

FR-4.3.4. Error Troubleshooting: Display clear and user-friendly error messages in case of failures, guiding users on how to troubleshoot.

FR-4.3.5. Summary on Completion: Summarize key steps and outputs upon successful pipeline completion.

4.4. Parallelism and Multi-threading

FR-4.4.1. Parallel Processing: Implement parallelism for computationally intensive steps like read mapping, transcript assembly, and feature extraction. Implement parallelism and multithreading for CPU-intensive tasks: Read mapping, Transcriptome assembly, Feature extraction.

FR-4.4.2. Multi-threading/Multi-processing: Utilize Python's multiprocessing or threading modules (or workflow management system capabilities) to leverage multi-core processors and improve execution speed. Utilize Python's multiprocessing or threading modules for parallel execution where appropriate.

FR-4.4.3. Configurable Parallelism: Allow users to configure the number of threads/processes to use for parallel steps, potentially through command-line arguments or the configuration file. Resource allocation (number of threads/processes, memory) should be configurable in config.yaml.

FR-4.4.4. Efficient Data Handling: Optimize data loading and processing to minimize bottlenecks and maximize the benefits of parallelism.

4.5. External Source Reachability and Data Download

FR-4.5.1. Source Availability Checks: Implement checks at the beginning of the pipeline to ensure that all required external sources (Ensembl, SRA, UCSC) are reachable and responsive. Report errors to the user if any source is unavailable. Implement checks at the beginning of the pipeline to verify the reachability of Ensembl, SRA, and UCSC Genome Browser.

FR-4.5.2. Robust Download Mechanisms: Use robust libraries like requests or urllib for downloading data from external sources, handling potential network issues and retries. Use robust HTTP request methods with timeouts and retry mechanisms to handle temporary network issues.

FR-4.5.3. Improved Download Speed: Optimize download speed from external resources.

FR-4.5.3.1. Asynchronous Downloads: Explore asynchronous download strategies (e.g., using asyncio or aiohttp) to download data from multiple sources concurrently. Utilize efficient HTTP libraries like requests with connection pooling and aiohttp for asynchronous requests if appropriate.

FR-4.5.3.2. Efficient Data Transfer: Utilize efficient data transfer protocols and methods offered by external sources (e.g., FTP, rsync, API endpoints). Consider using mirror sites or CDNs for faster downloads if available for Ensembl and UCSC.

FR-4.5.3.3. Caching (Optional but Recommended): Implement local caching of downloaded data (e.g., reference genome, annotations) to avoid redundant downloads in subsequent runs. Users should be able to clear the cache. Use parallel downloads for large files (e.g., reference genome, UCSC tracks).

FR-4.5.4. Error Reporting for Source Issues: Provide informative error messages if external sources are unreachable, suggesting potential causes (network issues, service downtime).

4.6. Testing with Provided Files

FR-4.6.1. Test Data Integration: Integrate the provided @metadata.csv and @test-list.txt files into the application's testing framework.

FR-4.6.2. Automated Testing: Develop automated tests that use these files to ensure the pipeline runs smoothly and produces expected outputs.

FR-4.6.3. Test Coverage: Aim for comprehensive test coverage of all pipeline stages and functionalities. Ensure test coverage for all critical functionalities and edge cases.

4.7. Improved Configuration (@config.yaml)

FR-4.7.1. Abstract Configuration: Redesign @config.yaml to be more abstract and user-friendly, separating configuration parameters into logical sections (e.g., input, output, genome, mapping, assembly, ml_model). YAML-based configuration file for job submission, offering an abstract way to define parameters.

FR-4.7.2. Parameter Descriptions: Include clear descriptions and default values for each parameter in the configuration file. Parameters should be clearly documented in the configuration file and in the application documentation. Provide sensible default parameter values.

FR-4.7.3. Command-line Overrides: Allow users to override configuration parameters through command-line arguments for flexibility.

FR-4.7.4. Validation: Implement configuration file validation to check for missing or invalid parameters and provide informative error messages.

FR-4.7.5. Configuration Parameters: Parameters should include:
    * Input data specification (FASTQ paths or SRA accessions).
    * Output directory path.
    * Reference genome source (Ensembl version, genome build).
    * Sample metadata file (optional, CSV format) for differential expression analysis.
    * Parameters for read mapping (e.g., aligner, parameters).
    * Parameters for transcriptome assembly (e.g., assembler, parameters).
    * Differential expression analysis settings (optional).
    * Feature extraction settings (UCSC tracks).
    * ML model path.
    * Resource allocation parameters (threads, memory).

4.8. Checkpoint Logic

FR-4.8.1. Checkpoint Implementation: Implement checkpointing at key stages of the pipeline (e.g., after read mapping, transcript assembly, feature extraction). Implement checkpointing at each major step of the pipeline (e.g., after read mapping, transcriptome assembly, feature extraction).

FR-4.8.2. Checkpoint Files: Store checkpoint data in files that can be easily loaded to resume the pipeline from the last successful checkpoint.

FR-4.8.3. Resumption Logic: Implement logic to automatically detect and load checkpoint files if available, allowing the pipeline to resume execution from the interrupted point. In case of failure during a step, allow the pipeline to be restarted from the last successful checkpoint, saving computational time and resources. Before starting a step, check if the output of that step already exists. If it does, skip the step and proceed to the next one.

FR-4.8.4. User Control: Provide users with options to:
    * Enable/disable checkpointing. Checkpointing logic should be configurable (e.g., enable/disable, checkpoint directory) in config.yaml.
    * Specify the checkpoint directory.
    * Force restart from the beginning, ignoring checkpoints.

4.9. Improved Output Files and Formats

FR-4.9.1. Human-Readable Output: Output files should be formatted for easy human consultation and analysis. Improve output file formats for human readability and downstream analysis.

FR-4.9.2. Standard Formats: Utilize standard bioinformatics file formats where applicable (e.g., BAM, GTF, CSV, TSV). Use standard bioinformatics file formats where applicable (BAM, GTF). For tabular data, use TSV (Tab-Separated Values) or CSV (Comma-Separated Values) formats with clear headers.

FR-4.9.3. Organized Output Directory: Organize output files into a structured directory hierarchy for clarity and easy navigation. Organize output files into a structured directory hierarchy based on the pipeline steps and samples.

FR-4.9.4. Summary Reports: Generate summary reports (e.g., in text or HTML format) at the end of the pipeline run, summarizing key results, statistics, and potential issues. Provide comprehensive summary files (e.g., HTML or Markdown reports) that summarize the pipeline execution, key parameters, and main results.

FR-4.9.5. Customizable Output: Consider allowing users to customize output file names and formats to some extent.

4.10. Parameter Tuning for Read Mapping and Transcript Assembly

FR-4.10.1. Expose Key Parameters: Identify and expose key parameters from the read mapping (e.g., STAR, HISAT2) and transcript assembly (e.g., StringTie, Cufflinks) tools in the configuration file. These parameters should be relevant for performance tuning and result customization. Examples include:
    * Number of threads for mapping and assembly.
    * Minimum and maximum intron size.
    * Minimum read alignment score.
    * Minimum transcript coverage.
    * Parameter Tweaking: Allow users to specify aligner and adjust key parameters (e.g., number of threads, mismatch tolerance) in the config.yaml for read mapping.
    * Parameter Tweaking: Allow users to specify assembler and adjust key parameters (e.g., minimum read coverage, minimum isoform fraction) in the config.yaml for transcriptome assembly.
    * Users should be able to tweak key parameters for read mapping and transcriptome assembly through the config.yaml file.

FR-4.10.2. Parameter Documentation: Provide clear documentation for each exposed parameter, explaining its purpose and potential impact on the results.

FR-4.10.3. Default Parameter Sets: Provide reasonable default parameter sets for common use cases.

4.11. Improved Bioinformatics Pipeline Management

FR-4.11.1. Clear Pipeline Stages: Define and clearly delineate the different stages of the bioinformatics pipeline (BI Stage) in the code and documentation.

FR-4.11.2. Pipeline Execution Control: Provide mechanisms to:
    * Run the entire pipeline from start to finish.
    * Run specific stages of the pipeline independently (for debugging or re-analysis).
    * Skip certain stages if needed.

FR-4.11.3. Dependency Management (If using Workflow System): If a workflow management system is integrated, leverage its dependency management capabilities to ensure steps are executed in the correct order and dependencies are met.

4.12. Bioinformatics (BI) Stage Modules Flow

This section details the modules within the Bioinformatics (BI) Stage and their flow.

4.12.1. Selector (Input Reads):
    * Python logic to determine the input type based on the configuration (config.yaml).
    * Branch execution flow based on whether input is local FASTQ files or SRA accessions.

4.12.2. Run Info & Map SRA Reads:
    * Python module to download SRA reads using SRA Toolkit (or equivalent Python library).
    * Extract relevant run metadata (if needed, specify what metadata is required).
    * Perform read mapping of SRA reads to the reference genome using a chosen aligner (e.g., STAR, HISAT2).
    * Parameter Tweaking: Allow users to specify aligner and adjust key parameters (e.g., number of threads, mismatch tolerance) in the config.yaml.
    * Output: Per sample BAM files (per_sample.BAM).

4.12.3. Map Local Reads:
    * Python module to perform read mapping of local FASTQ files to the reference genome using the same aligner as for SRA reads.
    * Parameter Tweaking: Parameter settings should be consistent with SRA read mapping and configurable via config.yaml.
    * Output: Per sample BAM files (per_sample.BAM).

4.12.4. Transcriptome Assembly:
    * Python module to assemble transcripts for each sample using a chosen assembler (e.g., StringTie, Cufflinks).
    * Parameter Tweaking: Allow users to specify assembler and adjust key parameters (e.g., minimum read coverage, minimum isoform fraction) in the config.yaml.
    * Output: Per sample GTF files (per_sample.GTF).

4.12.5. Merge Transcriptomes:
    * Python module to merge per-sample transcriptomes into a non-redundant transcriptome annotation using tools like StringTie merge or custom Python scripts.
    * Output: Non-redundant GTF file (non_redundant_transcriptomes.GTF).

4.12.6. Abundance Estimation:
    * Python module to estimate transcript abundances (e.g., FPKM - Fragments Per Kilobase of transcript per Million mapped reads) using tools like StringTie or RSEM.
    * Output: Coverage table in a human-readable format (e.g., TSV/CSV) with FPKM values (coverage_table_FPKM.tsv).

4.12.7. Coding Probability:
    * Integrate CPAT (Coding Potential Assessment Tool) using Python subprocess calls or a direct Python implementation if available.
    * Assess the coding probability of assembled transcripts.
    * Output: Coding probability scores added as attributes to the non-redundant GTF file or in a separate table.

4.12.8. Transcript Categorization:
    * Develop a Python module to categorize transcripts based on coding probability, gene model, and annotation information.
    * Categories should include: new gene isoforms, new microORFs (coding), new lincRNAs, intronic lncRNAs, antisense lncRNAs (non-coding).
    * Output: Per category GTF files (per_category.GTF) for each category of transcripts.

4.12.9. Selector (Differential Expression Analysis - Optional):
    * Python logic to conditionally execute differential expression analysis based on the presence of a sample metadata file in the configuration.
    * If metadata is provided, perform differential expression analysis using tools like DESeq2 or edgeR (via Python bindings or R subprocess calls).
    * Parameter Tweaking: Allow users to configure parameters for differential expression analysis (e.g., statistical method, p-value adjustment).
    * Output: Table with differential expression analysis results including q-values and fold-change (table_diff_expr.tsv).
    * Sample Metadata File (metadata.csv - optional): CSV file to map samples to experimental conditions for differential expression analysis. Columns should include: sample_id, condition.

4.13. Artificial Intelligence (AI) Stage

This stage will also be implemented in Python, focusing on feature extraction and ML prediction.

4.13.1. Feature Extraction (UCSC):
    * Python module to extract relevant features for candidate lncRNA transcripts from UCSC Genome Browser tracks (BigBed and BigWig files).
    * Features may include: conservation scores, chromatin accessibility, histone modifications, etc. (Specific features should be detailed in the paper and considered for configurability).
    * Improve download speed by parallelizing track data retrieval.
    * Output: Intermediate files storing extracted features per transcript.

4.13.2. Feature Matrix:
    * Python module to combine extracted features into a unified feature matrix.
    * Each row represents a candidate lncRNA transcript, and each column represents a feature.
    * Output: Feature matrix in a standardized format (e.g., CSV or Parquet) (feature_matrix.csv or feature_matrix.parquet).

4.13.3. ML Prediction:
    * Integrate the pre-trained ML model (specify model format and library, e.g., scikit-learn, TensorFlow, PyTorch).
    * Load the feature matrix and feed it into the ML model for prediction.
    * Classify transcripts as lncRNAs or not based on the model's output.
    * Output: Final results table in a human-readable format (e.g., TSV/CSV) with prediction scores and classifications (final_results_table.tsv). This table should include vital information for exploring, filtering, and analyzing candidate lncRNA genes, potentially including:
        * Transcript ID
        * Category (from Transcript Categorization)
        * Coding Probability
        * Differential Expression Analysis results (if performed)
        * ML Prediction Score
        * lncRNA classification (Yes/No)
        * Links to relevant genome browsers (e.g., UCSC, Ensembl)

4.14. Data Sources

4.14.1. Ensembl:
    * Python scripts to download reference genome sequence (FASTA) and annotation (GTF) files directly from Ensembl using their APIs or FTP/HTTP access.
    * Implement checks to ensure Ensembl source reachability before pipeline execution.
    * Improve download speed by using parallel downloads and efficient file handling.

4.14.2. Sequence Read Archive (SRA):
    * Utilize Python libraries (e.g., sra-tools Python bindings if available, or subprocess calls to sra-tools executables) for downloading reads based on SRA accession numbers.
    * Implement error handling for SRA download failures and network issues.

4.14.3. University of California Santa Cruz (UCSC) Genome Browser:
    * Develop Python scripts to access and download BigBed and BigWig files from UCSC Genome Browser using their public APIs or HTTP access.
    * Implement checks to ensure UCSC source reachability.
    * Improve download speed using parallel requests and efficient data parsing.

4.15. Input Data

4.15.1. Input Reads:
    * Support for local FASTQ files (gzip or plain text).
    * Support for SRA accession numbers for direct download from the Sequence Read Archive.

4.15.2. Configuration File (config.yaml):
    * YAML-based configuration file for job submission, offering an abstract way to define parameters.

4.15.3. Sample Metadata File (metadata.csv - optional):
    * CSV file to map samples to experimental conditions for differential expression analysis.
    * Columns should include: sample_id, condition.

5. Non-Functional Requirements

5.1. Performance

NFR-5.1.1. Speed: FLYNC v2 should execute significantly faster than v1, especially for large datasets, through parallelism and optimized data handling.

NFR-5.1.2. Resource Efficiency: The application should be resource-efficient in terms of CPU and memory usage, especially for memory-intensive steps like transcript assembly.

5.2. Reliability

NFR-5.2.1. Robustness: FLYNC v2 should be robust and handle errors gracefully, preventing crashes and providing informative error messages.

NFR-5.2.2. Data Integrity: Ensure data integrity throughout the pipeline, preventing data corruption or loss.

NFR-5.2.3. External Source Resilience: The application should be resilient to temporary unavailability or issues with external data sources.

5.3. Usability

NFR-5.3.1. Ease of Use: FLYNC v2 should be user-friendly and easy to use, even for users with limited command-line experience.

NFR-5.3.2. Clear Documentation: Provide comprehensive and well-structured documentation, including installation instructions, usage examples, parameter descriptions, and troubleshooting guides.

NFR-5.3.3. Informative Output: Output files and reports should be informative, well-organized, and easy to interpret.

5.4. Maintainability

NFR-5.4.1. Code Readability: The Python codebase should be well-structured, modular, and follow coding best practices for readability and maintainability.

NFR-5.4.2. Testability: The application should be designed for easy testing, with unit tests and integration tests covering all critical functionalities.

NFR-5.4.3. Extensibility: The architecture should be extensible, allowing for easy addition of new features, data sources, and ML models in the future.

5.5. Portability

NFR-5.5.1. Cross-Platform Compatibility: FLYNC v2 should be compatible with major operating systems (Linux, macOS). Windows compatibility is desirable but not mandatory.

NFR-5.5.2. Dependency Management: Use a robust dependency management system (e.g., pip, conda) to ensure easy installation and management of required Python packages and external tools. Containerization (Docker/Singularity) is highly recommended to further enhance portability and reproducibility.

6. Technical Specifications

6.1. Technology Stack

Programming Language: Python 3.x
Core Libraries:
    * logging: For application logging.
    * argparse or Click: For command-line interface (CLI) argument parsing.
    * PyYAML: For configuration file parsing.
    * requests or urllib: For downloading data from external sources.
    * Aiohttp: For efficient HTTP requests to external resources (Ensembl, UCSC, SRA).
    * Multiprocessing or threading: For parallelism and multi-threading.
    * pandas: For data manipulation and output formatting (e.g., for feature matrix and final results table).
    * NumPy: For numerical computations and array operations.
    * scikit-learn: For ML model prediction and potentially for feature preprocessing (if directly integrated).
    * Biopython: For bioinformatics file parsing and manipulation.
Bioinformatics Tools:
    * Use an appropriate read mapper (HISAT2). These may be called as subprocesses from Python.
    * Use an appropriate transcriptome assembler (StringTie, Cufflinks). These may be called as subprocesses from Python.
    * CPAT (Coding Potential Assessment Tool) - integrate as a subprocess call or reimplement in Python if feasible.
    * Differential expression analysis tools (e.g., DESeq2, edgeR) - integrate via Python bindings or R subprocess calls.
Workflow Management System (Recommended Future Enhancement): Snakemake or Nextflow.

6.2. Development Environment:

* Version Control: Git (GitHub, GitLab, or similar).
* Issue Tracking: GitHub Issues, Jira, or similar.
* Continuous Integration/Continuous Deployment (CI/CD): GitHub Actions, GitLab CI, or similar (Recommended).

6.3. Deployment:

* Package distribution via pip or conda (Recommended).
* Containerization (Docker/Singularity) for easy deployment and reproducibility (Highly Recommended).

7. Data Sources and Formats

7.1. Input Data:

* Sequencing Reads:
    * FASTQ files (gzipped or uncompressed).
    * SRA Accession Numbers.
* Sample Metadata (Optional): CSV file (@metadata.csv) mapping samples to conditions.
* Configuration File: YAML file (@config.yaml).
* Test List File: Text file (@test-list.txt) for testing.

7.2. External Data Sources:

* Ensembl: For reference genome sequence and annotation (GTF/GFF3, FASTA).
* Sequence Read Archive (SRA): For retrieving published transcriptomic reads (FASTQ).
* UCSC Genome Browser: For BigBed and BigWig files used for feature extraction.

7.3. Output Data:

* Per sample BAM files: Aligned reads in BAM format.
* Per sample GTF files: Transcriptome assembly for each sample in GTF format.
* Non-redundant GTF file: Merged transcriptome annotation in GTF format.
* Coverage table (FPKM): Transcript abundance estimates in a tabular format (CSV/TSV).
* Per category GTF file: Transcripts categorized into coding and non-coding subclasses in GTF format.
* Table with q-value & fold-change: Differential expression analysis results in a tabular format (CSV/TSV).
* Feature matrix: Feature values per candidate lncRNA transcript in a tabular format (CSV/TSV or Parquet).
* Final results table: ML prediction results, including lncRNA classification and relevant information, in a tabular format (CSV/TSV). This table should include:
    * Transcript ID
    * Category
    * Coding Probability
    * Differential Expression Results (if performed)
    * ML Prediction Score
    * lncRNA classification (Yes/No)
    * Links to genome browsers (e.g., UCSC, Ensembl)
* Log files: Application logs in text format.
* Summary reports: Pipeline execution summary in text or HTML format.

8. Error Handling and Checkpointing

8.1. Error Handling:

* Exception Handling: Implement comprehensive exception handling throughout the Python codebase to catch potential errors.
* Informative Error Messages: Provide clear and informative error messages to the user, indicating the type of error, location, and potential causes.
* Graceful Exit: In case of critical errors, the application should exit gracefully, logging the error and providing instructions to the user.
* Retry Mechanisms (Optional): For transient errors (e.g., network issues), consider implementing retry mechanisms with exponential backoff.

8.2. Checkpointing:

* Checkpoint Data: Store necessary data at checkpoint stages to allow for pipeline resumption. This may include intermediate files, progress status, and relevant variables.
* Checkpoint File Format: Use a suitable format for checkpoint files (e.g., serialized Python objects, JSON, or simple text files).
* Checkpoint Management: Implement logic to create, load, and manage checkpoint files.
* Checkpoint Frequency: Define appropriate checkpoint frequency to balance performance and recovery capability.

9. Configuration and Parameters

9.1. Configuration File (config.yaml):

* Sections: Organize parameters into logical sections (input, output, genome, mapping, assembly, ml_model, etc.).
* Parameters: Include parameters for:
    * Input data paths (FASTQ files, SRA accessions, metadata file).
    * Output directory.
    * Reference genome source (Ensembl release, genome assembly).
    * Read mapper selection and parameters.
    * Transcript assembler selection and parameters.
    * CPAT parameters (if configurable).
    * ML model path (if configurable).
    * Parallelism settings (number of threads/processes).
    * Logging level.
    * Checkpointing options.
    * UCSC data source parameters (genome browser track names).
    * Differential expression analysis settings (optional).
    * Feature extraction settings (UCSC tracks).

9.2. Command-line Arguments:

* Essential Arguments: Provide command-line arguments for essential parameters like configuration file path, input data paths, and output directory.
* Parameter Overrides: Allow command-line arguments to override parameters defined in the configuration file.
* Help Messages: Provide comprehensive help messages for command-line arguments using argparse or Click.

10. Testing and Validation

10.1. Unit Tests:

* Develop unit tests for individual Python modules and functions to ensure they function correctly in isolation.
* Use a testing framework like pytest or unittest.

10.2. Integration Tests:

* Develop integration tests to verify the correct interaction between different modules and pipeline stages.
* Use the provided @metadata.csv and @test-list.txt files as test datasets.
* Test different pipeline configurations and parameter settings.
* System Tests: Run end-to-end pipeline tests using provided test data (@metadata.csv, @test-list.txt) and realistic datasets.

10.3. Validation Data:

* Use known datasets and benchmarks to validate the accuracy and performance of FLYNC v2.
* Compare results with FLYNC v1 and other existing lncRNA identification tools.

10.4. Continuous Integration (CI):

* Implement CI using GitHub Actions or GitLab CI to automatically run tests on code changes and ensure code quality and stability.
* Implement a comprehensive testing strategy: Unit Tests, Integration Tests, System Tests.
* Use a testing framework (e.g., unittest, pytest) for automated test execution and reporting.

11. Future Enhancements

* Workflow Management System Integration: Migrate the pipeline to a workflow management system like Snakemake or Nextflow for improved scalability, reproducibility, and pipeline management. This would handle dependency management, parallel execution, and checkpointing more robustly.
* Containerization: Package FLYNC v2 as a Docker or Singularity container to ensure reproducibility and simplify deployment across different environments.
* Cloud Deployment (Optional): Explore cloud deployment options (e.g., AWS, Google Cloud, Azure) for scalability and accessibility. Enable deployment and execution of FLYNC v2 on cloud platforms (e.g., AWS, GCP, Azure) for scalability and accessibility.
* Web Interface/GUI (Optional): Consider developing a web-based or graphical user interface (GUI) for easier access and usage for non-command-line users. Develop a web interface or graphical user interface (GUI) to make FLYNC more accessible to users who are not comfortable with command-line tools.
* Expanded Data Source Support: Integrate support for additional data sources beyond Ensembl, SRA, and UCSC. Extend support to more external data sources for reference genomes, annotations, and functional genomics data.
* Advanced Visualization: Integrate tools for visualizing results, such as interactive plots of transcript expression, lncRNA classifications, and genome browser tracks.
* Plugin Architecture: Design a plugin architecture to allow users to easily extend FLYNC v2 with custom modules and functionalities.
* Model Retraining Pipeline: Develop a pipeline for retraining the ML model with new features and datasets. Provide a separate script or module for retraining the ML model. Explore automated model retraining pipelines and feature selection methods.
* Feature Selection and Engineering: Explore new features and feature selection methods to improve model accuracy.
* Automated Model Retraining Pipeline: Develop an automated pipeline for continuous model retraining with new data and features, including methods for feature selection and model evaluation.

12. Conclusion

FLYNC v2, with its Python-first approach and the enhancements outlined in this Product Specification Document, will represent a significant modernization and improvement over the original application. By addressing the identified requirements and incorporating future enhancements, FLYNC v2 will become a more robust, efficient, user-friendly, and maintainable tool for lncRNA discovery, empowering researchers in the field of non-coding RNA biology. This document provides a detailed roadmap for developers to create a significantly improved application ready for peer review and wider adoption. Regular communication and feedback are encouraged to ensure the final product meets the desired requirements and goals.