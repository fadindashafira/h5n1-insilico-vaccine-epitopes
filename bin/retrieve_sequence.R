#!/usr/bin/env Rscript

# Get command line arguments directly from args
args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
accession <- args[1]
protein_type <- args[2]
api_key <- args[3]
retry_delay <- as.numeric(args[4])
max_retries <- as.numeric(args[5])

# Install packages if not available
if (!requireNamespace("rentrez", quietly = TRUE)) {
    install.packages("rentrez", repos = "http://cran.us.r-project.org")
}

library(rentrez)

# Set up error handling
options(warn = 2)

# Set API key if provided
if (!is.null(api_key) && nchar(api_key) > 0) {
    set_entrez_key(api_key)
    cat("Using NCBI API key\n")
}

# Retrieve sequence for given accession number
cat("Retrieving", protein_type, "sequence for accession number:", accession, "\n")

# Add retry mechanism with exponential backoff
success <- FALSE

for (attempt in 1:max_retries) {
    if (attempt > 1) {
        cat("Attempt", attempt, "of", max_retries, "- waiting", retry_delay, "seconds...\n")
        Sys.sleep(retry_delay)
        retry_delay <- retry_delay * 2  # Exponential backoff
    }
    
    tryCatch({
        # First, check if the accession exists
        search_results <- entrez_search(db = "protein", term = accession)
        
        if (search_results$count == 0) {
            stop("No ", protein_type, " sequence found for accession: ", accession)
        }
        
        sequence_data <- entrez_fetch(
            db = "protein", 
            id = accession, 
            rettype = "fasta", 
            retmode = "text"
        )
        
        if (nchar(sequence_data) < 10) {
            stop("Retrieved ", protein_type, " sequence is too short or empty")
        }
        
        # Save the sequence to a FASTA file
        fasta_file <- paste0(protein_type, "_", accession, ".fasta")
        writeLines(sequence_data, fasta_file)
        
        cat(protein_type, "sequence downloaded and saved to", fasta_file, "\n")
        success <- TRUE
        break  # Exit the retry loop if successful
        
    }, error = function(e) {
        if (attempt == max_retries) {
            cat("ERROR: Failed to retrieve ", protein_type, " sequence after ", max_retries, " attempts. ", e$message, "\n", file = stderr())
            quit(status = 1)
        } else {
            cat("Warning: Attempt", attempt, "failed:", e$message, "\n")
        }
    })
}

if (!success) {
    cat("ERROR: Failed to retrieve sequence after all retry attempts\n", file = stderr())
    quit(status = 1)
}