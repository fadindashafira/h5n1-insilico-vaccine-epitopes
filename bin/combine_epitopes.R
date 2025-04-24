#!/usr/bin/env Rscript

# Parse command line arguments
library(optparse)

option_list <- list(
  make_option("--bcell", type="character", help="B-cell epitopes CSV file"),
  make_option("--tcell-i", type="character", help="T-cell Class I epitopes CSV file"),
  make_option("--tcell-ii", type="character", help="T-cell Class II epitopes CSV file"),
  make_option("--protein-type", type="character", help="Protein type (e.g., hemagglutinin, neuraminidase)"),
  make_option("--output", type="character", help="Output CSV file for combined epitopes")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check required arguments
required_args <- c("bcell", "tcell-i", "tcell-ii", "protein-type", "output")
missing_args <- required_args[!required_args %in% names(opt)]
if (length(missing_args) > 0) {
  stop("Missing required arguments: ", paste(missing_args, collapse=", "))
}

# Load required libraries
if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr", repos = "http://cran.us.r-project.org")
}
library(dplyr)

# Function to read and standardize CSV files
read_and_standardize_csv <- function(file_path) {
    # Read the CSV file
    df <- tryCatch({
        read.csv(file_path, stringsAsFactors = FALSE)
    }, error = function(e) {
        warning(paste("Failed to read file:", file_path, "-", e$message))
        # Return empty dataframe with correct structure
        return(data.frame(
            sequence = character(0),
            start = integer(0),
            end = integer(0),
            score = numeric(0),
            type = character(0),
            method = character(0),
            source = character(0),
            hla = character(0),
            ic50 = numeric(0),
            stringsAsFactors = FALSE
        ))
    })
    
    # Ensure consistent column names
    standard_columns <- c(
        "sequence", "start", "end", "score", "type", 
        "method", "source", "hla", "ic50"
    )
    
    # Create a dataframe with all standard columns
    result_df <- data.frame(
        sequence = df$sequence %||% NA_character_,
        start = df$start %||% NA_integer_,
        end = df$end %||% NA_integer_,
        score = df$score %||% NA_real_,
        type = df$type %||% NA_character_,
        method = df$method %||% NA_character_,
        source = df$source %||% NA_character_,
        hla = df$hla %||% NA_character_,
        ic50 = df$ic50 %||% NA_real_,
        stringsAsFactors = FALSE
    )
    
    return(result_df)
}

# Safe null coalescing operator
`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

# Log the protein type
protein_type <- opt$`protein-type`
cat("Combining epitopes for protein:", protein_type, "\n")

# Read all epitope files
cat("Reading B-cell epitopes from:", opt$bcell, "\n")
bcell_df <- read_and_standardize_csv(opt$bcell)

cat("Reading T-cell Class I epitopes from:", opt$`tcell-i`, "\n")
tcell_i_df <- read_and_standardize_csv(opt$`tcell-i`)

cat("Reading T-cell Class II epitopes from:", opt$`tcell-ii`, "\n")
tcell_ii_df <- read_and_standardize_csv(opt$`tcell-ii`)

# Combine all epitopes
combined_df <- rbind(bcell_df, tcell_i_df, tcell_ii_df)

# Make sure all epitopes are tagged with the correct protein type
combined_df$source <- protein_type

# Calculate a consensus score for ranking
combined_df <- combined_df %>%
    mutate(
        consensus_score = case_when(
            type == "B-cell" ~ score,
            !is.na(ic50) ~ 1 - (pmin(ic50, 5000) / 5000),
            TRUE ~ score
        )
    )

# Sort by consensus score (descending)
combined_df <- combined_df %>%
    arrange(desc(consensus_score))

# Save combined file
write.csv(combined_df, opt$output, row.names = FALSE)

# Print summary
cat("Epitopes combined from multiple prediction methods for", protein_type, "\n")
cat("Total epitopes:", nrow(combined_df), "\n")
cat("Epitope types:\n")
print(table(combined_df$type))