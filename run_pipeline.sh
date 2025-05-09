#!/bin/bash

# run_pipeline.sh - Script to run the in-silico vaccine design pipeline
# Usage: ./run_pipeline.sh [options]

# Exit on error
set -e

# Default values
ACCESSION1="ABK40530.1"
ACCESSION2="NP_740664.1"
OUTPUT_DIR="results"
EXPERIMENT_ID="exp1"
PROFILE="local"
RUN_MD=false
RESUME=false
HELP=false
OVERWRITE_DAG=true
VERBOSE=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --accession1)
      ACCESSION1="$2"
      shift 2
      ;;
    --accession2)
      ACCESSION2="$2"
      shift 2
      ;;
    --outdir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --experiment_id)
      EXPERIMENT_ID="$2"
      shift 2
      ;;
    --profile)
      PROFILE="$2"
      shift 2
      ;;
    --run_md)
      RUN_MD=true
      shift
      ;;
    --resume)
      RESUME=true
      shift
      ;;
    --no-dag)
      OVERWRITE_DAG=false
      shift
      ;;
    --verbose)
      VERBOSE=true
      shift
      ;;
    --help)
      HELP=true
      shift
      ;;
    *)
      echo "Unknown option: $1"
      HELP=true
      shift
      ;;
  esac
done

# Display help message
if [ "$HELP" = true ]; then
  echo "In-silico Vaccine Design Pipeline"
  echo "=================================="
  echo "Usage: ./run_pipeline.sh [options]"
  echo ""
  echo "Options:"
  echo "  --accession1 VALUE     accession number 1 (default: $ACCESSION1)"
  echo "  --accession2 VALUE     accession number 2 (default: $ACCESSION2)"
  echo "  --outdir VALUE         Base output directory (default: $OUTPUT_DIR)"
  echo "  --experiment_id VALUE  Experiment identifier (default: $EXPERIMENT_ID)"
  echo "  --profile VALUE        Configuration profile (default: $PROFILE)"
  echo "  --run_md               Run molecular dynamics simulation"
  echo "  --resume               Resume pipeline execution from the last checkpoint"
  echo "  --no-dag               Don't generate/overwrite DAG file"
  echo "  --verbose              Enable verbose logging"
  echo "  --help                 Display this help message"
  echo ""
  exit 0
fi

# Create base output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Also create the experiment-specific results directory and necessary subdirectories
EXPERIMENT_OUTPUT_DIR="${OUTPUT_DIR}/results_${EXPERIMENT_ID}"
mkdir -p "$EXPERIMENT_OUTPUT_DIR"
mkdir -p "$EXPERIMENT_OUTPUT_DIR/reports"
mkdir -p "$EXPERIMENT_OUTPUT_DIR/pipeline_info"

# Ensure bin directory exists for scripts
mkdir -p bin

# Check if necessary scripts exist
REQUIRED_SCRIPTS=(
  "bin/predict_bcell.py"
  "bin/predict_tcell_i.py"
  "bin/predict_tcell_ii.py"
  "bin/filter_epitopes.py"
  "bin/analyze_conservation.py"
  "bin/design_vaccine.py"
  "bin/evaluate_vaccine.py"
)

MISSING_SCRIPTS=()
for script in "${REQUIRED_SCRIPTS[@]}"; do
  if [ ! -f "$script" ]; then
    MISSING_SCRIPTS+=("$script")
  fi
done

if [ ${#MISSING_SCRIPTS[@]} -gt 0 ]; then
  echo "WARNING: The following required scripts are missing:"
  for script in "${MISSING_SCRIPTS[@]}"; do
    echo "  - $script"
  done
  echo "Please ensure all necessary scripts are in the bin directory."
  read -p "Do you want to continue anyway? (y/n) " -n 1 -r
  echo
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    exit 1
  fi
fi

# Create command string
CMD="nextflow run main.nf"
CMD+=" --accession1 '$ACCESSION1'"
CMD+=" --accession2 '$ACCESSION2'"
CMD+=" --outdir '$OUTPUT_DIR'"
CMD+=" --experiment_id '$EXPERIMENT_ID'"

# Add run_md flag if specified
if [ "$RUN_MD" = true ]; then
  CMD+=" --run_md"
fi

# Add resume flag if specified
if [ "$RESUME" = true ]; then
  CMD+=" -resume"
fi

# Add profile
CMD+=" -profile $PROFILE"

# Always add the reporting options to ensure consistent file locations
CMD+=" -with-report ${EXPERIMENT_OUTPUT_DIR}/reports/execution_report.html"
CMD+=" -with-timeline ${EXPERIMENT_OUTPUT_DIR}/reports/timeline.html"
CMD+=" -with-trace ${EXPERIMENT_OUTPUT_DIR}/reports/trace.txt"

# Add DAG generation parameters
if [ "$OVERWRITE_DAG" = true ]; then
  # Use the experiment-specific reports directory for DAG
  DAG_FILE="${EXPERIMENT_OUTPUT_DIR}/reports/dag_${EXPERIMENT_ID}_$(date +%Y%m%d_%H%M%S).html"
  CMD+=" -with-dag $DAG_FILE"
fi

# Print the command
echo "Running: $CMD"
echo "=================================================="
echo "Accession 1: $ACCESSION1"
echo "Accession 2: $ACCESSION2"
echo "Base output directory: $OUTPUT_DIR"
echo "Experiment ID: $EXPERIMENT_ID"
echo "Experiment output directory: $EXPERIMENT_OUTPUT_DIR"
echo "Run molecular dynamics: $RUN_MD"
echo "Resume execution: $RESUME"
echo "=================================================="

# Execute the command
eval $CMD

# Check execution status
if [ $? -eq 0 ]; then
  echo "=================================================="
  echo "Pipeline execution completed successfully!"
  echo "Results available in: $EXPERIMENT_OUTPUT_DIR"
  
  # List the generated output files
  echo "Generated files:"
  find "$EXPERIMENT_OUTPUT_DIR" -type f | sort
  
  if [ "$OVERWRITE_DAG" = true ]; then
    echo "DAG visualization available at: $DAG_FILE"
  fi
  
  echo "=================================================="
else
  echo "=================================================="
  echo "Pipeline execution failed."
  echo "Check the logs for more details."
  echo "Nextflow work directory: $(pwd)/work"
  echo "=================================================="
  exit 1
fi