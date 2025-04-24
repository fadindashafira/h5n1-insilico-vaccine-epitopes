# In Silico Vaccine Design Automation

A NextFlow-based pipeline for automated computational vaccine design targeting avian influenza viruses, based on epitope prediction and molecular dynamics simulations.

## Overview

This pipeline automates the process of designing epitope-based vaccines for H5N1 avian influenza virus using computational methods. The workflow implements the methodology described in:

> Tambunan et al. Vaccine Design for H5N1 Based on B- and T-cell Epitope Predictions. *Bioinformatics and Biology Insights* 2016:10 27-35.

The pipeline incorporates the following steps:
1. Retrieve protein sequences from NCBI
2. Predict B-cell epitopes
3. Predict T-cell epitopes (MHC class I and II binding)
4. Filter and select high-quality epitopes
5. Design a multi-epitope vaccine construct with appropriate linkers
6. Evaluate vaccine properties
7. Optional molecular dynamics simulation

## Requirements

- Nextflow (>=21.10)
- Container technology (optional but recommended):
  - Docker or Singularity
- R (>=4.0) packages:
  - Biostrings
  - rentrez
  - seqinr
  - httr
- Python (>=3.8) packages:
  - biopython
  - pandas
  - numpy
  - requests

## Quick Start

```bash
# Clone this repository
git clone [https://github.com/username/in-silico-vaccine-design.git](https://github.com/fadindashafira/h5n1-insilico-vaccine-epitopes.git)
cd h5n1-insilico-vaccine-epitopes

# Run with default parameters
./run_pipeline.sh

# Run with custom parameters
./run_pipeline.sh --ha_accession "BAL61222.1" --na_accession "BAL61230.1" --experiment_id "exp2"

# Run on a SLURM cluster
./run_pipeline.sh --profile slurm
```

## Pipeline Parameters

You can run the pipeline with different parameters:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--ha_accession` | Hemagglutinin protein accession number | BAL61222.1 |
| `--na_accession` | Neuraminidase protein accession number | BAL61230.1 |
| `--outdir` | Base output directory | results |
| `--experiment_id` | Unique identifier for the experiment | exp1 |
| `--profile` | Execution profile (local, slurm, etc.) | local |
| `--run_md` | Run molecular dynamics simulation | false |
| `--resume` | Resume pipeline execution from last checkpoint | false |
| `--verbose` | Enable verbose logging | false |

## Output Structure

The pipeline generates a structured output directory for each experiment:

```
results/
└── results_<experiment_id>/
    ├── reports/                    # Pipeline execution reports
    │   ├── dag.html                # Workflow visualization
    │   ├── execution_report.html   # Execution statistics
    │   ├── timeline.html           # Execution timeline
    │   └── trace.txt               # Execution trace
    ├── sequences/                  # Retrieved protein sequences
    ├── epitopes/                   # Predicted epitopes
    │   ├── hemagglutinin/
    │   │   ├── bcell/              # B-cell epitopes
    │   │   ├── tcell_i/            # T-cell MHC class I epitopes  
    │   │   ├── tcell_ii/           # T-cell MHC class II epitopes
    │   │   └── filtered/           # Filtered epitopes
    │   └── neuraminidase/
    │       └── ...
    ├── vaccine_constructs/         # Designed vaccine sequences
    │   ├── hemagglutinin_vaccine.fasta
    │   ├── neuraminidase_vaccine.fasta
    │   └── combined_vaccine.fasta
    ├── vaccine_evaluation/         # Evaluation reports
    └── molecular_dynamics/         # MD simulation results (if enabled)
```

## Customization

You can customize the pipeline by editing the `nextflow.config` file:

- HLA alleles for T-cell epitope prediction
- Epitope prediction methods and parameters
- Linker sequences for vaccine construction
- Resource allocation for high-performance computing

```bash
# Example of running with custom configuration
./run_pipeline.sh --experiment_id custom_alleles
```

## Advanced Usage

### Using different prediction methods

The pipeline supports multiple epitope prediction methods:

- B-cell: Bepipred-1.0, Bepipred-2.0, Chou-Fasman, Emini, Karplus-Schulz, Kolaskar-Tongaonkar, Parker
- T-cell: NetMHCpan, NetMHCIIpan

You can modify these in the nextflow.config file.

### Running on HPC clusters

The pipeline is designed to work with high-performance computing clusters:

```bash
# Run on SLURM cluster
./run_pipeline.sh --profile slurm --experiment_id hpc_run
```

## Troubleshooting

If you encounter issues:

1. Check that all required Python scripts exist in the `bin/` directory
2. Verify that the necessary Python and R packages are installed
3. Check the log files in `results/results_<experiment_id>/reports/` for errors
4. Use the `--verbose` flag for more detailed logging
