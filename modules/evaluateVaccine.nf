#!/usr/bin/env nextflow

/*
 * Module for vaccine evaluation processes
 */

// Evaluate the vaccine construct
process evaluateVaccineConstruct {
    tag "${protein_type}:evaluate_vaccine"
    publishDir "${params.experiment_output}/evaluation/${protein_type}", mode: params.publish_dir_mode, overwrite: true
    
    input:
    path vaccine
    val protein_type
    
    output:
    path "${protein_type}_vaccine_evaluation.txt", emit: evaluation_report
    path "${protein_type}_vaccine_properties.csv", emit: properties
    path "${protein_type}_vaccine_for_colabfold.fasta", emit: colabfold_fasta
    
    script:
    """
    python ${workflow.projectDir}/bin/evaluate_vaccine.py \\
        --vaccine=${vaccine} \\
        --protein-type=${protein_type} \\
        --linker=${params.linker ?: 'GPGPG'} \\
        --iedb-api-url=${params.iedb_api_url ?: 'http://tools-api.iedb.org/tools_api/'} \\
        --output-evaluation=${protein_type}_vaccine_evaluation.txt \\
        --output-properties=${protein_type}_vaccine_properties.csv \\
        --output-colabfold=${protein_type}_vaccine_for_colabfold.fasta
    """
}