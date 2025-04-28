#!/usr/bin/env nextflow

/*
 * Module for vaccine design processes
 */

// Combine epitopes from different prediction methods
process combineEpitopes {
    tag "${protein_type}:combine_epitopes"
    publishDir "${params.experiment_output}/epitopes/${protein_type}", mode: params.publish_dir_mode
    
    input:
    path bcell_epitopes
    path tcell_i_epitopes
    path tcell_ii_epitopes
    val protein_type
    
    output:
    path "${protein_type}_combined_epitopes.csv", emit: combined_epitopes
    
    script:
    """
    Rscript ${workflow.projectDir}/bin/combine_epitopes.R \\
        --bcell=${bcell_epitopes} \\
        --tcell-i=${tcell_i_epitopes} \\
        --tcell-ii=${tcell_ii_epitopes} \\
        --protein-type=${protein_type} \\
        --output=${protein_type}_combined_epitopes.csv
    """
}

// Design vaccine construct from top epitopes
process designVaccineConstruct {
    tag "${protein_type}:design_vaccine"
    publishDir "${params.experiment_output}/vaccine/${protein_type}", mode: params.publish_dir_mode, overwrite: true
    
    input:
    path combined_epitopes
    val protein_type
    
    output:
    path "${protein_type}_vaccine_construct.fasta", emit: vaccine
    path "${protein_type}_vaccine_report.html", emit: report
    
    script:
    // If multiple files are provided, create a single combined file first
    if (combined_epitopes instanceof List && combined_epitopes.size() > 1) {
        """
        # First, concatenate the epitope files
        head -n 1 ${combined_epitopes[0]} > combined_epitopes_all.csv
        for f in ${combined_epitopes.join(" ")}; do
            tail -n +2 \$f >> combined_epitopes_all.csv
        done
        
        # Now run with the combined file
        python ${workflow.projectDir}/bin/design_vaccine.py \
            --combined-epitopes=combined_epitopes_all.csv \
            --protein-type=${protein_type} \
            --linker=${params.linker ?: 'GPGPG'} \
            --max-epitopes=${params.max_epitopes ?: 10} \
            --min-epitopes=${params.min_epitopes ?: 1} \
            --leading-seq=${params.leading_seq ?: ''} \
            --trailing-seq=${params.trailing_seq ?: ''} \
            --similarity-threshold=${params.similarity_threshold ?: 0.7} \
            --output-fasta=${protein_type}_vaccine_construct.fasta \
            --output-report=${protein_type}_vaccine_report.html
        """
    } else {
        // Original script for single file
        """
        python ${workflow.projectDir}/bin/design_vaccine.py \
            --combined-epitopes=${combined_epitopes} \
            --protein-type=${protein_type} \
            --linker=${params.linker ?: 'GPGPG'} \
            --max-epitopes=${params.max_epitopes ?: 10} \
            --min-epitopes=${params.min_epitopes ?: 1} \
            --leading-seq=${params.leading_seq ?: ''} \
            --trailing-seq=${params.trailing_seq ?: ''} \
            --similarity-threshold=${params.similarity_threshold ?: 0.7} \
            --output-fasta=${protein_type}_vaccine_construct.fasta \
            --output-report=${protein_type}_vaccine_report.html
        """
    }
}