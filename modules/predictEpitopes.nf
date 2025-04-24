#!/usr/bin/env nextflow
/*
 * Module for epitope prediction with IEDB integration
 * Optimized for H5N1 vaccine design
 */

// Utility function to safely parse alleles
def parseAlleles(allelesToParse) {
    def parsedAlleles = []
    
    // Handle different input types
    if (allelesToParse instanceof String) {
        // Remove brackets, quotes, and split
        parsedAlleles = allelesToParse
            .replaceAll(/[\[\]\'\"]/, '')
            .split(',')
            .collect { it.trim() }
            .findAll { it }
    } else if (allelesToParse instanceof List) {
        parsedAlleles = allelesToParse
    }
    
    return parsedAlleles
}

// B-cell epitope prediction using IEDB-python
process predictBCellEpitopes {
    tag "${protein_type}:${fasta.baseName}"
    publishDir "${params.experiment_output}/epitopes/${protein_type}/bcell", mode: params.publish_dir_mode, overwrite: true
    
    input:
    path fasta
    val protein_type
    
    output:
    path "${protein_type}_${fasta.baseName}_bcell_epitopes.csv", emit: bcell_epitopes
    
    script:
    // Define B-cell prediction parameters with safer defaults for H5N1
    def method = params.bcell_method ?: 'Bepipred-2.0'
    def threshold = params.bcell_threshold ?: 0.5
    def window_size = params.bcell_length ?: 12  // Larger window for H5N1 epitopes
    
    """
    python ${workflow.projectDir}/bin/predict_bcell.py \\
        --fasta=${fasta} \\
        --protein-type=${protein_type} \\
        --method=${method} \\
        --threshold=${threshold} \\
        --window-size=${window_size} \\
        --output=${protein_type}_${fasta.baseName}_bcell_epitopes.csv
    """
}

// T-cell epitope prediction for MHC Class I
process predictTCellEpitopesI {
    tag "${protein_type}:${fasta.baseName}"
    publishDir "${params.experiment_output}/epitopes/${protein_type}/tcell_i", mode: params.publish_dir_mode, overwrite: true
    
    input:
    path fasta
    val protein_type
    
    output:
    path "${protein_type}_${fasta.baseName}_tcell_i_epitopes.csv", emit: tcell_i_epitopes
    
    script:
    // Convert alleles to a comma-separated string
    def allelesList = parseAlleles(params.mhci_alleles)
    def alleleString = allelesList.join(',')
    
    // Ensure method and threshold are defined with optimal defaults for H5N1
    def method = params.mhci_method ?: 'NetMHCpan'
    def threshold = params.mhci_threshold ?: 500
    def length = params.mhci_length ?: 9
    
    """
    python ${workflow.projectDir}/bin/predict_tcell_i.py \\
        --fasta=${fasta} \\
        --protein-type=${protein_type} \\
        --method=${method} \\
        --threshold=${threshold} \\
        --length=${length} \\
        --alleles='${alleleString}' \\
        --output=${protein_type}_${fasta.baseName}_tcell_i_epitopes.csv
    """
}

// T-cell epitope prediction for MHC Class II
process predictTCellEpitopesII {
    tag "${protein_type}:${fasta.baseName}"
    publishDir "${params.experiment_output}/epitopes/${protein_type}/tcell_ii", mode: params.publish_dir_mode, overwrite: true
    
    input:
    path fasta
    val protein_type
    
    output:
    path "${protein_type}_${fasta.baseName}_tcell_ii_epitopes.csv", emit: tcell_ii_epitopes
    
    script:
    // Convert alleles to a comma-separated string
    def allelesList = parseAlleles(params.mhcii_alleles)
    def alleleString = allelesList.join(',')
    
    // Ensure method and threshold are defined with optimal defaults for H5N1
    def method = params.mhcii_method ?: 'NetMHCIIpan'
    def threshold = params.mhcii_threshold ?: 500
    def length = params.mhcii_length ?: 15
    
    """
    python ${workflow.projectDir}/bin/predict_tcell_ii.py \\
        --fasta=${fasta} \\
        --protein-type=${protein_type} \\
        --method=${method} \\
        --threshold=${threshold} \\
        --length=${length} \\
        --alleles='${alleleString}' \\
        --output=${protein_type}_${fasta.baseName}_tcell_ii_epitopes.csv
    """
}

// Epitope filtering and selection for H5N1 vaccine design
process filterEpitopes {
    tag "${protein_type}"
    publishDir "${params.experiment_output}/epitopes/${protein_type}/filtered", mode: params.publish_dir_mode, overwrite: true
    
    input:
    path bcell_epitopes
    path tcell_i_epitopes
    path tcell_ii_epitopes
    val protein_type
    
    output:
    path "${protein_type}_filtered_epitopes.csv", emit: filtered_epitopes
    
    script:
    """
    python ${workflow.projectDir}/bin/filter_epitopes.py \\
        --bcell=${bcell_epitopes} \\
        --tcell-i=${tcell_i_epitopes} \\
        --tcell-ii=${tcell_ii_epitopes} \\
        --protein-type=${protein_type} \\
        --output=${protein_type}_filtered_epitopes.csv
    """
}

// Optional: Conservation analysis for H5N1 strains
process analyzeConservation {
    tag "${protein_type}"
    publishDir "${params.experiment_output}/conservation/${protein_type}", mode: params.publish_dir_mode, overwrite: true
    
    input:
    path filtered_epitopes
    path alignment
    val protein_type
    
    output:
    path "${protein_type}_conserved_epitopes.csv", emit: conserved_epitopes
    
    script:
    """
    python ${workflow.projectDir}/bin/analyze_conservation.py \\
        --epitopes=${filtered_epitopes} \\
        --alignment=${alignment} \\
        --protein-type=${protein_type} \\
        --threshold=${params.conservation_threshold ?: 0.9} \\
        --output=${protein_type}_conserved_epitopes.csv
    """
}