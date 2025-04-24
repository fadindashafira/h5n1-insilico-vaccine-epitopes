#!/usr/bin/env nextflow

/*
 * Module for molecular dynamics simulation
 */

// Molecular dynamics simulation of vaccine construct
process molecularDynamics {
    tag "${protein_type}:md_simulation"
    publishDir "${params.experiment_output}/molecular_dynamics/${protein_type}", mode: params.publish_dir_mode, overwrite: true
    
    input:
    path vaccine
    val protein_type
    
    output:
    path "md_results", type: 'dir', emit: md_results
    path "md_report.txt", emit: md_report
    
    script:
    """
    python3 ${projectDir}/bin/md_simulation.py \\
        --vaccine ${vaccine} \\
        --protein-type ${protein_type} \\
        --temperature ${params.md_temperature ?: '310'} \\
        --simulation_time ${params.md_time ?: '10'} \\
        --force_field "${params.md_forcefield ?: 'amber99sb-ildn'}" \\
        --t_helper "${params.t_helper ?: ''}" \\
        --linkers ${params.linker ? "\"${params.linker}\"" : '""'}
    """
}