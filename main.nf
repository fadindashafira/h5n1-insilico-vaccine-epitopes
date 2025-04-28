#!/usr/bin/env nextflow

/*
========================================================================================
    In Silico Vaccine Design Pipeline
========================================================================================
*/

nextflow.enable.dsl = 2

// Set default parameters
params.outdir = 'results'
params.experiment_id = 'exp1'
params.ha_accession = ''
params.na_accession = ''
params.experiment_output = "${params.outdir}/results_${params.experiment_id}"
params.run_md = false
params.help = false

// NCBI API and retrieval parameters
params.ncbi_api_key = ''
params.retry_delay = 3
params.max_retries = 5

// Default parameters for epitope prediction
params.bcell_method = 'Bepipred-1.0'  // Changed from 2.0 to 1.0
params.bcell_threshold = 0.5
params.bcell_length = 9  // Changed from 12 to 9 to match the paper

params.mhci_alleles = ['HLA-A*02:01', 'HLA-B*07:02', 'HLA-B*35:01', 'HLA-A*11:01']
params.mhci_method = 'NetMHCpan'
params.mhci_threshold = 500
params.mhci_length = 9

params.mhcii_alleles = ['HLA-DRB1*03:01', 'HLA-DRB1*04:01', 'HLA-DRB1*01:01', 'HLA-DRB1*07:01']
params.mhcii_method = 'NetMHCIIpan'
params.mhcii_threshold = 500
params.mhcii_length = 15

params.publish_dir_mode = 'copy'
params.linker = 'GGGGS'
params.max_epitopes = 20
params.min_epitopes = 5

// Additional pipeline configuration
params.conservation_threshold = 0.9

// Import modules
include { retrieveSequence as retrieveHA } from './modules/retrieveSequence'
include { retrieveSequence as retrieveNA } from './modules/retrieveSequence'
include { predictBCellEpitopes as predictHABCellEpitopes } from './modules/predictEpitopes'
include { predictTCellEpitopesI as predictHATCellEpitopesI } from './modules/predictEpitopes'
include { predictTCellEpitopesII as predictHATCellEpitopesII } from './modules/predictEpitopes'
include { predictBCellEpitopes as predictNABCellEpitopes } from './modules/predictEpitopes'
include { predictTCellEpitopesI as predictNATCellEpitopesI } from './modules/predictEpitopes'
include { predictTCellEpitopesII as predictNATCellEpitopesII } from './modules/predictEpitopes'
include { combineEpitopes as combineHAEpitopes } from './modules/designVaccine'
include { combineEpitopes as combineNAEpitopes } from './modules/designVaccine'
include { combineEpitopes as combinedEpitopes } from './modules/designVaccine'
include { designVaccineConstruct as designHAVaccine } from './modules/designVaccine'
include { designVaccineConstruct as designNAVaccine } from './modules/designVaccine'
include { designVaccineConstruct as designCombinedVaccine } from './modules/designVaccine'
include { evaluateVaccineConstruct as evaluateHAVaccine } from './modules/evaluateVaccine'
include { evaluateVaccineConstruct as evaluateNAVaccine } from './modules/evaluateVaccine'
include { evaluateVaccineConstruct as evaluateCombinedVaccine } from './modules/evaluateVaccine'
include { molecularDynamics as runHAMD } from './modules/molecularDynamics'
include { molecularDynamics as runNAMD } from './modules/molecularDynamics'
include { molecularDynamics as runCombinedMD } from './modules/molecularDynamics'

// Function to display help message
def helpMessage() {
    log.info"""
    =============================================
    IN SILICO VACCINE DESIGN PIPELINE
    =============================================
    Usage:
    nextflow run main.nf --ha_accession <ha_accession> --na_accession <na_accession> 
    
    Required Parameters:
      --ha_accession      Hemagglutinin protein accession number from NCBI
      --na_accession      Neuraminidase protein accession number from NCBI
    
    Optional Parameters:
      --outdir            Base output directory (default: results)
      --experiment_id     Identifier for the experiment (default: exp1)
      --run_md            Run molecular dynamics simulation (default: false)
      
    Epitope Prediction Options:
      --bcell_method      B-cell epitope prediction method (default: Bepipred-2.0)
      --bcell_threshold   B-cell epitope score threshold (default: 0.5)
      --mhci_method       MHC Class I epitope prediction method (default: NetMHCpan)
      --mhci_threshold    MHC Class I prediction threshold (default: 500)
      --mhcii_method      MHC Class II epitope prediction method (default: NetMHCIIpan)
      --mhcii_threshold   MHC Class II prediction threshold (default: 500)
      
    Advanced Options:
      --ncbi_api_key      NCBI API key for sequence retrieval
      --retry_delay       Delay between retrieval attempts (default: 3)
      --max_retries       Maximum number of retrieval attempts (default: 5)
      --help              Show this help message
    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

// Validate required parameters
if (!params.ha_accession || !params.na_accession) {
    log.error "ERROR: Both --ha_accession and --na_accession must be provided"
    helpMessage()
    exit 1
}

// Warn about missing NCBI API key
if (!params.ncbi_api_key) {
    log.warn "WARNING: No NCBI API key provided. This may limit sequence retrieval capabilities."
}

// Ensure output directory exists
file(params.experiment_output).mkdirs()

// Log pipeline parameters
log.info"""
=============================================
 IN SILICO VACCINE DESIGN PIPELINE
=============================================
 Hemagglutinin accession: ${params.ha_accession}
 Neuraminidase accession: ${params.na_accession}
 Base output directory: ${params.outdir}
 Experiment ID: ${params.experiment_id}
 Experiment output directory: ${params.experiment_output}
 Run molecular dynamics: ${params.run_md}
 NCBI API key provided: ${params.ncbi_api_key ? 'Yes' : 'No'}
"""

// Main workflow
workflow {
    // Step 1: Retrieve sequences
    retrieveHA(params.ha_accession, 'hemagglutinin')
    ha_seq = retrieveHA.out.fasta

    retrieveNA(params.na_accession, 'neuraminidase')
    na_seq = retrieveNA.out.fasta
    
    // Step 2: Predict epitopes for both proteins
    ha_bcell = predictHABCellEpitopes(ha_seq, 'hemagglutinin')
    ha_tcell_i = predictHATCellEpitopesI(ha_seq, 'hemagglutinin')
    ha_tcell_ii = predictHATCellEpitopesII(ha_seq, 'hemagglutinin')
    
    na_bcell = predictNABCellEpitopes(na_seq, 'neuraminidase')
    na_tcell_i = predictNATCellEpitopesI(na_seq, 'neuraminidase')
    na_tcell_ii = predictNATCellEpitopesII(na_seq, 'neuraminidase')
    
    // Step 3: Combine epitopes and design vaccine for each protein
    ha_combined = combineHAEpitopes(
        ha_bcell.bcell_epitopes,
        ha_tcell_i.tcell_i_epitopes,
        ha_tcell_ii.tcell_ii_epitopes,
        'hemagglutinin'
    )
    
    na_combined = combineNAEpitopes(
        na_bcell.bcell_epitopes,
        na_tcell_i.tcell_i_epitopes,
        na_tcell_ii.tcell_ii_epitopes,
        'neuraminidase'
    )
    
    ha_vaccine = designHAVaccine(ha_combined.combined_epitopes, 'hemagglutinin')
    na_vaccine = designNAVaccine(na_combined.combined_epitopes, 'neuraminidase')
    
    // Optional: Combine HA and NA epitopes for a multi-target vaccine
    all_epitopes = Channel.empty()
    all_epitopes = all_epitopes.mix(ha_combined.combined_epitopes, na_combined.combined_epitopes)
    combined_vaccine = designCombinedVaccine(all_epitopes.collect(), 'combined')
    
    // Step 4: Evaluate vaccines
    ha_eval = evaluateHAVaccine(ha_vaccine.vaccine, 'hemagglutinin')
    na_eval = evaluateNAVaccine(na_vaccine.vaccine, 'neuraminidase')
    combined_eval = evaluateCombinedVaccine(combined_vaccine.vaccine, 'combined')
    
    // Step 5: Molecular dynamics (optional)
    if (params.run_md) {
        ha_md = runHAMD(ha_vaccine.vaccine, 'hemagglutinin')
        na_md = runNAMD(na_vaccine.vaccine, 'neuraminidase')
        combined_md = runCombinedMD(combined_vaccine.vaccine, 'combined')
    }
}

// On completion
workflow.onComplete {
    // Check if results were generated
    def outputDir = file(params.experiment_output)
    def fileCount = outputDir.listFiles()?.size() ?: 0
    
    log.info"""
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Experiment  : ${params.experiment_id}
    Output dir  : ${params.experiment_output}
    Files generated: ${fileCount}
    """
    
    // Print a warning if no files were generated
    if (workflow.success && fileCount == 0) {
        log.warn "WARNING: The pipeline completed successfully but no files were found in the output directory: ${params.experiment_output}"
        log.warn "Check that your module processes have publishDir directives with the correct path."
    }
}

workflow.onError {
    log.error "Pipeline execution stopped with error: ${workflow.errorMessage}"
    log.info "You can resume the pipeline with: nextflow run main.nf -resume"
}