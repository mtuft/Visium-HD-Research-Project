#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * CRC Visium HD Spatial Clustering Pipeline
 * Processes P1, P2, P5 samples independently through Phase 1 and BANKSY clustering
 */

// Parameters
params.data_dir = "/scratch/prj/mmg_grp_single_cell_lab/projects/michael/data/CRC_VisiumHD"
params.outdir = "/scratch/prj/mmg_grp_single_cell_lab/projects/michael/pipeline_output/CRC_Nextflow"
params.container = "/scratch/prj/mmg_grp_single_cell_lab/projects/michael/bioc_3.19.sif"
params.r_libs = "/scratch/prj/mmg_grp_single_cell_lab/projects/michael/R_libs_4.4"

// Sample definitions
params.samples = ['P1', 'P2', 'P5']

// Create sample channel
Channel
    .from(params.samples)
    .map { sample_id ->
        def data_path = file("${params.data_dir}/${sample_id}/binned_outputs/square_008um")
        def deconv_path = file("${params.data_dir}/${sample_id}/DeconvolutionResults_${sample_id}CRC.csv.gz")
        tuple(sample_id, data_path, deconv_path)
    }
    .set { samples_ch }

/*
 * Process: Phase 1 Non-Spatial Clustering
 */
process PHASE1_CLUSTER {
    tag "${sample_id}"
    publishDir "${params.outdir}/Phase1/${sample_id}", mode: 'copy'
    
    cpus 4
    memory '256 GB'
    time '2h'
    
    input:
    tuple val(sample_id), path(data_dir), path(deconv_file)
    
    output:
    tuple val(sample_id), path("${sample_id}_phase1.rds"), path("plots/*"), emit: results
    path "phase1_${sample_id}.log", emit: log
    
    script:
    """
    mkdir -p plots
    
    export SINGULARITYENV_R_LIBS_USER="${params.r_libs}"
    
    singularity exec \
        --bind ${params.data_dir}:${params.data_dir} \
        --bind ${params.outdir}:${params.outdir} \
        --bind ${params.r_libs}:${params.r_libs} \
        --bind /cephfs:/cephfs \
        ${params.container} \
        Rscript ${projectDir}/scripts/phase1_cluster.R \
            ${data_dir} \
            ${deconv_file} \
            ${sample_id} \
            ${sample_id}_phase1.rds \
            plots \
        > phase1_${sample_id}.log 2>&1
    """
}

/*
 * Process: BANKSY Spatial Clustering
 */
process BANKSY_CLUSTER {
    tag "${sample_id}"
    publishDir "${params.outdir}/BANKSY/${sample_id}", mode: 'copy'
    
    cpus 4
    memory '256 GB'
    time '2h'
    
    input:
    tuple val(sample_id), path(data_dir), path(deconv_file)
    
    output:
    tuple val(sample_id), path("${sample_id}_banksy.rds"), path("plots/*"), emit: results
    path "banksy_${sample_id}.log", emit: log
    
    script:
    """
    mkdir -p plots
    
    export SINGULARITYENV_R_LIBS_USER="${params.r_libs}"
    
    singularity exec \
        --bind ${params.data_dir}:${params.data_dir} \
        --bind ${params.outdir}:${params.outdir} \
        --bind ${params.r_libs}:${params.r_libs} \
        --bind /cephfs:/cephfs \
        ${params.container} \
        Rscript ${projectDir}/scripts/banksy_cluster.R \
            ${data_dir} \
            ${deconv_file} \
            ${sample_id} \
            ${sample_id}_banksy.rds \
            plots \
        > banksy_${sample_id}.log 2>&1
    """
}

/*
 * Workflow
 */
workflow {
    // Run both pipelines in parallel for each sample
    PHASE1_CLUSTER(samples_ch)
    BANKSY_CLUSTER(samples_ch)
    
    // Combine outputs for summary
    PHASE1_CLUSTER.out.results
        .join(BANKSY_CLUSTER.out.results)
        .set { combined_results }
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'SUCCESS' : 'FAILED' }"
    println "Results directory: ${params.outdir}"
}
