#!/usr/bin/env nextflow

gtf = params.reference.gtf
resultsRoot = params.resultsRoot

expressionLevel = params.level
expressionScaling = params.scaling

REFERENCE_GTF = Channel.fromPath(gtf, checkIfExists: true)

// Make a transcript-to-gene mapping from the GTF file

process transcript_to_gene {

    conda 'bioconductor-rtracklayer'
    
    memory { 5.GB * task.attempt }
    errorStrategy { task.attempt<=3 ? 'retry' : 'finish' } 
    
    input:
        file gtf from REFERENCE_GTF
    output:
        file tx2gene into TRANSCRIPT_TO_GENE

    """
        $SCRIPTS_DIR/transcriptToGene.R ${gtf} transcript_id gene_id tx2gene
    """
}

// Note: we can call tximport in different ways to create different matrix types 

process kallisto_gene_count_matrix {
    
    conda "${baseDir}/envs/kallisto_matrix.yml"

    memory { 5.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'finish' }
    maxRetries 20

    publishDir "$resultsRoot/matrices", mode: 'move', overwrite: true
    
    input:
        file tx2Gene from TRANSCRIPT_TO_GENE

    output:
        file("${expressionLevel}_${expressionScaling}_counts.zip") into KALLISTO_COUNT_MATRIX
        file("${expressionLevel}_${expressionScaling}_tpm.zip") into KALLISTO_ABUNDANCE_MATRIX
        file("${expressionLevel}_${expressionScaling}.stats.tsv") into KALLISTO_STATS

    script:

        def txOut
        if ( expressionLevel == 'transcript' ){
            txOut = 'TRUE'
        }else{
            txOut = 'FALSE'
        }

        """
        find -L ${resultsRoot} -name abundance.h5 > kallisto_results.txt
        tximport.R --files=kallisto_results.txt --type=kallisto --tx2gene=$tx2Gene \
            --countsFromAbundance=$expressionScaling --ignoreTxVersion=TRUE --txOut=$txOut \
            --outputCountsFile=${expressionLevel}_${expressionScaling}_counts/matrix.mtx \
            --outputAbundancesFile=${expressionLevel}_${expressionScaling}_tpm/matrix.mtx \
            --outputStatsFile=${expressionLevel}_${expressionScaling}.stats.tsv
        
        rm -f kallisto_results.txt
        zip -r ${expressionLevel}_${expressionScaling}_counts.zip ${expressionLevel}_${expressionScaling}_counts
        zip -r ${expressionLevel}_${expressionScaling}_tpm.zip ${expressionLevel}_${expressionScaling}_tpm
        """
}
