#!/usr/bin/env nextflow

resultsRoot = params.resultsRoot
quantDir = params.quantDir
expressionLevel = params.level
expressionScaling = params.scaling
referenceGtf = params.referenceGtf

REFERENCE_GTF = Channel.fromPath("${resultsRoot}/${referenceGtf}", checkIfExists: true)

// Make a transcript-to-gene mapping from the GTF file

process transcript_to_gene {

    conda "${baseDir}/envs/bioconductor-rtracklayer.yml"
    
    cache 'lenient'

    memory { 5.GB * task.attempt }
    errorStrategy { task.attempt<=3 ? 'retry' : 'finish' } 
    
    input:
        file gtf from REFERENCE_GTF
    output:
        file tx2gene into TRANSCRIPT_TO_GENE

    """
        transcriptToGene.R ${gtf} transcript_id gene_id tx2gene
    """
}

// Load Kallisto outputs into a channel, split in to chunks of the size
// specified in the parameters

Channel.fromPath( "${quantDir}/*/abundance.h5" )
    .map { "${it}" }
    .collectFile(sort: true, name: 'kallisto_results.txt', newLine: true)
    .splitText( by: params.chunkSize )
    .set{
        KALLISTO_CHUNKS
    } 

// Note: we can call tximport in different ways to create different matrix types 

process kallisto_gene_count_matrix {
    
    conda "${baseDir}/envs/kallisto_matrix.yml"

    cache 'deep'

    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20

    input:
        file tx2Gene from TRANSCRIPT_TO_GENE.first()
        file(kallistoChunk) from KALLISTO_CHUNKS        

    output:
        file("${expressionLevel}_${expressionScaling}_counts") into KALLISTO_CHUNK_COUNT_MATRICES
        file("${expressionLevel}_${expressionScaling}_tpm") into KALLISTO_CHUNK_ABUNDANCE_MATRICES
        file("${expressionLevel}_${expressionScaling}.stats.tsv") into KALLISTO_CHUNK_STATS

    script:

        def txOut
        if ( expressionLevel == 'transcript' ){
            txOut = 'TRUE'
        }else{
            txOut = 'FALSE'
        }

        """
        tximport.R --files=${kallistoChunk} --type=kallisto --tx2gene=$tx2Gene \
            --countsFromAbundance=$expressionScaling --ignoreTxVersion=TRUE --txOut=$txOut \
            --outputCountsFile=${expressionLevel}_${expressionScaling}_counts/matrix.mtx \
            --outputAbundancesFile=${expressionLevel}_${expressionScaling}_tpm/matrix.mtx \
            --outputStatsFile=${expressionLevel}_${expressionScaling}.stats.tsv
        """
}

// Merge the chunks into one matrix

process merge_count_matrices {
    
    conda "${baseDir}/envs/kallisto_matrix.yml"

    cache 'lenient'
    
    memory { 5.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20
    
    publishDir "$resultsRoot/matrices", mode: 'copy', overwrite: true
    
    input:
        file('dir??/*') from KALLISTO_CHUNK_COUNT_MATRICES.collect()

    output:
        file("${expressionLevel}_${expressionScaling}_counts.zip")

    """
        find . -name '${expressionLevel}_${expressionScaling}_counts' > dirs.txt
        mergeMtx.R dirs.txt ${expressionLevel}_${expressionScaling}_counts
        rm -f kallisto_results.txt
        zip -r ${expressionLevel}_${expressionScaling}_counts.zip ${expressionLevel}_${expressionScaling}_counts
    """
}

process merge_tpm_matrices {

    conda "${baseDir}/envs/kallisto_matrix.yml"
    
    memory { 5.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20
    
    publishDir "$resultsRoot/matrices", mode: 'copy', overwrite: true
    
    input:
        file('dir??/*') from KALLISTO_CHUNK_ABUNDANCE_MATRICES.collect()

    output:
        file("${expressionLevel}_${expressionScaling}_tpm.zip")

    """
        find . -name '${expressionLevel}_${expressionScaling}_tpm' > dirs.txt
        mergeMtx.R dirs.txt ${expressionLevel}_${expressionScaling}_tpm
        rm -f kallisto_results.txt
        zip -r ${expressionLevel}_${expressionScaling}_tpm.zip ${expressionLevel}_${expressionScaling}_tpm
    """
}

KALLISTO_CHUNK_STATS
    .collectFile( sort: true, name: "${expressionLevel}_${expressionScaling}.stats.tsv", storeDir: "${resultsRoot}/matrices", keepHeader: true )

