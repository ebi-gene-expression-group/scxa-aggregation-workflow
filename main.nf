#!/usr/bin/env nextflow

resultsRoot = params.resultsRoot
quantDir = params.quantDir
expressionLevel = params.level
expressionScaling = params.scaling

// Find results from all the quantification subdirectories

Channel
    .fromPath("$quantDir/*", type: 'dir')
    .set { QUANT_DIRS }

Channel
    .fromPath("$quantDir/*/*.gtf.gz", checkIfExists: true )
    .set { REFERENCE_GTFS }

// Look at the results dirs and work out what quantification methods and
// protocols have been used

process gather_results {
    
    executor 'local'
    
    input:
        file(quantDir) from QUANT_DIRS

    output:
        set file('protocol'), file('quantType'), file('quantResults') into ALL_RESULTS

    """
        cp -p $quantDir/protocol protocol

        if [ -e $quantDir/kallisto ]; then
            echo kallisto > quantType
            cp -p $quantDir/kallisto quantResults
        elif [ -e $quantDir/alevin ]; then
            echo alevin > quantType
            cp -p $quantDir/alevin quantResults
        else
            echo "cannot determine quantificaiton type from \$(pwd)" 1>&2
            exit 1
        fi
    """

}

// Convert file outputs to strings

ALL_RESULTS
    .map{ row-> tuple( row[0].text, row[1].text, row[3]) }        
    .set{ ALL_RESULTS_VALS }


// Move Kallisto and Alevin results to different channels

ALEVIN_RESULTS = Channel.create()
KALLISTO_RESULTS = Channel.create()

ALL_RESULTS_VALS.choice( KALLISTO_RESULTS, ALEVIN_RESULTS ) {a -> 
    a[2] == 'kallisto' ? 1 : 0
}
    
// Make a transcript-to-gene mapping from the GTF file

process transcript_to_gene {

    conda "${baseDir}/envs/bioconductor-rtracklayer.yml"
    
    cache 'lenient'

    memory { 5.GB * task.attempt }
    errorStrategy { task.attempt<=3 ? 'retry' : 'finish' } 
    
    input:
        file gtf from REFERENCE_GTFS
    output:
        file tx2gene into TRANSCRIPT_TO_GENE_MANY

    """
        transcriptToGene.R ${gtf} transcript_id gene_id tx2gene
    """
}

// Allowing for the possibility of multiple sub-experiment2 in future,
// so creating a joint GTF. But there's probably only 1....

process merge_transcript_to_gene {

    input:
        file('??/tx2gene') from TRANSCRIPT_TO_GENE_MANY

    output:
        file tx2gene into TRANSCRIPT_TO_GENE

    """
    cat \$(ls */tx2gene | head -n 1) | head -n 1 > tx2gene
    tail -q -n +2 */tx2gene | sort | uniq >> tx2gene
    """    
}

// Generate the sets of files for each Kallisto sub-directory

process find_kallisto_results {
    
    executor 'local'
    
    input:
        set val(protocol), val(quantType), file('kallisto') from KALLISTO_RESULTS

    output:
        set val(protocol), file("kallisto_results.txt") into KALLISTO_RESULT_SETS

    """
        dir=\$(readlink kallisto)
        ls kallisto/*/abundance.h5 | while read -r l; do
            echo \${dir}/\$l >> kallisto_results.txt
        done
    """
}

// Split each result set into smaller chunks

process chunk_kallisto {

    executor 'local'

    input:
        set val(protocol), file(kallistoResults) from KALLISTO_RESULT_SETS

    output: 
        file("$protocol/chunks/*") into KALLISTO_CHUNKS

    """
        mkdir -p $protocol/chunks
        split -l ${params.chunkSize} ${kallistoResults} $protocol/chunks/
    """

}

// Flatten the chunk list

KALLISTO_CHUNKS
    .collect()
    .flatten()
    .set { FLATTENED_KALLISTO_CHUNKS }

// Re-associate the chunks with their protocols of origin

process associate_kallisto_chunks {

    executor 'local'
    
    input:
        file(kallistoChunk) from FLATTENED_KALLISTO_CHUNKS

    output:
        set stdout, file('out/kallisto_chunk') into READY_KALLISTO_CHUNKS 

    """
        mkdir -p out
        cp -p kallistoChunk out
        
        readlink \$kallistoChunk | awk -F'/' '{print \$(NF-2)}' 
    """
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
        set val(protocol), file(kallistoChunk) from READY_KALLISTO_CHUNKS        

    output:
        set val(protocol), file("counts_mtx") into KALLISTO_CHUNK_COUNT_MATRICES
        set val(protocol), file("tpm_mtx") into KALLISTO_CHUNK_ABUNDANCE_MATRICES
        set val(protocol), file("stats.tsv") into KALLISTO_CHUNK_STATS

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
            --outputCountsFile=counts_mtx/matrix.mtx \
            --outputAbundancesFile=tpm_mtx/matrix.mtx \
            --outputStatsFile=stats.tsv
        """
}

// Convert Alevin output to MTX. There will be one of these for every run, or
// technical replicate group of runs

process alevin_to_mtx {

    conda "${baseDir}/envs/parse_alevin.yml"

    input:
        set val(protocol), val(quantType), file('alevin') from ALEVIN_RESULTS

    output:
        set val(protocol), file("counts_mtx") into ALEVIN_CHUNK_COUNT_MATRICES

    """
    dir=\$(readlink alevin)
    alevinToMtx.py alevin counts_mtx
    """ 
} 

// Merge the chunks for each protocol into one matrix. For Kallisto
// results this will be sub-matrices generated due the costs of running
// tximport on 10s of 1000s of runs. For Alevin this will be the matrices
// generated for each library

KALLISTO_CHUNK_COUNT_MATRICES
    .concat(ALEVIN_CHUNK_COUNT_MATRICES)
    .groupTuple()
    .set { PROTOCOL_COUNT_CHUNKS }

KALLISTO_CHUNK_ABUNDANCE_MATRICES
    .groupTuple()
    .set { PROTOCOL_KALLISTO_ABUNDANCE_CHUNKS }

process merge_count_chunk_matrices {
    
    conda "${baseDir}/envs/kallisto_matrix.yml"

    cache 'lenient'
    
    memory { 5.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20
    
    input:
        set val(protocol), file('dir??/*') from PROTOCOL_COUNT_CHUNKS

    output:
        file("counts_mtx_${protocol}") into PROTOCOL_COUNT_MATRICES

    """
        find . -name 'counts_mtx' > dirs.txt
        mergeMtx.R dirs.txt counts_mtx_${protocol}
        rm -f dirs.txt
    """
}

// Merge the sub-experiments corresponding to different protocols

process merge_protocol_count_matrices {
    
    conda "${baseDir}/envs/kallisto_matrix.yml"

    cache 'lenient'
    
    memory { 5.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20
    
    publishDir "$resultsRoot/matrices", mode: 'copy', overwrite: true
    
    input:
        file('*') from PROTOCOL_COUNT_MATRICES.collect()

    output:
        file("counts_mtx.zip") into EXP_COUNT_MATRICES

    """
        find . -name 'counts_mtx_*' > dirs.txt
        mergeMtx.R dirs.txt counts_mtx
        rm -f dirs.txt
        zip -r counts_mtx.zip counts_mtx
    """
}

process merge_tpm_chunk_matrices {

    conda "${baseDir}/envs/kallisto_matrix.yml"
    
    memory { 5.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20
    
    publishDir "$resultsRoot/matrices", mode: 'copy', overwrite: true
    
    input:
        set val(protocol), file('dir??/*') from PROTOCOL_KALLISTO_ABUNDANCE_CHUNKS

    output:
        set val(protocol), file("tpm_mtx.zip")

    """
        find . -name 'tpm_mtx' > dirs.txt
        mergeMtx.R dirs.txt tpm_mtx
        rm -f dirs.txt
        zip -r tpm_mtx.zip tpm_mtx
    """
}

KALLISTO_CHUNK_STATS
    .collectFile( sort: true, name: "stats.tsv", storeDir: "${resultsRoot}/matrices", keepHeader: true )

