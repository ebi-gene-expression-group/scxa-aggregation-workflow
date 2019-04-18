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
            echo -n kallisto > quantType
            cp -rp $quantDir/kallisto quantResults
        elif [ -e $quantDir/alevin ]; then
            echo -n alevin > quantType
            cp -rp $quantDir/alevin quantResults
        else
            echo "cannot determine quantification type from \$(pwd)" 1>&2
            exit 1
        fi
    """

}

// Convert file outputs to strings

ALL_RESULTS
    .map{ row-> tuple( row[0].text, row[1].text, row[2]) }        
    .set{ ALL_RESULTS_VALS }


// Move Kallisto and Alevin results to different channels

ALEVIN_RESULTS = Channel.create()
KALLISTO_RESULTS = Channel.create()

ALL_RESULTS_VALS.choice( KALLISTO_RESULTS, ALEVIN_RESULTS ) {a -> 
    a[1] == 'kallisto' ? 0 : 1
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
            echo \$(dirname \${dir})/\$l >> kallisto_results.txt
        done
    """
}

// Split each result set into smaller chunks

process chunk_kallisto {

    executor 'local'

    input:
        set val(protocol), file(kallistoResults) from KALLISTO_RESULT_SETS

    output: 
        set val(protocol), file("chunks/*") into KALLISTO_CHUNKS

    """
        mkdir -p chunks
        split -l ${params.chunkSize} ${kallistoResults} chunks/
    """

}

// Flatten the chunk list

KALLISTO_CHUNKS
    .transpose()
    .set { FLATTENED_KALLISTO_CHUNKS }

// Note: we can call tximport in different ways to create different matrix types 

process kallisto_gene_count_matrix {
    
    conda "${baseDir}/envs/kallisto_matrix.yml"

    cache 'deep'

    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20

    input:
        file tx2Gene from TRANSCRIPT_TO_GENE.first()
        set val(protocol), file(kallistoChunk) from FLATTENED_KALLISTO_CHUNKS        

    output:
        set val(protocol), file("counts_mtx") into KALLISTO_CHUNK_COUNT_MATRICES
        set val(protocol), file("tpm_mtx") into KALLISTO_CHUNK_ABUNDANCE_MATRICES
        file("kallisto_stats.tsv") into KALLISTO_CHUNK_STATS

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
            --outputStatsFile=kallisto_stats.tsv
        """
}

// Get the run-wise Alevin results. In the case of Alevin, we'll have one
// matrix from each library. We can just copy the symlink to the 'alevin'
// folder that contains the library-wise Alevin runs. Nextflow will then put
// each result set into the output channel.

process alevin_runs {

    executor 'local'
    
    input:
        set val(protocol), val(quantType), file('alevin') from ALEVIN_RESULTS

    output:
         set val(protocol), file("alevin_runs/*") into ALEVIN_RESULTS_BY_LIB
    
    """
    cp -rp alevin alevin_runs
    """
}

// Flatten Alevin channel

ALEVIN_RESULTS_BY_LIB
    .transpose()
    .into{
        FLATTENED_ALEVIN_RESULTS_BY_LIB
        FLATTENED_ALEVIN_RESULTS_BY_LIB_FOR_STATS
    }

// Convert Alevin output to MTX. There will be one of these for every run, or
// technical replicate group of runs

process alevin_to_mtx {

    conda "${baseDir}/envs/parse_alevin.yml"
    
    memory { 10.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20

    input:
        set val(protocol), file('alevin_run') from FLATTENED_ALEVIN_RESULTS_BY_LIB

    output:
        set val(protocol), file("counts_mtx") into ALEVIN_CHUNK_COUNT_MATRICES

    """
    runId=\$(basename \$(readlink alevin_run))
    alevinToMtx.py --cell_prefix \${runId}- alevin_run counts_mtx
    """ 
} 
  
// Extract the output stats for each run and store to a tsv for later collation
 
process alevin_stats {
    
    conda 'r-rjson'

    input:
        set val(protocol), file('alevin_run') from FLATTENED_ALEVIN_RESULTS_BY_LIB_FOR_STATS

    output:
        set val(protocol), file("alevin_stats.tsv") into ALEVIN_CHUNK_STATS

    """
    #!/usr/bin/env Rscript
    
    suppressPackageStartupMessages(library(rjson))    
    
    json <- fromJSON(file = "alevin_run/aux_info/alevin_meta_info.json") 
    stats <- t(data.frame(unlist(lapply(json, function(j) paste(j, collapse = ' ')))))
    run <- basename(Sys.readlink("alevin_run"))
    
    write.table(data.frame(cbind(run=run, stats)), file = 'alevin_stats.tsv', quote = FALSE, sep="\\t", row.names=FALSE)
    """
    
}
 
// Merge the chunks for each protocol into one matrix. For Kallisto
// results this will be sub-matrices generated due the costs of running
// tximport on 10s of 1000s of runs. For Alevin this will be the matrices
// generated for each library

ALEVIN_CHUNK_COUNT_MATRICES
    .concat(KALLISTO_CHUNK_COUNT_MATRICES)
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
        find \$(pwd) -name 'counts_mtx' > dirs.txt
        ndirs=\$(cat dirs.txt | wc -l)
        if [ "\$ndirs" -gt 1 ]; then 
            mergeMtx.R dirs.txt counts_mtx_${protocol}
        else
            ln -s \$(cat dirs.txt) counts_mtx_${protocol}
        fi
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
        find \$(pwd) -name 'counts_mtx_*' > dirs.txt
        
        ndirs=\$(cat dirs.txt | wc -l)
        if [ "\$ndirs" -gt 1 ]; then 
            mergeMtx.R dirs.txt counts_mtx
        else
            ln -s \$(cat dirs.txt) counts_mtx
        fi
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
    .collectFile( sort: true, name: "kallisto_stats.tsv", storeDir: "${resultsRoot}/matrices", keepHeader: true )

ALEVIN_CHUNK_STATS
    .collectFile( sort: true, name: "alevin_stats.tsv", storeDir: "${resultsRoot}/matrices", keepHeader: true )
