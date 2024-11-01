#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.resultsRoot = ''
params.quantDir = ''
params.level = ''
params.scaling = ''
params.chunkSize = ''
params.reference = [ignoreTxVersion: '']

// Process definitions

process gather_results {
    executor 'local'
    
    input:
    path quantDir

    output:
    tuple path('protocol'), path('quantType'), path('quantResults')

    script:
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

process merge_transcript_to_gene {
    input:
    path('??/tx2gene')

    output:
    path 'tx2gene'

    """
    cat \$(ls */tx2gene | head -n 1) | head -n 1 > tx2gene
    tail -q -n +2 */tx2gene | sort | uniq >> tx2gene
    """    
}

process find_kallisto_results {
    executor 'local'
    
    input:
    tuple val(protocol), val(quantType), path('kallisto')

    output:
    tuple val(protocol), path("kallisto_results.txt")

    """
    dir=\$(readlink kallisto)
    ls kallisto/*/abundance.h5 | while read -r l; do
        echo \$(dirname \${dir})/\$l >> kallisto_results.txt
    done
    """
}

process chunk_kallisto {
    executor 'local'

    input:
    tuple val(protocol), path(kallistoResults)

    output: 
    tuple val(protocol), path("chunks/*")

    """
    mkdir -p chunks
    split -l ${params.chunkSize} ${kallistoResults} chunks/
    """
}

process kallisto_gene_count_matrix {
    conda "${baseDir}/envs/kallisto_matrix.yml"

    cache 'deep'

    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20

    input:
    path tx2Gene
    tuple val(protocol), path(kallistoChunk)        

    output:
    tuple val(protocol), path("counts_mtx"), emit: counts
    tuple val(protocol), path("tpm_mtx"), emit: tpm
    path "kallisto_stats.tsv", emit: stats

    script:
    def txOut = (params.level == 'transcript') ? 'TRUE' : 'FALSE'
    """
    ignoreTxVersion=${params.reference.ignoreTxVersion}
    example_file=\$(head -n 1 ${kallistoChunk})
    example_id=\$(sed '2q;d'  \${example_file/\\.h5/.tsv} | awk '{print \$1}')
    grep -P "^\$example_id\t" tx2gene > /dev/null

    if [ \$? -eq 0 ]; then
        ignoreTxVersion=FALSE
    fi

    sed -e 's/\t/,/g' ${tx2Gene} > ${tx2Gene}.csv
    tximport.R --files=${kallistoChunk} --type=kallisto --tx2gene=${tx2Gene}.csv \
        --countsFromAbundance=${params.scaling} --ignoreTxVersion=\$ignoreTxVersion --txOut=$txOut \
        --outputCountsFile=counts_mtx/matrix.mtx \
        --outputAbundancesFile=tpm_mtx/matrix.mtx \
        --outputStatsFile=kallisto_stats.tsv
    """
}

process alevin_runs {
    executor 'local'
    
    input:
    tuple val(protocol), val(quantType), path('alevin')

    output:
    tuple val(protocol), path("alevin_runs/*")
    
    """
    cp -P alevin alevin_runs
    """
}

process alevin_to_mtx {
    conda "${baseDir}/envs/parse_alevin.yml"
    
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 || task.exitStatus == 141 ? 'retry' : 'finish' }
    maxRetries 10

    input:
    tuple val(protocol), path('alevin_run')

    output:
    tuple val(protocol), path("counts_mtx")

    """
    ln -s alevin_run/alevin/mtx/counts_mtx_nonempty counts_mtx
    """ 
}

process alevin_stats {
    conda 'r-rjson'

    input:
    tuple val(protocol), path('alevin_run')

    output:
    tuple val(protocol), path("alevin_stats.tsv")

    """
    #!/usr/bin/env Rscript
    
    suppressPackageStartupMessages(library(rjson))    
    
    json <- fromJSON(file = "alevin_run/meta_info.json") 
    stats <- t(data.frame(unlist(lapply(json, function(j) paste(j, collapse = ' ')))))
    run <- basename(Sys.readlink("alevin_run"))
    
    write.table(data.frame(cbind(run=run, stats)), file = 'alevin_stats.tsv', quote = FALSE, sep="\\t", row.names=FALSE)
    """
}

process merge_count_chunk_matrices {
    conda "${baseDir}/envs/kallisto_matrix.yml"

    cache 'lenient'
    
    memory { 5.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20
    
    input:
    tuple val(protocol), path('dir??/*')

    output:
    path "counts_mtx_${protocol}"

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

process merge_protocol_count_matrices {
    conda "${baseDir}/envs/kallisto_matrix.yml"

    cache 'lenient'
    
    memory { 5.GB * task.attempt }
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 20
    
    publishDir "${params.resultsRoot}/matrices", mode: 'copy', overwrite: true
    
    input:
    path '*'

    output:
    path "counts_mtx.zip"

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
    
    publishDir "${params.resultsRoot}/matrices", mode: 'copy', overwrite: true
    
    input:
    tuple val(protocol), path('dir??/*')

    output:
    tuple val(protocol), path("tpm_mtx.zip")

    """
    find . -name 'tpm_mtx' > dirs.txt
    mergeMtx.R dirs.txt tpm_mtx
    rm -f dirs.txt
    zip -r tpm_mtx.zip tpm_mtx
    """
}

// Workflow definition

workflow {
    // Input channels
    quant_dirs_ch = Channel.fromPath("${params.quantDir}/*", type: 'dir')
    transcript_to_gene_ch = Channel.fromPath("${params.quantDir}/*/transcript_to_gene.txt", checkIfExists: true)

    // Process execution
    gather_results(quant_dirs_ch)
    merge_transcript_to_gene(transcript_to_gene_ch.collect())

    // Split results into Kallisto and Alevin
    gather_results.out
        .map { it -> [it[0].text, it[1].text, it[2]] }
        .branch {
            kallisto: it[1] == 'kallisto'
            alevin: it[1] == 'alevin'
        }
        .set { all_results }

    // Kallisto workflow
    find_kallisto_results(all_results.kallisto)
    chunk_kallisto(find_kallisto_results.out)
    kallisto_gene_count_matrix(merge_transcript_to_gene.out, chunk_kallisto.out.transpose())

    // Alevin workflow
    alevin_runs(all_results.alevin)
    alevin_to_mtx(alevin_runs.out.transpose())
    alevin_stats(alevin_runs.out.transpose())

    // Merge count matrices
    merge_count_chunk_matrices(
        kallisto_gene_count_matrix.out.counts.mix(alevin_to_mtx.out).groupTuple()
    )
    merge_protocol_count_matrices(merge_count_chunk_matrices.out.collect())

    // Merge TPM matrices (Kallisto only)
    merge_tpm_chunk_matrices(kallisto_gene_count_matrix.out.tpm.groupTuple())

    // Collect stats
    kallisto_gene_count_matrix.out.stats
        .collectFile(name: "kallisto_stats.tsv", storeDir: "${params.resultsRoot}/matrices", keepHeader: true)
    
    alevin_stats.out
        .collectFile(name: "alevin_stats.tsv", storeDir: "${params.resultsRoot}/matrices", keepHeader: true)
}
