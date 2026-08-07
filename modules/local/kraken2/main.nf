process kraken2 {

    tag { meta.id }

    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/kraken2", pattern: "${meta.id}_kraken2.txt", mode: 'copy'

    input:
    tuple val(meta), path(reads_1), path(reads_2), path(kraken2_db)

    output:
    tuple val(meta), path("${meta.id}_kraken2.txt")

    script:
    // After running kraken2, check if output file is empty. If it is, return exit code 1 to cause the process to error.
    """
    kraken2 --db ${kraken2_db} --threads 8 --output "-" --report ${meta.id}_kraken2.txt --paired ${reads_1} ${reads_2}
    [ -s ${meta.id}_kraken2.txt ]
    """
}
