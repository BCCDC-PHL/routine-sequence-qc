process kraken2 {

    tag { sample_id }

    errorStrategy 'ignore'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sample_id}_kraken2.txt", mode: 'copy'


    input:
      tuple val(sample_id), path(reads_1), path(reads_2), path(kraken2_db)

    output:
      tuple val(sample_id), path("${sample_id}_kraken2.txt")

    script:
    // After running kraken2, check if output file is empty. If it is, return exit code 1 to cause the process to error.
    """
    kraken2 --db ${kraken2_db} --threads 8 --output "-" --report ${sample_id}_kraken2.txt --paired ${reads_1} ${reads_2}
    [ -s ${sample_id}_kraken2.txt ]
    """
}
