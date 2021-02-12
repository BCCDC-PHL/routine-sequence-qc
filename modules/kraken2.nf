process kraken2 {

    tag { sample_id }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sample_id}_kraken2.txt", mode: 'copy'

    cpus 8

    input:
      tuple val(grouping_key), path(reads), path(kraken2_db)

    output:
      tuple val(sample_id), path("${sample_id}_kraken2.txt")

    script:
    if (grouping_key =~ '_S[0-9]+_') {
        sample_id = grouping_key.split("_S[0-9]+_")[0]
    } else if (grouping_key =~ '_') {
      sample_id = grouping_key.split("_")[0]
    } else {
      sample_id = grouping_key
    }
    """
    kraken2 --db ${kraken2_db} --threads 8 --output "-" --report ${sample_id}_kraken2.txt --paired ${reads}
    """
}
