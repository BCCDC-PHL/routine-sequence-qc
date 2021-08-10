process seqkit_stats {

    tag { sample_id }

    // errorStrategy 'ignore'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sample_id}_seqkit_stats.csv", mode: 'copy'


    input:
      tuple val(grouping_key), path(reads)

    output:
      tuple val(sample_id), path("${sample_id}_seqkit_stats.csv")

    script:
    if (grouping_key =~ '_S[0-9]+_') {
        sample_id = grouping_key.split("_S[0-9]+_")[0]
    } else if (grouping_key =~ '_') {
      sample_id = grouping_key.split("_")[0]
    } else {
      sample_id = grouping_key
    }
    """
    echo 'file,format,type,num_seqs,sum_len,min_len,avg_len,max_len,q1_len,q2_len,q3_len,sum_gap,n50,percent_bases_greater_than_q20,percent_bases_greater_than_q30' > header.csv
    seqkit stats --all -j ${task.cpus} --tabular ${reads} | tr \$'\\t' ',' > seqkit_stats.csv
    cat header.csv <(tail -n+2 seqkit_stats.csv) > ${sample_id}_seqkit_stats.csv
    """
}
