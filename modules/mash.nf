process mash {
    /**
    * 
    * @input tuple val(sample_id), path(forward)
    * @output "${params.outdir}/mash/${sample_id}.msh.results"
    */

    tag { sample_id }


    publishDir "${params.outdir}/mash", mode: 'copy'

    cpus 4

    input:
      tuple val(grouping_key), path(reads)

    output:
        path ("${sample_id}.msh.results")

    script:
      if (grouping_key =~ '_S[0-9]+_') {
        sample_id = grouping_key.split("_S[0-9]+_")[0]
      } else if (grouping_key =~ '_') {
        sample_id = grouping_key.split("_")[0]
      } else {
        sample_id = grouping_key
      }
      """
        mash sketch -r ${reads[0]} -m 5 -k 21 -C ${sample_id} -o ${sample_id} &> ${sample_id}.msh.results
      """
}

process mashTable {

    publishDir "${params.outdir}/mash", mode: 'copy'

    input:
        file mash_outputs

    output:
        file 'mash_summary.csv'

    script:
    """
    echo "sample,estGenomeSize,estCoverage" > mash_summary.csv
    for mash_file in *.results; 
        do awk 'BEGIN{OFMT="%f"; print ARGV[1]} {if(NR==1) {print ","\$4} if(NR==2) {print ","\$3*2}}' \$mash_file | tr -d '\\n' | sed 's/.msh.results//' >> mash_summary.csv;
        echo "" >> mash_summary.csv
    done 
    """

} 
