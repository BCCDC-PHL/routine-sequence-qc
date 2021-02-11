process fastqc {
    /**
    * 
    * @input tuple val(sample_id), path(forward), path(reverse)
    * @output 
    */

    tag { sample_id }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sample_id}_*_fastqc", mode: 'copy'

    cpus 2

    input:
      tuple val(grouping_key), path(reads)

    output:
      tuple val(sample_id), path("${sample_id}_R1_fastqc"), path("${sample_id}_R2_fastqc")

    script:
      if (grouping_key =~ '_S[0-9]+_') {
        sample_id = grouping_key.split("_S[0-9]+_")[0]
      } else if (grouping_key =~ '_') {
        sample_id = grouping_key.split("_")[0]
      } else {
        sample_id = grouping_key
      }
      """
      fastqc -t 2 ${reads}
      for d in *.zip; do unzip \$d; done
      mv ${sample_id}_*_R1_*_fastqc ${sample_id}_R1_fastqc
      mv ${sample_id}_*_R2_*_fastqc ${sample_id}_R2_fastqc
      """
}
