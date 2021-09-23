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
      tuple val(sample_id), path(reads_1), path(reads_2)

    output:
      tuple val(sample_id), path("${sample_id}_R1_fastqc"), path("${sample_id}_R2_fastqc")

    script:
      """
      fastqc -t 2 ${reads_1} ${reads_2}
      for d in *.zip; do unzip \$d; done
      if [ ! -d "${sample_id}_R1_fastqc" ]
      then
        mv ${sample_id}*_R1*_fastqc ${sample_id}_R1_fastqc
      fi
      if [ ! -d "${sample_id}_R2_fastqc" ]
      then
        mv ${sample_id}*_R2*_fastqc ${sample_id}_R2_fastqc
      fi
      """
}
