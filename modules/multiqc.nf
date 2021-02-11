process multiqc {
    /**
    * 
    * @input tuple val(run_id), [path(fastqc_outdir)]
    * @output tuple path(multiqc_report), path(multiqc_data)
    */

    tag { run_id }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "multiqc_*", mode: 'copy'

    cpus 2

    input:
      tuple val(run_id), path("fastqc_outdir")

    output:
      tuple path("multiqc_report.html"), path("multiqc_data")

    script:
      """
      multiqc .
      """
}
