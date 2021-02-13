process multiqc {
    /**
    * 
    * @input tuple val(run_id), [path(qc_outdir)]
    * @output tuple path(multiqc_report), path(multiqc_data)
    */

    tag { run_id }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${run_id}_multiqc_*", mode: 'copy'

    cpus 2

    input:
      tuple val(run_id), path(qc_outdir)

    output:
      tuple path("${run_id}_multiqc_report.html"), path("${run_id}_multiqc_report_data")

    script:
      """
      multiqc --title "${run_id}" .
      """
}
