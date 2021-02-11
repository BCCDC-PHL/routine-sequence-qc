process multiqc {
    /**
    * 
    * @input tuple val(sample_id), path(forward), path(reverse)
    * @output 
    */

    tag { run_id }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "multiqc_*", mode: 'copy'

    cpus 2

    input:
      path("fastqc_outdir")

    output:
      tuple path("multiqc_report.html"), path("multiqc_data")

    script:
      """
      multiqc .
      """
}
