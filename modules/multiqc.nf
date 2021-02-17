process multiqc {
    /**
    * 
    * @input tuple val(run_id), [path(qc_outdir)]
    * @output tuple path(multiqc_report), path(multiqc_data)
    */

    tag { run_id }
    
    validExitStatus 0,1

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "multiqc_*", mode: 'copy'

    cpus 2

    input:
      tuple path(multiqc_config), val(run_id), path(qc_outdir)

    output:
      tuple path("multiqc_report.html"), path("multiqc_data")

    script:
      """
      cp ${multiqc_config} multiqc_config.yaml
      echo "report_header_info:" >> multiqc_config.yaml
      echo "  - Run ID: ${run_id}" >> multiqc_config.yaml
      multiqc .
      """
}
