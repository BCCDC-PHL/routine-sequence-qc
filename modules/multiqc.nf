process multiqc {
    /**
    * 
    * @input tuple val(run_id), [path(qc_outdir)]
    * @output tuple path(multiqc_report), path(multiqc_data)
    */

    errorStrategy 'finish'

    tag { run_id }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "multiqc_*", mode: 'copy'
    publishDir "${params.outdir}", pattern: "pipeline_complete.json", mode: 'copy'

    cpus 2

    input:
      tuple path(multiqc_config), val(run_id), path(qc_outdir)

    output:
      tuple path("multiqc_report.html"), path("multiqc_data"), path("pipeline_complete.json")

    script:
      """
      cp ${multiqc_config} multiqc_config.yaml
      echo "report_header_info:" >> multiqc_config.yaml
      echo "  - Run ID: ${run_id}" >> multiqc_config.yaml
      multiqc .
      echo "{\\"pipeline_name\\": \\"BCCDC-PHL/routine-sequence-qc\\", \\"timestamp_completed\\": \\"\$(date --iso=seconds)\\"}" > pipeline_complete.json
      """
}
