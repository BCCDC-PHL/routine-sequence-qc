process interop_summary {
    /**
    * 
    * @input tuple val(run_id), path(run_dir)
    * @output path(interop_summary_json)
    */

    tag { run_id }

    executor 'local'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "interop*.{csv,json}", mode: 'copy'

    cpus 1

    input:
      tuple val(run_id), path(run_dir)

    output:
      tuple val(run_id), path("interop_summary.csv"), path("interop_summary.json"), path("interop_index-summary.csv")

    script:
      """
      interop_summary ${run_dir} --csv=1 > interop_summary.csv
      interop_index-summary ${run_dir} --csv=1 > interop_index-summary.csv
      parse_run_summary.py --summary interop_summary.csv > interop_summary.json
      """
}
