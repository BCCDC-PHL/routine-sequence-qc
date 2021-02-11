process interop_summary {
    /**
    * 
    * @input tuple val(run_id), path(run_dir)
    * @output path(interop_summary_json)
    */

    tag { run_id }

    executor 'local'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "interop_summary.json", mode: 'copy'

    cpus 1

    input:
      tuple val(run_id), path(run_dir)

    output:
      tuple path("interop_summary.json")

    script:
      """
      interop_summary ${run_dir} > interop_summary.txt
      parse_run_summary.py --summary interop_summary.txt | python -m json.tool > interop_summary.json
      """
}
