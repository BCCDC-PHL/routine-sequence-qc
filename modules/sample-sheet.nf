process parse_sample_sheet {
    /**
    * 
    * @input tuple val(run_id), path(sample_sheet_csv)
    * @output path(sample_sheet_json)
    */

    tag { run_id }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "sample_sheet.json", mode: 'copy'

    executor 'local'

    cpus 1

    input:
      tuple val(run_id), path(sample_sheet)

    output:
      path("sample_sheet.json")

    script:
    def parser_script = params.instrument_type == "miseq" ? "samplesheet_parser_miseq.py" : "samplesheet_parser_nextseq.py"
      """
      ${parser_script} ${sample_sheet} | python -m json.tool > sample_sheet.json
      """
}
