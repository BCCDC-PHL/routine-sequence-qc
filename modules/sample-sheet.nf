process parse_sample_sheet {

    tag { run_id }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "sample_sheet.json", mode: 'copy'

    executor 'local'

    cpus 1

    input:
      tuple val(run_id), path(sample_sheet)

    output:
      path("sample_sheet.json")

    script:
    def parser_script = "samplesheet_parser_nextseq.py"
    if (params.instrument_type == "miseq"){
        parser_script = "samplesheet_parser_miseq.py"
    } else if (params.instrument_type == "i100") {
        parser_script = "samplesheet_parser_i100.py"
    }
    """
    ${parser_script} ${sample_sheet} | python -m json.tool > sample_sheet.json
    """
}
