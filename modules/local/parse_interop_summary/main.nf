process parse_interop_summary {

    tag { run_id }

    executor 'local'
    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/interop", pattern: "interop_summary.json", mode: 'copy'

    cpus 1

    input:
    tuple val(run_id), path(interop_summary_csv)

    output:
    tuple val(run_id), path("interop_summary.json")

    script:
    """
    parse_interop_summary.py --summary interop_summary.csv > interop_summary.json
    """
}
