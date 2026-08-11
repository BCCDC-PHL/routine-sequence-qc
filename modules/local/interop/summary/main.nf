process interop_summary {

    tag { run_id }

    executor 'local'
    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/interop", pattern: "interop*.{csv,json}", mode: 'copy'

    cpus 1

    input:
    tuple val(run_id), path(run_dir)

    output:
    tuple val(run_id), path("interop_summary.csv"), path("interop_index-summary.csv")

    script:
    """
    interop_summary ${run_dir} --csv=1 > interop_summary.csv
    interop_index-summary ${run_dir} --csv=1 > interop_index-summary.csv
    """
}
