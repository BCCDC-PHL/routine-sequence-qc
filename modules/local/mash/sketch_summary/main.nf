process mash_sketch_summary {

    tag { meta.id }

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(mash_sketch_output)

    output:
    tuple val(meta), path("${meta.id}_mash_sketch_summary.csv")

    script:
    """
    summarize_mash_sketch.py ${mash_sketch_output} --sample-id ${meta.id} > ${meta.id}_mash_sketch_summary.csv
    """
}
