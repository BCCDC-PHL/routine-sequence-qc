process mash_sketch {

    tag { output_basename }

    errorStrategy 'ignore'
    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/mash_sketch", pattern: "${output_basename}_mash_sketch.txt", mode: 'copy'

    input:
    tuple val(meta), path(reads_1), path(reads_2)

    output:
    tuple val(meta), path ("${output_basename}_mash_sketch.txt")

    script:
    output_basename = reads_1.baseName.split('\\.')[0]
    """
    mash sketch \
      -r ${reads_1} \
      -m ${params.mash_sketch_minimum_copies} \
      -k ${params.mash_sketch_kmer_size} \
      -C ${output_basename} \
      -o ${output_basename} \
      &> ${output_basename}_mash_sketch.txt \
      || echo -e "Estimated genome size: 0.0\nEstimated coverage:    0.0" \
      > ${output_basename}_mash_sketch.txt
    """
}
