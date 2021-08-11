process mash_sketch {

    tag { output_basename }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${output_basename}_mash_sketch.txt", mode: 'copy'

    input:
      tuple val(grouping_key), path(reads)

    output:
      tuple val(sample_id), path ("${output_basename}_mash_sketch.txt")

    script:
      if (grouping_key =~ '_S[0-9]+_') {
        sample_id = grouping_key.split("_S[0-9]+_")[0]
      } else if (grouping_key =~ '_') {
        sample_id = grouping_key.split("_")[0]
      } else {
        sample_id = grouping_key
      }
      output_basename = reads[0].baseName.split('\\.')[0]
      """
      mash sketch -r ${reads[0]} -m ${params.mash_sketch_minimum_copies} -k ${params.mash_sketch_kmer_size} -C ${output_basename} -o ${output_basename} &> ${output_basename}_mash_sketch.txt
      """
}

process mash_sketch_summary {

    tag { sample_id }

    cpus 1

    executor 'local'

    input:
      tuple val(sample_id), path(mash_sketch_output)

    output:
      tuple val(sample_id), path("${sample_id}_mash_sketch_summary.csv")

    script:
      """
      summarize_mash_sketch.py ${mash_sketch_output} --sample-id ${sample_id} > ${sample_id}_mash_sketch_summary.csv
      """
}
