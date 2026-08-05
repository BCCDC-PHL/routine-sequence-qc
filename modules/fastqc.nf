process fastqc {

    tag { meta.id }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${meta.id}_*_fastqc", mode: 'copy'

    cpus 2

    input:
      tuple val(meta), path(reads_1), path(reads_2)

    output:
      tuple val(meta), path("${meta.id}_R1_fastqc"), path("${meta.id}_R2_fastqc")

    script:
    """
    mkdir -p ./tmp
    
    fastqc \
	-t ${task.cpus} \
	--dir ./tmp \
	${reads_1} \
	${reads_2}

    for d in *.zip; do unzip \$d; done
    if [ ! -d "${meta.id}_R1_fastqc" ]; then
      mv ${meta.id}*_R1*_fastqc ${meta.id}_R1_fastqc
    fi
    if [ ! -d "${meta.id}_R2_fastqc" ]; then
      mv ${meta.id}*_R2*_fastqc ${meta.id}_R2_fastqc
    fi
    """
}
