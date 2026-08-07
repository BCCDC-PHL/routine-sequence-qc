process seqtk_fqchk {

    tag { meta.id }

    cpus 1
    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/seqtk_fqchk", pattern: "${meta.id}_seqtk_fqchk*_position*.csv", mode: 'copy'

    input:
    tuple val(meta), path(reads_1), path(reads_2)

    output:
    tuple val(meta), path("${meta.id}_seqtk_fqchk_all_positions.csv"), path("${meta.id}_seqtk_fqchk_by_position.csv")

    script:
    """
    echo 'filename,position,num_bases,percent_a,percent_c,percent_g,percent_t,percent_n,average_q,error_q,percent_bases_below_q${params.seqtk_fqchk_threshold},percent_bases_above_q${params.seqtk_fqchk_threshold}' > header.csv
    seqtk fqchk -q ${params.seqtk_fqchk_threshold} ${reads_1} | tr \$'\\t' ',' | tail -n+3 | awk -F ',' 'BEGIN {OFS=FS}; {print "${reads_1}", \$0}' > ${meta.id}_R1_seqtk_fqchk_data.csv
    grep 'ALL' ${meta.id}_R1_seqtk_fqchk_data.csv > ${meta.id}_R1_seqtk_fqchk_data_all_positions.csv
    grep -v 'ALL' ${meta.id}_R1_seqtk_fqchk_data.csv > ${meta.id}_R1_seqtk_fqchk_data_by_position.csv || touch ${meta.id}_R1_seqtk_fqchk_data_by_position.csv
    seqtk fqchk -q ${params.seqtk_fqchk_threshold} ${reads_2} | tr \$'\\t' ',' | tail -n+3 | awk -F ',' 'BEGIN {OFS=FS}; {print "${reads_2}", \$0}' > ${meta.id}_R2_seqtk_fqchk_data.csv
    grep 'ALL' ${meta.id}_R2_seqtk_fqchk_data.csv > ${meta.id}_R2_seqtk_fqchk_data_all_positions.csv
    grep -v 'ALL' ${meta.id}_R2_seqtk_fqchk_data.csv > ${meta.id}_R2_seqtk_fqchk_data_by_position.csv || touch ${meta.id}_R2_seqtk_fqchk_data_by_position.csv
    cat header.csv ${meta.id}_R1_seqtk_fqchk_data_all_positions.csv ${meta.id}_R2_seqtk_fqchk_data_all_positions.csv > ${meta.id}_seqtk_fqchk_all_positions.csv
    cat header.csv ${meta.id}_R1_seqtk_fqchk_data_by_position.csv ${meta.id}_R2_seqtk_fqchk_data_by_position.csv > ${meta.id}_seqtk_fqchk_by_position.csv
    """
}

