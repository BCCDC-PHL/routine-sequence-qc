process copy_to_nginx {

    tag { run_id }

    input: 
        tuple val(run_id), path(nginx_server_location), path(multiqc_file), path(basic_qc_stats), path(abundance_top_5)

    script:
    """
        cp ${abundance_top_5} ${nginx_server_location}/data/species-abundance/${run_id}_top_5_abundances_species.csv
        cp ${basic_qc_stats} ${nginx_server_location}/data/basic-qc/${run_id}_basic_qc_stats.csv
        cp ${multiqc_file} ${nginx_server_location}/data/multi-qc/${run_id}.html
    """

}