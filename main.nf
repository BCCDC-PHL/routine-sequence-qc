#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BCCDC-PHL/routine-sequence-qc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/BCCDC-PHL/routine-sequence-qc
----------------------------------------------------------------------------------------
*/

include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_routine_sequence_qc_pipeline'
include { ROUTINE_SEQUENCE_QC     } from './workflows/routine_sequence_qc'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_routine_sequence_qc_pipeline'


workflow BCCDCPHL_ROUTINE_SEQUENCE_QC {

  

    main:

    ROUTINE_SEQUENCE_QC (
        params.run_dir,
        params.outdir,
    )

    emit:
    multiqc_report = ROUTINE_SEQUENCE_QC.out.multiqc_report
}

workflow {

    main:

    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.run_dir,
        params.help,
        params.help_full,
        params.show_hidden
    )

    BCCDCPHL_ROUTINE_SEQUENCE_QC ()

    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
	BCCDCPHL_ROUTINE_SEQUENCE_QC.out.multiqc_report
    )
}
