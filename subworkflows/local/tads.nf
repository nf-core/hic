params.options = [:]

include { TADS_HICEXPLORER } from '../../modules/local/tads_hicexplorer' addParams( options: params.options )
include { TADS_INSULATION } from '../../modules/local/tads_insulation' addParams( options: params.options )

workflow TADS {

    take:


    main:
    

    emit:
    
}