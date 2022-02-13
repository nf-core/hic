params.options = [:]

include { COOLER_RAW } from '../../modules/local/cooler_raw' addParams( options: params.options )
include { COOLER_BALANCE } from '../../modules/local/cooler_balance' addParams( options: params.options )
include { COOLER_ZOOMIFY } from '../../modules/local/cooler_zoomify' addParams( options: params.options )

workflow COOLER {

    take:


    main:
    

    emit:
    
}