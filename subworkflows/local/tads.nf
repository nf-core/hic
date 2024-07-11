include { COOLTOOLS_INSULATION } from '../../modules/local/cooltools/insulation'
include { HIC_FIND_TADS } from '../../modules/local/hicexplorer/hicFindTADs'

workflow TADS {

    take:
    cool

    main:
    ch_versions = Channel.empty()
    ch_tads = Channel.empty()

    if (params.tads_caller =~ 'insulation'){
        COOLTOOLS_INSULATION(cool)
        ch_versions = ch_versions.mix(COOLTOOLS_INSULATION.out.versions)
        ch_tads = ch_tads.mix(COOLTOOLS_INSULATION.out.tsv)
    }

    if (params.tads_caller =~ 'hicexplorer'){
        HIC_FIND_TADS(cool)
        ch_versions = ch_versions.mix(HIC_FIND_TADS.out.versions)
        ch_tads = ch_tads.mix(HIC_FIND_TADS.out.results)
    }

    emit:
    tads = ch_tads
    versions = ch_versions
}
