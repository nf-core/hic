include { COOLTOOLS_EIGSCIS } from '../../modules/local/cooltools/eigscis'
include { CALDER2 } from '../../modules/nf-core/calder2/main'

workflow COMPARTMENTS {

    take:
    cool
    fasta
    chrsize

    main:
    ch_versions = Channel.empty()

    if (params.compartments_caller =~ 'cooltools'){
        COOLTOOLS_EIGSCIS(
            cool,
            fasta.map{it -> it[1]}.collect(),
            chrsize.map{it -> it[1]}.collect()
        )
        ch_versions = ch_versions.mix(COOLTOOLS_EIGSCIS.out.versions)
        ch_comp = COOLTOOLS_EIGSCIS.out.results
    }

    if (params.compartments_caller =~ 'calder2'){
        CALDER2(
            cool.map{meta, cool, res -> [meta, cool] },
            Channel.value([])
        )
        ch_versions = ch_versions.mix(CALDER2.out.versions)
        ch_comp = CALDER2.out.output_folder
    }

    emit:
    versions = ch_versions
    compartments = ch_comp
}
