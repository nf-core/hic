process FITHICHIP {
    tag "$meta.id"
    label 'process_medium'
    label 'process_long'

    // FitHiChIP has no Bioconda recipe. The image is pinned by digest: FitHiChIP 11.0 (December 2022)
    // on the HiC-Pro 3.0.0 environment (python3, R, bedtools, samtools, bc).
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://riasc/fithichip@sha256:af69295105cecfd2001e4f8f090a4ff60c0f62c234a92f771863dd32ad22d802' :
        'riasc/fithichip@sha256:af69295105cecfd2001e4f8f090a4ff60c0f62c234a92f771863dd32ad22d802' }"

    input:
    tuple val(meta), path(bins), path(matrix), path(peaks), path(chromsize)
    val lower_dist
    val upper_dist
    val qvalue
    val int_type
    val bias_type
    val p2p_background

    output:
    tuple val(meta), path("${prefix}/FitHiChIP_*/P2PBckgr_*/*/FitHiC_BiasCorr/*.interactions_FitHiC_Q*.bed")      , emit: interactions
    tuple val(meta), path("${prefix}/FitHiChIP_*/P2PBckgr_*/*/FitHiC_BiasCorr/*_IGV.bedpe")                       , emit: bedpe   , optional: true
    tuple val(meta), path("${prefix}/FitHiChIP_*/P2PBckgr_*/*/FitHiC_BiasCorr/*_WashU.bed.gz*")                   , emit: washu   , optional: true
    tuple val(meta), path("${prefix}/FitHiChIP_*/P2PBckgr_*/*/FitHiC_BiasCorr/Merge_Nearby_Interactions")         , emit: merged  , optional: true
    tuple val(meta), path("${prefix}/FitHiChIP_*/P2PBckgr_*/*/FitHiC_BiasCorr/Statistics_Dist.txt")               , emit: stats   , optional: true
    tuple val(meta), path("${prefix}/FitHiChIP_*/P2PBckgr_*/*/FitHiC_BiasCorr/Plots")                             , emit: plots   , optional: true
    tuple val(meta), path("${prefix}/Summary_results_FitHiChIP.html")                                              , emit: summary , optional: true
    tuple val(meta), path("${prefix}/Parameters.txt")                                                              , emit: params  , optional: true
    tuple val(meta), path("${prefix}.fithichip.config")                                                            , emit: config
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FITHICHIP module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    prefix = task.ext.prefix ?: "${meta.id}"
    def bin_size = meta.resolution
    def p2p = p2p_background ? 1 : 0
    // FitHiChIP reads a key=value configuration file (FitHiChIP_HiCPro.sh -C) with absolute paths.
    // Interval/Matrix: HiC-Pro style bins and raw contact counts at the loop resolution; PeakFile: the 1D anchors.
    // The script changes into its own code directory and R writes Rplots.pdf there, which fails in a read-only
    // container (Singularity): a writable copy of the code is used instead (test data and docs left out).
    """
    cat > ${prefix}.fithichip.config <<-FITHICHIP_CONFIG
    Interval=\$PWD/${bins}
    Matrix=\$PWD/${matrix}
    PeakFile=\$PWD/${peaks}
    ChrSizeFile=\$PWD/${chromsize}
    OutDir=\$PWD/${prefix}
    PREFIX=${prefix}
    IntType=${int_type}
    BINSIZE=${bin_size}
    LowDistThr=${lower_dist}
    UppDistThr=${upper_dist}
    UseP2PBackgrnd=${p2p}
    BiasType=${bias_type}
    MergeInt=1
    QVALUE=${qvalue}
    CircularGenome=0
    OverWrite=1
    FITHICHIP_CONFIG

    mkdir -p fithichip_code
    cp -r /root/FitHiChIP/src /root/FitHiChIP/Analysis /root/FitHiChIP/Imp_Scripts /root/FitHiChIP/UtilScript fithichip_code/
    cp /root/FitHiChIP/*.sh /root/FitHiChIP/Release_Notes.txt /root/FitHiChIP/interact.as /root/FitHiChIP/bedToBigBed fithichip_code/
    export PATH=\$PWD/fithichip_code:\$PATH

    bash fithichip_code/FitHiChIP_HiCPro.sh -C \$PWD/${prefix}.fithichip.config

    rm -rf fithichip_code

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fithichip: \$(grep -m1 -oE 'Version [0-9.]+' /root/FitHiChIP/Release_Notes.txt | sed 's/Version //')
    END_VERSIONS
    """
}
