process FITHICHIP {
    tag "$meta.id"
    label 'process_medium'
    label 'process_long'

    // FitHiChIP has no Bioconda recipe. The image is pinned by digest: FitHiChIP 11.0 (December 2022) with its
    // dependencies from conda-forge/Bioconda, built from https://github.com/birkiy/fithichip-container
    // (ghcr.io/birkiy/fithichip:11.0; the Dockerfile, the conda lock file and the build workflow are in that repository).
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/birkiy/fithichip@sha256:ae01888ea5c881470b8e0a225c8b7ba48a9e1c40aba6e984c61c87d3222b16d3' :
        'ghcr.io/birkiy/fithichip@sha256:ae01888ea5c881470b8e0a225c8b7ba48a9e1c40aba6e984c61c87d3222b16d3' }"

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
    // container (Singularity): a writable copy of the code is used instead (test data left out). The code
    // directory is taken from FITHICHIP_HOME when the image sets it, else from the /root/FitHiChIP of older images.
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

    FITHICHIP_HOME=\${FITHICHIP_HOME:-/root/FitHiChIP}
    mkdir -p fithichip_code
    for entry in "\$FITHICHIP_HOME"/*; do
        [ "\$(basename "\$entry")" = "TestData" ] || cp -r "\$entry" fithichip_code/
    done
    export PATH=\$PWD/fithichip_code:\$PATH
    # writable font cache for the R plots (the home directory can be read-only in a container)
    export XDG_CACHE_HOME=\$PWD/.cache

    bash fithichip_code/FitHiChIP_HiCPro.sh -C \$PWD/${prefix}.fithichip.config

    rm -rf fithichip_code

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fithichip: \$(grep -m1 -oE 'Version [0-9.]+' "\$FITHICHIP_HOME"/Release_Notes.txt | sed 's/Version //')
    END_VERSIONS
    """
}
