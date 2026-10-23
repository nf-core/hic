/*
 * cooltools - call_compartments
 */

process COOLTOOLS_EIGSCIS {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::cooltools=0.7.1 bioconda::ucsc-bedgraphtobigwig=377"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c81d8d6b6acf4714ffaae1a274527a41958443f6:cc7ea58b8cefc76bed985dcfe261cb276ed9e0cf-0' :
        'biocontainers/mulled-v2-c81d8d6b6acf4714ffaae1a274527a41958443f6:cc7ea58b8cefc76bed985dcfe261cb276ed9e0cf-0' }"

    input:
    tuple val(meta), path(cool), val(resolution)
    path(fasta)
    path(chrsize)

    output:
    path("*compartments*"), emit: results
    path("versions.yml"), emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cooltools genome binnify --all-names ${chrsize} ${resolution} > genome_bins.txt
    cooltools genome gc genome_bins.txt ${fasta} > genome_gc.txt

    # cooler >= 0.9 annotate() raises IndexError on an empty pixel selection when the
    # bins index does not start at 0; cooltools queries per-region bins with a global
    # index. Restore the pre-0.9 annotate semantics (graceful on empty selections) first.
    python <<'NFCORE_COOLTOOLS'
    import cooler
    import pandas as pd
    import sys

    def _annotate_compat(tiles, bins, *a, **k):
        anns = []
        for col, suffix in (("bin1_id", "1"), ("bin2_id", "2")):
            if col in tiles.columns:
                ids = tiles[col].to_numpy()
                base = bins.iloc[ids - bins.index[0]] if len(ids) else bins.iloc[0:0]
                anns.append(base.rename(columns=lambda x: x + suffix).reset_index(drop=True))
        if not anns:
            return tiles
        out = pd.concat([*anns, tiles.reset_index(drop=True)], axis=1)
        out.index = tiles.index
        return out

    cooler.annotate = _annotate_compat
    sys.argv = ["cooltools", "eigs-cis"] + "${args}".split() + ["--phasing-track", "genome_gc.txt", "-o", "${prefix}_compartments", "${cool}"]
    from cooltools.cli import cli
    sys.exit(cli())
    NFCORE_COOLTOOLS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooltools: \$(cooltools --version | grep 'cooltools, version ' | sed 's/cooltools, version //')
    END_VERSIONS
    """
}
