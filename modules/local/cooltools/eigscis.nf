/*
 * cooltools - call_compartments
 */

process COOLTOOLS_EIGSCIS {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::cooltools=0.7.1 bioconda::ucsc-bedgraphtobigwig=377"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/27/274f5820b192395eea6566d62cbbe1d41653d6e4458d358cf677fb0dbacb6150/data' :
        'community.wave.seqera.io/library/cooltools_ucsc-bedgraphtobigwig:cbb762ca30048223' }"

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

    if not hasattr(cooler.Cooler, "_load_attrs"):
        cooler.Cooler._load_attrs = lambda self, path: dict(self.open("r")[path].attrs)
    if tuple(int(x) for x in cooler.__version__.split(".")[:2]) >= (0, 9):
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
