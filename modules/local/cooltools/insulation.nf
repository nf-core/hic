/*
 * Cooltools - diamond-insulation
 */

process COOLTOOLS_INSULATION {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::cooltools=0.7.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooltools:0.5.1--py37h37892f8_0' :
        'biocontainers/cooltools:0.5.1--py37h37892f8_0' }"

    input:
    tuple val(meta), path(cool)

    output:
    path("*tsv"), emit:tsv
    path("versions.yml"), emit:versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # cooler >= 0.9 annotate() raises IndexError on an empty pixel selection when the
    # bins index does not start at 0. cooltools insulation queries per-region bins with
    # a global index and hits this on regions without pixels inside the diagonal window.
    # Restore the pre-0.9 annotate semantics (graceful on empty selections) first.
    python <<'NFCORE_COOLTOOLS' > ${prefix}_insulation.tsv
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
    sys.argv = ["cooltools", "insulation", "${cool}"] + "${args}".split()
    from cooltools.cli import cli
    sys.exit(cli())
    NFCORE_COOLTOOLS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooltools: \$(cooltools --version | grep 'cooltools, version ' | sed 's/cooltools, version //')
    END_VERSIONS
    """
}
