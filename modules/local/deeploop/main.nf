// TODO Split this up into seperate modules
process DEEPLOOP {
    tag "$meta.id"
    label 'process_high_memory'

    // TODO There's no conda package 🤠
    container 'nf-core/deeploop:1.0.1'

    input:
    tuple val(meta), path(tar_file)

    output:
    path("versions.yml"), emit: versions

    script:
    """
    tar -xvf $tar_file
    HiCorr_path=HiCorr_output
    chr=chr11
    python3 /DeepLoop-master/prediction/predict_chromosome.py --full_matrix_dir \$HiCorr_path/ \
                                                --input_name anchor_2_anchor.loop.\$chr.p_val \
                                                --h5_file /DeepLoop-master/DeepLoop_models/CPGZ_trained/LoopDenoise.h5 \
                                                --json_file /DeepLoop-master/DeepLoop_models/CPGZ_trained/LoopDenoise.json \
                                                --out_dir output/ \
                                                --anchor_dir /DeepLoop-master/DeepLoop_models/ref/hg19_HindIII_anchor_bed/ \
                                                --chromosome \$chr \
                                                --small_matrix_size 128 \
                                                --step_size 128 \
                                                --dummy 5 \
                                                --val_cols obs exp pval

    ls output
    head output/\$chr.denoised.anchor.to.anchor

    # Visulaize contact heatmaps from raw, HiCorr, and DeepLoop given a genomic location chr start end
    chr=chr11
    start=130000000
    end=130800000
    outplot="./test"
    /DeepLoop-master/lib/generate.matrix.from_HiCorr.pl DeepLoop/DeepLoop_models/ref/hg19_HindIII_anchor_bed/\$chr.bed \$HiCorr_path/anchor_2_anchor.loop.\$chr \$chr \$start \$end ./\${chr}_\${start}_\${end}
    /DeepLoop-master/lib/generate.matrix.from_DeepLoop.pl DeepLoop/DeepLoop_models/ref/hg19_HindIII_anchor_bed/\$chr.bed output/\$chr.denoised.anchor.to.anchor \$chr \$start \$end ./\${chr}_\${start}_\${end}
    /DeepLoop-master/lib/plot.multiple.r \$outplot 1 3 \${chr}_\${start}_\${end}.raw.matrix \${chr}_\${start}_\${end}.ratio.matrix \${chr}_\${start}_\${end}.denoise.matrix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/Python //')
    END_VERSIONS
    """
}
