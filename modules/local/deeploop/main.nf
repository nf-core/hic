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
    python3 /DeepLoop-master/prediction/predict_chromosome.py --full_matrix_dir \$HiCorr_path/ \\
                                                --input_name anchor_2_anchor.loop.\$chr.p_val \\
                                                --h5_file /DeepLoop-master/DeepLoop_models/CPGZ_trained/LoopDenoise.h5 \\
                                                --json_file /DeepLoop-master/DeepLoop_models/CPGZ_trained/LoopDenoise.json \\
                                                --out_dir output/ \\
                                                --anchor_dir /DeepLoop-master/DeepLoop_models/ref/hg19_HindIII_anchor_bed/ \\
                                                --chromosome \$chr \\
                                                --small_matrix_size 128 \\
                                                --step_size 128 \\
                                                --dummy 5 \\
                                                --val_cols obs exp pval

    python3 /DeepLoop-master/utils/convert_to_cooler.py --anchor_dir /DeepLoop-master/training_data/anchor_bed/ \
                                        --loop_dir output/H9_denoised/ \
                                        --out_file coolers/H9_denoise_chr11.cool \
                                        --col_names a1 a2 denoise \
                                        --cooler_col denoise \
                                        --single_chrom chr11;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/Python //')
    END_VERSIONS
    """
}
