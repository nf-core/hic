process MERGE_BOWTIE2{
  tag "$prefix"
  label 'process_medium'

  input:
  tuple val(meta), path(bam1), path(bam2) 

  output:
  tuple val(meta), path("${prefix}_bwt2merged.bam"), emit: bam
  tuple val(meta), path("${prefix}.mapstat"), emit: stats

  script:
  prefix = meta.id + "_" + meta.mates
  tag = meta.mates
  """
  samtools merge -@ ${task.cpus} \\
                 -f ${prefix}_bwt2merged.bam \\
                  ${bam1} ${bam2}

  samtools sort -@ ${task.cpus} -m 800M \\
                -n  \\
                -o ${prefix}_bwt2merged.sorted.bam \\
                ${prefix}_bwt2merged.bam

  mv ${prefix}_bwt2merged.sorted.bam ${prefix}_bwt2merged.bam

  echo "## ${prefix}" > ${prefix}.mapstat
  echo -n "total_${tag}\t" >> ${prefix}.mapstat
  samtools view -c ${prefix}_bwt2merged.bam >> ${prefix}.mapstat
  echo -n "mapped_${tag}\t" >> ${prefix}.mapstat
  samtools view -c -F 4 ${prefix}_bwt2merged.bam >> ${prefix}.mapstat
  echo -n "global_${tag}\t" >> ${prefix}.mapstat
  samtools view -c -F 4 ${bam1} >> ${prefix}.mapstat
  echo -n "local_${tag}\t"  >> ${prefix}.mapstat
  samtools view -c -F 4 ${bam2} >> ${prefix}.mapstat
  """
}
