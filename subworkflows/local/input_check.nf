//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    if (params.split_fastq){

      SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
	.map { create_fastq_channels(it) }
	.splitFastq( by: params.fastq_chunks_size, pe:true, file: true, compress:true)
	.map { it -> [it[0], [it[1], it[2]]]}
	.groupTuple(by: [0])
        .flatMap { it -> setMetaChunk(it) }
        .collate(2)
	//.map { it ->
	//  def meta = it[0].clone()
	//  meta.chunk = it[1].baseName - ~/.fastq(.gz)?/
	//  return [meta, [it[1], it[2]]]
	//}
        .set { reads }

    }else{
      SAMPLESHEET_CHECK ( samplesheet )
      	.csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channels(it) }
	.map { it -> [it[0], [it[1], it[2]]]}
	.groupTuple(by: [0])
        .flatMap { it -> setMetaChunk(it) }
        .collate(2)
        .set { reads }
   }

    emit:
    reads // channel: [ val(meta), [ reads ] ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channels(LinkedHashMap row) {
  def meta = [:]
  meta.id = row.sample
  meta.single_end = false

  def array = []
  if (!file(row.fastq_1).exists()) {
    exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
  }
  if (!file(row.fastq_2).exists()) {
    exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
  }
  array = [ meta, file(row.fastq_1), file(row.fastq_2) ]
  return array
}

// Set the meta.chunk value in case of technical replicates
def setMetaChunk(row){
  def map = []
  row[1].eachWithIndex() { file,i ->
    meta = row[0].clone()
    meta.chunk = i
    map += [meta, file]
  }
  return map
}