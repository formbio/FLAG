process trinity {
  label 'trinity'
  publishDir "$params.output/transcripts", mode: 'copy'
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  input:
    tuple val(nameType),path(bam)
    val filetype
  output:
    path("*rinity.fasta")
  script:
  """
  bash ${params.repoDir}/scripts/geneannotation/trinity.sh -f ${bam} -s ${filetype} -n ${nameType} -t 15 -m 30 -i 10000
  """
  stub:
  """
  touch trinity_out_dir.Trinity.fasta
  touch trinity_out_dir.Trinity.fasta.gene_trans_map
  """
}

process windowmasker {
  label 'ncbiclibraries'
  publishDir "$params.output/mask", mode: 'copy'
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  input:
    path genome
  output:
    path("maskedGenome.fa")
  script:
  """
  bash ${params.repoDir}/scripts/geneannotation/windowmasker.sh -g ${genome}
  """
  stub:
  """
  touch trinity_out_dir.Trinity.fasta
  """
}

process exonerate_p2g {
  time '24d'
  label 'exonerate_p2g'
  publishDir "$params.output/exonerate_p2g", mode: 'copy'
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  input:
    tuple path(protein),path(genome)
  output:
    path("*xonerate*.out")
  script:
  """
  bash ${params.repoDir}/scripts/geneannotation/exonerate_p2g.sh -g ${genome} -p ${protein} -n 6270 -t 15
  """
  stub:
  """
  touch trinity_out_dir.Trinity.fasta
  """
}

process liftoff {
  label 'liftoff'
  publishDir "$params.output/liftoff", mode: 'copy'
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  input:
    path query
    path target
    path reference
    val gaps
  output:
    path("liftoff.gff3")
  script:
  """
  bash ${params.repoDir}/scripts/geneannotation/liftoff.sh -q ${query} -t ${target} -r ${reference} -g ${gaps}
  """
  stub:
  """
  touch trinity_out_dir.Trinity.fasta
  """
}

process prosplign {
  time '24d'
  label 'prosplign'
  publishDir "$params.output/prosplign", mode: 'copy'
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  input:
    tuple path(genome),path(protein)
  output:
    path("*.gff3")
  script:
  """
  bash ${params.repoDir}/scripts/geneannotation/tblastn_prosplign.sh -g ${genome} -p ${protein} -n 1 -t 15
  """
  stub:
  """
  touch trinity_out_dir.Trinity.fasta
  """
}

process splign {
  label 'ncbiclibraries'
  publishDir "$params.output/splign", mode: 'copy'
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  stageInMode = 'copy'
  stageOutMode = 'copy'
  input:
    path genome
    path rna
  output:
    path("*.gff3")
  script:
  """
  bash ${params.repoDir}/scripts/geneannotation/splign.sh -g ${genome} -r ${rna}
  """
  stub:
  """
  touch trinity_out_dir.Trinity.fasta
  """
}

process transdecoder {
  label 'transdecoder'
  publishDir "$params.output/transdecoder", mode: 'copy'
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  input:
    path transcript
    path annotation
    path genome
    val threads
  output:
    path("*ransdecoder.gff3")
  script:
  """
  bash ${params.repoDir}/scripts/geneannotation/transdecoder.sh -i ${transcript} -a ${annotation} -g ${genome} -t ${threads}
  """
  stub:
  """
  touch trinity_out_dir.Trinity.fasta
  """
}

process pasa {
  label 'pasa'
  publishDir "$params.output/pasa", mode: 'copy'
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  stageInMode = 'copy'
  stageOutMode = 'copy'
  input:
    path genome
    path transcript
  output:
    path("*_assembly.gff3"), emit: assembly
    path("sample_mydb_pasa.sqlite.valid_blat_alignments.gtf"), emit: alignments
    path("*asa_pblat.pslx"), emit: pslx
  script:
  """
  bash ${params.repoDir}/scripts/geneannotation/pasa.sh -g ${genome} -r ${transcript} -t 15
  """
  stub:
  """
  touch trinity_out_dir.Trinity.fasta
  """
}

process genome2protein {
  label 'genome2protein'
  publishDir "$params.output/refseq_protein", mode: 'copy'
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  input:
    tuple path(genome),path(database)
  output:
    path("referenceProtein*.fa")
  script:
  """
  bash ${params.repoDir}/scripts/geneannotation/genome2protein.sh -i ${genome} -d ${database} -a diamond -t 15
  """
  stub:
  """
  touch referenceProtein.fa
  """
}

process genome2transcriptome {
  label 'ncbiclibraries'
  publishDir "$params.output/refseq_transcriptome", mode: 'copy'
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  input:
    path genome
    path database
  output:
    path("referenceRNA.fa")
  script:
  """
  bash ${params.repoDir}/scripts/geneannotation/genome2transcriptome.sh -i ${genome} -d ${database} -t 15
  """
  stub:
  """
  touch referenceRNA.fa
  """
}

process trnascan {
  label 'ncbiclibraries'
  publishDir "$params.output/tRNA", mode: 'copy'
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  input:
    path genome
  output:
    path("tRNAScan.gff3")
  script:
  """
  bash ${params.repoDir}/scripts/geneannotation/tRNAScan.sh -g ${genome} -t 15
  """
  stub:
  """
  touch tRNAScan.gff3
  """
}
process maskingNoLibrary {
  label "${ params.algo == 'WindowMasker' ? 'ncbiclibraries' : 'tetools' }"
  publishDir "$params.output/mask", mode: 'copy'
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  input:
    tuple path(genome),val(algo)
  output:
    path("*masked.fa"), emit: genome
    path("*.out.gff"), emit: gff optional true
  script:
  """
  bash ${params.repoDir}/scripts/geneannotation/mask.sh -i ${genome} -a ${algo} -t 15
  """
  stub:
  """
  touch trinity_out_dir.Trinity.fasta
  """
}
process entap {
  label 'entap'
  publishDir "$params.output/entap", mode: 'copy'
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  input:
    tuple path(input),path(annotation),path(databases)
  output:
    path("final_annotations_lvl0.tsv")
  script:
  """
  bash ${params.repoDir}/scripts/geneannotation/entap.sh -i ${input} -a ${annotation} -d ${databases}
  """
  stub:
  """
  touch final_annotations_lvl0.tsv
  """
}
process CombineAndFilter {
  label 'evm'
  publishDir "$params.output/StructuralAnnotation", mode: 'copy'
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  stageInMode = 'copy'
  stageOutMode = 'copy'
  input:
    path genome
    path files
    val lineage
    val size
  output:
    path("weights*"), emit: weights
    path("FinalStructuralAnnotationStrictFilter.gff3"), emit: gff3Strict
    path("FinalStructuralAnnotationLenientFilter.gtf"), emit: gtfLenient
    path("FinalStructuralAnnotationLenientFilter.gff3"), emit: gff3Lenient
    path("FinalStructuralAnnotation*Filter.gff3"), emit: gff3s
    path("FinalStructuralAnnotation*.gff3.stats"), emit: stats
    path("short_summary.specific.*.txt"), emit: buscos
  script:
  """
  bash ${params.repoDir}/scripts/geneannotation/filterEukGenes.sh -i ${genome} -l ${lineage} -s ${size}
  """
  stub:
  """
  touch testfinalstructannot.gff3
  """
}
process Helixer {
  label 'helixer'
  publishDir "$params.output/Helixer", mode: 'copy'
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  input:
    path genome
    val name
    val model
    val size
  output:
    path("Helixer.gff3"), emit: gff3
  script:
  """
  bash ${params.repoDir}/scripts/geneannotation/helixer.sh -i ${genome} -n ${name} -m ${model} -s ${size}
  """
  stub:
  """
  touch testfinalstructannot.gff3
  """
}
process augustus {
  time '24d'
  label 'augustus'
  publishDir "$params.output/augustus", mode: 'copy'
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  stageInMode = 'copy'
  stageOutMode = 'copy'
  input:
    path genome
    path transcript
    val species
    path exoneratehints
    val pretrained
    path extrahints
  output:
    path("*.gtf")
  script:
  """
  bash ${params.repoDir}/scripts/geneannotation/augustus_pasa.sh -g ${genome} -r ${transcript} -s ${species} -n 15 -e ${exoneratehints} -p ${pretrained}
  """
  stub:
  """
  touch augustus.gtf
  """
}

process otherToolTrainAugustus {
  time '24d'
  label 'augustus'
  publishDir "$params.output/augustus", mode: 'copy'
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  stageInMode = 'copy'
  stageOutMode = 'copy'
  input:
    path genome
    path transcript
    val species
    path extrahints
  output:
    path("*_augustus.gtf")
  script:
  """
  bash ${params.repoDir}/scripts/geneannotation/augustus_pasa.sh -g ${genome} -r ${transcript} -s ${species} -n 15 -p False
  """
  stub:
  """
  touch augustus.gtf
  """
}
process genomethreader {
  time '24d'
  label 'exonerate_p2g'
  publishDir "$params.output/genomethreader", mode: 'copy'
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  input:
    tuple path(protein),path(genome)
  output:
    path("*enomethreader*.gff3")
  script:
  """
  bash ${params.repoDir}/scripts/geneannotation/genomethreader.sh -g ${genome} -p ${protein} -t 15
  """
  stub:
  """
  touch trinity_out_dir.Trinity.fasta
  """
}
process miniprot {
  time '24d'
  label 'exonerate_p2g'
  publishDir "$params.output/miniprot", mode: 'copy'
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  input:
    tuple path(protein),path(genome)
  output:
    path("*iniprot*.gtf")
  script:
  """
  bash ${params.repoDir}/scripts/geneannotation/miniprot.sh -g ${genome} -p ${protein} -t 15
  """
  stub:
  """
  touch trinity_out_dir.Trinity.fasta
  """
}
process combinestructwfunct {
  label 'evm'
  publishDir "$params.output/finalAnnots", mode: 'copy'
  publishDir "${params.output}", pattern: '*html', mode: 'copy'
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  input:
    tuple path(entap),path(annotation),val(species),val(lineage),path(genome)
  output:
    path("final*.gtf")
    path("final*.gff3")
    path("final*.AGAT.stats")
    path("cdna_*.fa")
    path("proteins_*.fa")
    path("short_summary.*.buscoout.txt")
    path("results.html")
    path("*.gb") optional true
  script:
  """
  bash ${params.repoDir}/scripts/geneannotation/parseEntap.sh -i ${entap} -a ${annotation} -s ${species} -l ${lineage} -g ${genome}
  """
  stub:
  """
  touch final_annotations_lvl0.tsv
  """
}
