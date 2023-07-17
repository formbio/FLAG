process splitgenomebysize {
  label 'formatting'
  publishDir "$params.output/splitgenomebysize", mode: 'copy'
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
   input:
    path fa
    val num
    val size
  output:
    file("*.fasta")
  script:
  """
  bash ${params.repoDir}/scripts/formatting/genomesplitter.sh -i ${fa} -n ${num} -s ${size}
  """
  stub:
  """
  touch test.fasta
  """
}

process combiner {
  label 'formatting'
  publishDir "$params.output/combined", mode: 'copy'
   input:
    path fa
    val type
  output:
    path("*xonerate.out"), emit: exoner optional true
    path("*rosplign.gff3"), emit: prosplign optional true
    path("*enomethreader.gff3"), emit: gth optional true
    path("referenceProtein.fa"), emit: referenceprotein optional true
    path("basecall.*"), emit: basecallcombined optional true
    path("maskedGenome.fa"), emit: masked optional true
  script:
  """
  bash ${params.repoDir}/scripts/formatting/combiner.sh -p ${type}
  """
  stub:
  """
  touch test.fasta
  """
}

process removedupsfa {
  label 'formatting'
  publishDir "$params.output/fafiles", mode: 'copy'
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
  input:
    path fa
  output:
    file("*.fa")
  script:
  """
  bash ${params.repoDir}/scripts/formatting/removedupsinfa.sh -i ${fa}
  """
  stub:
  """
  touch formatted_referenceRNA.fa
  """
}
