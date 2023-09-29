#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

repoDir=workflow.projectDir

if (params.repoDir) {
  repoDir=params.repoDir
}
//input data
genome = Channel.fromPath(params.genome)
masker = Channel.from(params.masker)
entapDB = Channel.fromPath(params.entapDB)
// Either RepeatMasker or WindowMasker. Repeat for highly characterised species with Dfam libraries and Window for others. Window is ncbi
 
log.info """\
 Test - N F   P I P E L I N E
 ===================================
 outdir               : ${params.output}
 masker               : ${params.masker}
 genome               : ${params.genome}
 proteins             : ${params.proteins}
 rna                  : ${params.rna}
 reference_genome     : ${params.fafile}
 reference_annotation : ${params.gtffile}
 protein database     : ${params.blastdb}
 rna database         : ${params.rnaDB}
 transcriptIn         : ${params.transcriptIn}
 Busco Lineage        : ${params.lineage}
 entapDB              : ${params.entapDB}
 Augustus Pretrained Species   : ${params.pretrainedAugustusSpecies}
 Helixer Model                 : ${params.helixerModel}
 Helixer Models Available      : verterbrate, invertebrate, land_plant, fungi
 Genome Size                   : ${params.size}
 Species Scientific Name       : ${params.speciesScientificName}
 
 all annotation algos options  : Helixer, Liftoff, denovo_augustus, related_species_augustus, augustus_pretrained, liftoff_trained_augustus, helixer_trained_augustus, transdecoder
 chosen annotation algos       : ${params.annotationalgo}
 
 all external algos options    : input_transcript, input_proteins, transcript_from_database, proteins_from_database
 chosen external algos         : ${params.externalalgo}
 
 all protein algos             : exonerate, genomethreader, prosplign, miniprot
 chosen protein algos          : ${params.proteinalgo}
"""

include { trinity; exonerate_p2g; exonerate_p2g as homology_exonerate_p2g; augustus; augustus as related_species_augustus; augustus as augustus_pretrained; liftoff; windowmasker; prosplign; prosplign as homology_prosplign; splign; splign as homology_splign; transdecoder; transdecoder as homology_transdecoder; pasa; pasa as homology_pasa; genome2protein; genome2transcriptome; CombineAndFilter; trnascan; entap; Helixer; genomethreader; genomethreader as homology_genomethreader; otherToolTrainAugustus as liftoff_trained_augustus; otherToolTrainAugustus as helixer_trained_augustus; miniprot; miniprot as homology_miniprot; combinestructwfunct; maskingNoLibrary } from './modules/nf-modules/annotation/main.nf'
include { removedupsfa as protein_removedupsfa; removedupsfa as rna_removedupsfa; removedupsfa as rprotein_removedupsfa; removedupsfa as rrna_removedupsfa; splitgenomebysize; splitgenomebysize as splitgenomebysizeMasking; combiner as rexoncombiner; combiner as exoncombiner; combiner as rprocombiner; combiner as procombiner; combiner as g2pcombiner; combiner as gthcombiner; combiner as rgthcombiner; combiner as maskingcombiner } from './modules/nf-modules/formatting/main.nf'

annotationalgos = params.annotationalgo.tokenize(',')
externalalgos = params.externalalgo.tokenize(',')
proteinalgos = params.proteinalgo.tokenize(',')

file("${params.output}/emptyProteinPlaceHolder.txt").text = "\n"
file("${params.output}/emptyTranscriptPlaceHolder.txt").text = "\n"
file("${params.output}/emptyPlaceHolder.txt").text = "\n"

workflow {
  //Mask the genomic sequence
 if ( params.masker == 'WindowMasker' ) {
   windowmasker(genome)
   maskedGenome = windowmasker.out
 } else if ( params.masker == 'RepeatMasker_with_RepeatModeler' ) {
  splitgenomebysizeMasking(genome,1200,20000000)
  splitgenomebysizeMasking.out
   .flatten()
   .set { g2p_split_genome_unmasked }
  maskingNoLibrary(g2p_split_genome_unmasked.combine(Channel.from("RepeatMasker_with_RepeatModeler"))).genome
  maskingcombiner(maskingNoLibrary.out.genome.flatten().toList(),"masker")
    maskingcombiner.out.masked
    .flatten()
    .set { maskedGenome }
 } else if ( params.masker == 'RepeatMasker' ) {
  splitgenomebysizeMasking(genome,1200,20000000)
  splitgenomebysizeMasking.out
   .flatten()
   .set { g2p_split_genome_unmasked }
  maskingNoLibrary(g2p_split_genome_unmasked.combine(Channel.from("RepeatMasker"))).genome
  maskingcombiner(maskingNoLibrary.out.genome.flatten().toList(),"masker")
  maskingcombiner.out.masked
    .flatten()
    .set { maskedGenome }
 } else if (params.masker == 'skip' ) {
   maskedGenome = genome
 }  else{
   //repeatmasker(genome)
 }
 trnascan(maskedGenome)
 
 if (externalalgos =~ /proteins/) {
    if ((proteinalgos =~ /exonerate/) || (proteinalgos =~ /prosplign/) || (proteinalgos =~ /genomethreader/) || (externalalgos =~ /proteins_from_database/)) {
        splitgenomebysize(maskedGenome,1200,20000000)
        splitgenomebysize.out
         .flatten()
         .set { g2p_split_genome }
    }
 }

 if (( params.proteins ) && (externalalgos =~ /input_proteins/)) {
   proteins = Channel.fromPath(params.proteins)
   formatted_proteins = protein_removedupsfa(proteins)
   if (proteinalgos =~ /exonerate/) {
     exonerate_p2g(formatted_proteins.combine(g2p_split_genome))
     exoncombiner(exonerate_p2g.out.flatten().toList(),"exonerate")
     exoncombiner.out.exoner
      .flatten()
      .set { exonerate }
   } else {
     Channel.empty().set { exonerate }
   }
   if (proteinalgos =~ /miniprot/) {
     mprot = miniprot(formatted_proteins.combine(maskedGenome))
    } else {
     Channel.empty().set { mprot }
   }
   if (proteinalgos =~ /genomethreader/) {
    genomethreader(formatted_proteins.combine(g2p_split_genome))
    gthcombiner(genomethreader.out.flatten().toList(),"genomethreader")
      gthcombiner.out.gth
        .flatten()
        .set { gth }
   } else {
     Channel.empty().set { gth }
   }
   if (proteinalgos =~ /prosplign/) {
    prosplign(g2p_split_genome.combine(formatted_proteins))
    procombiner(prosplign.out.flatten().toList(),"prosplign")
    procombiner.out.prosplign
        .flatten()
        .set { pro }
   } else {
     Channel.empty().set { pro }
   }
   Channel.empty().mix(pro, exonerate, gth, mprot).flatten().toList().set { proteinFile }
 } else {
   Channel.empty().set { proteinFile }
 }
 
 
 if (externalalgos =~ /proteins_from_database/) {
     proteinDB = Channel.fromPath(params.blastdb)
     genome2protein(g2p_split_genome.combine(proteinDB))
     g2pcombiner(genome2protein.out.flatten().toList(),"genome2protein")
     g2pcombiner.out.referenceprotein
       .flatten()
       .set { rProtein }
     formatted_rproteins = rprotein_removedupsfa(rProtein)
   if (proteinalgos =~ /exonerate/) {
     homology_exonerate_p2g(formatted_rproteins.combine(g2p_split_genome))
     rexoncombiner(homology_exonerate_p2g.out.flatten().toList(),"rexonerate")
     rexoncombiner.out.exoner
      .flatten()
      .set { rexonerate }
   } else {
     Channel.empty().set { rexonerate }
   }
   if (proteinalgos =~ /miniprot/) {
     rmprot = homology_miniprot(formatted_rproteins.combine(maskedGenome))
    } else {
     Channel.empty().set { rmprot }
    }
    if (proteinalgos =~ /genomethreader/) {
     homology_genomethreader(formatted_rproteins.combine(g2p_split_genome))
     rgthcombiner(genomethreader.out.flatten().toList(),"rgenomethreader")
     rgthcombiner.out.gth
      .flatten()
      .set { rgth }
    } else {
     Channel.empty().set { rgth }
    }
    if (proteinalgos =~ /prosplign/) {
      homology_prosplign(g2p_split_genome.combine(formatted_rproteins))
      rprocombiner(homology_prosplign.out.flatten().toList(),"rprosplign")
      rprocombiner.out.prosplign
       .flatten()
       .set { rpro }
     } else {
       Channel.empty().set { rpro }
     }
     Channel.empty().mix(rpro, rexonerate, rgth, rmprot).flatten().toList().set { referenceProteinFile }
  } else {
    Channel.empty().set { referenceProteinFile }
  }
  
  
  if (( params.rna ) && (externalalgos =~ /input_transcript/)) {
   rna = Channel.fromPath(params.rna)
   formatted_rna = rna_removedupsfa(rna)
   //formatted_rna = rna
   if ( params.transcriptIn != false ) {
     transcript = formatted_rna
   } else {
     transcript = trinity(Channel.from('input').combine(formatted_rna),'fa')
   }
   pasa(maskedGenome,transcript)
   if (annotationalgos =~ /transdecoder/) {
      trans = transdecoder(transcript,pasa.out.alignments,maskedGenome,'16')
   } else {
      Channel.empty().set { trans }
   }
   lign = splign(genome,transcript)
   Channel.empty().mix(lign, trans, pasa.out.assembly, pasa.out.pslx).flatten().toList().set { rnaFile }
 } else {
   Channel.empty().set { rnaFile }
 }
 if (externalalgos =~ /transcript_from_database/) {
   rnaDB = Channel.fromPath(params.rnaDB)
   rRNA = genome2transcriptome(maskedGenome,rnaDB)
   formatted_rrna = rrna_removedupsfa(rRNA)
   rtranscript = formatted_rrna
 
   //refseq is comprised of transcripts
   homology_pasa(maskedGenome,rtranscript)
   if (annotationalgos =~ /transdecoder/) {
     rtrans = homology_transdecoder(rtranscript,homology_pasa.out.alignments,maskedGenome,'16')
   } else {
     Channel.empty().set { rtrans }
   }
   rlign = homology_splign(maskedGenome,rtranscript)
   Channel.empty().mix(rlign, rtrans, homology_pasa.out.assembly, homology_pasa.out.pslx).flatten().toList().set { referenceRnaFile }
 } else {
   Channel.empty().set { referenceRnaFile }
 }

 if (annotationalgos =~ /denovo_augustus/) {
      if (( params.proteins ) && (externalalgos =~ /input_proteins/)) {
        exoneratehints = exonerate
      } else if (externalalgos =~ /proteins_from_database/) {
        exoneratehints = rexonerate
      } else {
        exoneratehints = file("${params.output}/emptyProteinPlaceHolder.txt")
      }
      emptyAugPlaceholder = file("${params.output}/emptyTranscriptPlaceHolder.txt")
      if (( params.rna ) && (externalalgos =~ /input_transcript/)) {
        augustus(maskedGenome,transcript,'sample',exoneratehints,"false",emptyAugPlaceholder)
        Channel.empty().mix(augustus.out).flatten().toList().set { denovoAugFile }
      } else if (externalalgos =~ /transcript_from_database/) {
        augustus(maskedGenome,rtranscript,'sample',exoneratehints,"false",emptyAugPlaceholder)
        Channel.empty().mix(augustus.out).flatten().toList().set { denovoAugFile }
      } else {
        Channel.empty().set { denovoAugFile }
      }
 } else {
     Channel.empty().set { denovoAugFile }
 }
 
 //make external hints channel
 emptyPlaceholder = Channel.fromPath(file("${params.output}/emptyPlaceHolder.txt"))
 Channel
  .empty()
  .mix(referenceRnaFile, rnaFile, referenceProteinFile, proteinFile, emptyPlaceholder)
  .flatten()
  .toList()
  .set { externalAnnotFiles }

 if (annotationalgos =~ /augustus_pretrained/) {
    pretrainedexoneratehints = file("${params.output}/emptyProteinPlaceHolder.txt")
    pretrainedtranscript = file("${params.output}/emptyTranscriptPlaceHolder.txt")
    augustus_pretrained(maskedGenome,pretrainedtranscript,params.pretrainedAugustusSpecies,pretrainedexoneratehints,"true",externalAnnotFiles)
    Channel.empty().mix(augustus_pretrained.out).flatten().toList().set { pretrainedAugFile }
 } else {
    Channel.empty().set { pretrainedAugFile }
 }

 if (annotationalgos =~ /Helixer/) {
    Helixer(maskedGenome,params.speciesScientificName,params.helixerModel,params.size)
    Channel.empty().mix(Helixer.out.gff3).flatten().toList().set { helixerFile }
    if (annotationalgos =~ /helixer_trained_augustus/) {
        helixer_trained_augustus(maskedGenome,helixerFile,params.lineage,externalAnnotFiles)
        Channel.empty().mix(helixer_trained_augustus.out).flatten().toList().set { helixtrainedAugFile }
    } else {
        Channel.empty().set { helixtrainedAugFile }
    }
 } else {
    Channel.empty().set { helixerFile }
    Channel.empty().set { helixtrainedAugFile }
 }
 
 if ( params.gtffile ) {
   if (annotationalgos =~ /related_species_augustus/) {
      related_species_augustus(maskedGenome,Channel.fromPath(params.gtffile),'sample',Channel.fromPath(params.fafile),"false",externalAnnotFiles)
      Channel.empty().mix(related_species_augustus.out).flatten().toList().set { relatedAugFile }
   } else {
      Channel.empty().set { relatedAugFile }
   }
   if (( params.fafile ) && (annotationalgos =~ /Liftoff/)) {
      lift = liftoff(maskedGenome,Channel.fromPath(params.fafile),Channel.fromPath(params.gtffile),Channel.from(params.gaps))
      Channel.empty().mix(lift).flatten().toList().set { liftFile }
      
      if (annotationalgos =~ /liftoff_trained_augustus/) {
        liftoff_trained_augustus(maskedGenome,liftFile,params.lineage,pretrainedtranscript)
        Channel.empty().mix(liftoff_trained_augustus.out).flatten().toList().set { externalAnnotFiles }
       } else {
        Channel.empty().set { lifttrainedAugFile }
       }
   } else {
      Channel.empty().set { lifttrainedAugFile }
      Channel.empty().set { liftFile }
   }
 } else {
   Channel.empty().set { liftFile }
   Channel.empty().set { relatedAugFile }
   Channel.empty().set { lifttrainedAugFile }
 }
 
 //making channels
 Channel
  .empty()
  .mix(denovoAugFile, helixerFile, pretrainedAugFile, relatedAugFile, liftFile, lifttrainedAugFile, helixtrainedAugFile)
  .flatten()
  .toList()
  .set { geneAnnotFiles }
 Channel
  .empty()
  .mix(externalAnnotFiles, geneAnnotFiles, trnascan.out)
  .flatten()
  .toList()
  .set { filesToFilter }
  
 CombineAndFilter(maskedGenome,filesToFilter,params.lineage,params.size)
 entap(maskedGenome.combine(CombineAndFilter.out.gtfLenient).combine(entapDB))
 combinestructwfunct(entap.out.combine(CombineAndFilter.out.gtfLenient).combine(Channel.from(params.speciesScientificName)).combine(Channel.from(params.lineage)).combine(maskedGenome))
}
