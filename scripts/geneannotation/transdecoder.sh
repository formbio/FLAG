#!/bin/bash
#transdecoder.sh

usage() {
  echo "-h Help documentation for liftoff.sh"
  echo "-i  --input transcript file in fasta format"
  echo "-a  --input annotation file in gtf format"
  echo "-g  --input genome file in fasta format"
  echo "-t  --threads"
  echo "Example: bash transdecoder.sh -i transcript.fa -t threads"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:a:g:t:h opt
do
    case $opt in
        i) transcript=$OPTARG;;
        a) annotation=$OPTARG;;
        g) genome=$OPTARG;;
        t) threads=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${transcript} ]]
then
    usage
fi

if [[ -z ${threads} ]]
then
    threads=16
fi

mv ${transcript} transcripts.fasta
if [ -z $genome ] || [ -z $annotation ]; then
    #first step
    TransDecoder.LongOrfs -t transcripts.fasta

    #optional 2nd steps of doing a blast search and a pfam search (in this case we are making it manditory)
    blastp -query transcripts.fasta.transdecoder_dir/longest_orfs.pep -db /seq/RNASEQ/DBs/TRINOTATE_RESOURCES/TRINOTATE_V3/uniprot_sprot.pep -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads $threads > blastp.outfmt6

    hmmscan --cpu $threads --domtblout pfam.domtblout /seq/RNASEQ/DBs/TRINOTATE_RESOURCES/TRINOTATE_V3/Pfam-A.hmm transcripts.fasta.transdecoder_dir/longest_orfs.pep

    #step3 and final step
    TransDecoder.Predict -t transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6

    if [ "$transcript" == "reference_trinity.fasta" ]
    then
        mv transcripts.fasta.transdecoder.gff3 reference_transdecoder.gff3
        sed -i -e 's/transdecoder/HOMOLOGY_transdecoder/g' reference_transdecoder.gff3
    else
        mv transcripts.fasta.transdecoder.gff3 transdecoder.gff3
    fi
else
    mv ${genome} genome.fasta
    mv ${annotation} transcripts.gtf

    #Construct the transcript fasta file using the genome and the transcripts.gtf file like so:
    perl /opt/TransDecoder-TransDecoder-v5.5.0/util/gtf_genome_to_cdna_fasta.pl transcripts.gtf genome.fasta > transcripts.fasta

    #Next, convert the transcript structure GTF file to an alignment-GFF3 formatted file (this is done only because our processes operate on gff3 rather than the starting gtf file - nothing of great consequence).
    perl /opt/TransDecoder-TransDecoder-v5.5.0/util/gtf_to_alignment_gff3.pl transcripts.gtf > transcripts.gff3

    #Now, run the process described above to generate your best candidate ORF predictions:
    #first step
    TransDecoder.LongOrfs -t transcripts.fasta

    #optional 2nd steps of doing a blast search and a pfam search (in this case we are making it manditory)
    blastp -query transcripts.fasta.transdecoder_dir/longest_orfs.pep -db /seq/RNASEQ/DBs/TRINOTATE_RESOURCES/TRINOTATE_V3/uniprot_sprot.pep -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads $threads > blastp.outfmt6

    hmmscan --cpu $threads --domtblout pfam.domtblout /seq/RNASEQ/DBs/TRINOTATE_RESOURCES/TRINOTATE_V3/Pfam-A.hmm transcripts.fasta.transdecoder_dir/longest_orfs.pep

    #step3 and final step
    TransDecoder.Predict -t transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6

    And finally, generate a genome-based coding region annotation file:
    perl /opt/TransDecoder-TransDecoder-v5.5.0/util/cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3

    if [ "$transcript" == "reference_trinity.fasta" ] || [ "$transcript" == "formatted_referenceRNA.fa" ]
    then
        mv transcripts.fasta.transdecoder.genome.gff3 reference_transdecoder.gff3
        sed -i -e 's/transdecoder/HOMOLOGY_transdecoder/g' reference_transdecoder.gff3
    else
        mv transcripts.fasta.transdecoder.genome.gff3 transdecoder.gff3
    fi
fi
rm transcripts.fasta.transdecoder.gff3
