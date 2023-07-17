#!/bin/bash
#splign.sh

usage() {
  echo "-h Help documentation for blastn_splign.sh, this is mainly for miRNA"
  echo "-g  --input genome file in fasta format"
  echo "-r  --input cDNA or rna file in fasta format"
  echo "-t  --threads"
  echo "Example: bash splign.sh -g GCF_000001905.1_Loxafr3.0_genomic.fna -p GCF_000001905.1_Loxafr3.0_rna.fna"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :g:r:t:h opt
do
    case $opt in
        g) genome=$OPTARG;;
        r) input=$OPTARG;;
        t) threads=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${genome} ]] || [[ -z ${input} ]]
then
    usage
fi
if [[ -z $threads ]]
then
    threads=`nproc`
fi

mv ${genome} genome.fa
mv ${input} cdna.fa

# Create LDS index that splign will use to access your FASTA sequences
splign -mklds .


# Generate preliminary cDNA-to-genomic alignments using compart. Genomes dont need to be masked
makeblastdb -dbtype nucl -parse_seqids -in cdna.fa
makeblastdb -dbtype nucl -parse_seqids -in genome.fa

tblastn -query cdna.fa -db genome.fa -num_threads $threads -outfmt 6  | sort -k 2,2 -k 1,1 > blast.hit

# Run splign with the index and the files generated above
splign -ldsdir . -hits blast.hit -asn splign.asn > splign.out

if [ "$input" == "reference_trinity.fasta" ] || [ "$input" == "formatted_referenceRNA.fa" ]
then
    annotwriter -i splign.asn -format gff3 -o reference_splign.gff3
    sed -i -e 's/RefSeq/HOMOLOGY_RefSeq/g' reference_splign.gff3
    #cat reference_spligninter.gff3 | awk '{ if ($2 == ".") {$2 = "HOMOLOGY_RefSeq"}; print }' > reference_splign.gff3
else
    annotwriter -i splign.asn -format gff3 -o splign.gff3
    #cat spligninter.gff3 | awk '{ if ($2 == ".") {$2 = "RefSeq"}; print }' > splign.gff3
fi
