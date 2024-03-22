#!/bin/bash
#splign.sh

usage() {
  echo "-h Help documentation for liftoff.sh"
  echo "-g  --input genome file in fasta format"
  echo "-r  --input cDNA or rna file in fasta format"
  echo "-n  --number to split"
  echo "-p  --parallel"
  echo "Example: bash splign.sh -g GCF_000001905.1_Loxafr3.0_genomic.fna -p GCF_000001905.1_Loxafr3.0_rna.fna"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :g:r:n:p:h opt
do
    case $opt in
        g) genome=$OPTARG;;
        r) input=$OPTARG;;
        n) numSeqsSplit=$OPTARG;;
        p) parallelNum=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${genome} ]] || [[ -z ${input} ]]
then
    usage
fi
if [[ -z ${numSeqsSplit} ]]
then
    numSeqsSplit="50000"
fi
if [[ -z ${parallelNum} ]]
then
    parallelNum="5"
fi

#checking if the fasta is a .gz
faEnd="${genome: -3}"
if [[ "${faEnd}" == ".gz" ]]; then
    cp $genome genome.fa.gz
    gunzip genome.fa.gz
else
    cp $genome genome.fa
fi
cp ${input} cdna.fa

seqkit split cdna.fa -s ${numSeqsSplit}
mv cdna.fa.split/* .
rm -rf cdna.fa.split
makeblastdb -dbtype nucl -parse_seqids -in genome.fa
numQueries=0
for file in cdna.part_*.fa; do
    numQueries=$(( $numQueries + 1 ))
    mkdir "${numQueries}_folder"
    padded_number=$(printf "%03d" "$numQueries")
    mv cdna.part_${padded_number}.fa ${numQueries}_folder/cdna.fa
    echo "cd ${numQueries}_folder" >> parallel_${padded_number}.txt
    echo "cp ../genome.fa* ." >> parallel_${padded_number}.txt
    echo "splign -mklds ." >> parallel_${padded_number}.txt
    echo "makeblastdb -dbtype nucl -parse_seqids -in cdna.fa" >> parallel_${padded_number}.txt
    #echo "makeblastdb -dbtype nucl -parse_seqids -in genome.fa" >> parallel_${padded_number}.txt
    echo "compart -qdb cdna.fa -sdb genome.fa > cdna.compartments" >> parallel_${padded_number}.txt
    echo "splign -ldsdir . -comps cdna.compartments -asn splign.asn > splign.out" >> parallel_${padded_number}.txt
    echo "annotwriter -i splign.asn -format gff3 -o splign.gff3" >> parallel_${padded_number}.txt
    echo "splign -mklds ." >> parallel_${padded_number}.txt
    echo "cd .." >> parallel_${padded_number}.txt
    echo "rm -rf ${numQueries}_folder" >> parallel_${padded_number}.txt
    echo "bash parallel_${padded_number}.txt" >> parallel_commands.txt
done



parallel -j $parallelNum < parallel_commands.txt

numFolders=0
for folder in *_folder; do
    numFolders=$(( $numFolders + 1 ))
    sed -i "s/ID=aln/ID=aln_${numFolders}_/g" ${numFolders}_folder/splign.gff3
    cat ${numFolders}_folder/splign.gff3 >> splign.gff3
done

## Create LDS index that splign will use to access your FASTA sequences
#splign -mklds .
#
#
## Generate preliminary cDNA-to-genomic alignments using compart. Genomes dont need to be masked
#makeblastdb -dbtype nucl -parse_seqids -in cdna.fa
#makeblastdb -dbtype nucl -parse_seqids -in genome.fa
#compart -qdb cdna.fa -sdb genome.fa > cdna.compartments
#
#
## Run splign with the index and the files generated above
#splign -ldsdir . -comps cdna.compartments -asn splign.asn > splign.out

if [ "$input" == "reference_trinity.fasta" ] || [ "$input" == "formatted_referenceRNA.fa" ]
then
#    annotwriter -i splign.asn -format gff3 -o reference_splign.gff3
    mv splign.gff3 reference_splign.gff3
    sed -i -e 's/RefSeq/HOMOLOGY_RefSeq/g' reference_splign.gff3
    #cat reference_spligninter.gff3 | awk '{ if ($2 == ".") {$2 = "HOMOLOGY_RefSeq"}; print }' > reference_splign.gff3
#else
#    annotwriter -i splign.asn -format gff3 -o splign.gff3
    #cat spligninter.gff3 | awk '{ if ($2 == ".") {$2 = "RefSeq"}; print }' > splign.gff3
fi
