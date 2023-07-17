#!/bin/bash
#Splitting the fasta file into multiple fasta files based off of proteins
usage() {
  echo "-h Help documentation for combiner.sh"
  echo "-i  --input files"
  echo "-p  --type (exonerate or prosplign)"
  echo "Example: bash fasta_splitter.sh -i seq.fasta"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:p:h opt
do
    case $opt in
        i) fa=$OPTARG;;
        p) param=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
count=0
if [[ $param == "exonerate" ]]; then
    for FILE in ./*exonerate*.out; do
        count=$(( $count + 1 ))
    done
    for (( i=1; i<=$count; i++ )); do
        cat exonerate${i}.out >> exonerate.out
        rm exonerate${i}.out
    done
    
elif [[ $param == "rexonerate" ]]; then
    for FILE in ./*exonerate*.out; do
        count=$(( $count + 1 ))
    done
    for (( i=1; i<=$count; i++ )); do
        cat reference_exonerate${i}.out >> reference_exonerate.out
        rm reference_exonerate${i}.out
    done
elif [[ $param == "genomethreader" ]]; then
    for FILE in ./*genomethreader*.gff3; do
        count=$(( $count + 1 ))
    done
    for (( i=1; i<=$count; i++ )); do
        cat genomethreader${i}.gff3 >> genomethreader.gff3
        rm genomethreader${i}.gff3
    done
elif [[ $param == "rgenomethreader" ]]; then
    for FILE in ./*genomethreader*.gff3; do
        count=$(( $count + 1 ))
    done
    for (( i=1; i<=$count; i++ )); do
        cat reference_genomethreader${i}.gff3 >> reference_genomethreader.gff3
        rm reference_genomethreader${i}.gff3
    done
elif [[ $param == "prosplign" ]]; then
    for FILE in ./*prosplign*.gff3; do
        count=$(( $count + 1 ))
    done
    for (( i=1; i<=$count; i++ )); do
        cat prosplign${i}.gff3 >> prosplign.gff3
        rm prosplign${i}.gff3
    done
    
elif [[ $param == "rprosplign" ]]; then
    for FILE in ./*prosplign*.gff3; do
        count=$(( $count + 1 ))
    done
    for (( i=1; i<=$count; i++ )); do
        cat reference_prosplign${i}.gff3 >> reference_prosplign.gff3
        rm reference_prosplign${i}.gff3
    done
elif [[ "$param" == "genome2protein" ]]; then
    for FILE in ./referenceProtein*.fa; do
        count=$(( $count + 1 ))
        echo "$count"
    done
    for (( i=1; i<=$count; i++ )); do
        cat referenceProtein${i}.fa >> referenceProtein.fa
        rm referenceProtein${i}.fa
    done
elif [[ "$param" == "sequencerreads" ]]; then
    for FILE in ./*.fasta; do
        cat $FILE >> basecall.fasta
        rm $FILE
    done
    for FILE in ./*.fa; do
        cat $FILE >> basecall.fa
        rm $FILE
    done
    for FILE in ./*.fastq; do
        cat $FILE >> basecall.fastq
        rm $FILE
    done
    for FILE in ./*.fq; do
        cat $FILE >> basecall.fq
        rm $FILE
    done
    #delete empty files
    find . -type f -empty -delete
elif [[ "$param" == "masker" ]]; then
    for FILE in ./splitGenome_*.fa*; do
        count=$(( $count + 1 ))
        echo "$count"
    done
    for (( i=1; i<=$count; i++ )); do
        cat splitGenome_${i}_*.fa* >> maskedGenome.fa
        rm splitGenome_${i}_*.fa*
    done
fi
