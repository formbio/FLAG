#!/bin/bash
#Splitting the fasta file into multiple fasta files based off of proteins
usage() {
  echo "-h Help documentation for tasser.sh"
  echo "-i  --input fasta sequence"
  echo "-n  --number of sequences per fa"
  echo "-s  --size of files"
  echo "Example: bash fasta_splitter.sh -i seq.fasta"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:n:s:h opt
do
    case $opt in
        i) fa=$OPTARG;;
        n) num=$OPTARG;;
        s) size=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $fa ]]
then
    usage
fi
if [[ -z $num ]]
then
    num=40
fi
if [[ -z $size ]]
then
    size=99999999999999999999999999999
fi

mv $fa genomeIn.fa
declare -i count=$(grep -o ">" genomeIn.fa | wc -l)
if [ "$count" != "1" ]; then
    awk '/>/{n++}{print >"firstSplitGenome"  n ".fa" }' "genomeIn.fa"
    rm genomeIn.fa
fi

countCurrent=0
countTotal=0
declare -i fsize=$(( $size + 0 ))
for (( i=1; i<=$count; i++ )); do
    countCurrent=$(( $countCurrent + 1 ))
    FILENAME=firstSplitGenome${i}.fa
    declare -i FILESIZE=$(stat -c%s "$FILENAME")
    if [ $countCurrent == 1 ]; then
        countTotal=$(( $countTotal + 1 ))
        touch splitGenome_${countTotal}.fasta
    elif [ $countCurrent == ${num} ]; then
        countCurrent=0
    elif [ $FILESIZE -ge $fsize ]; then
        countCurrent=0
        countTotal=$(( $countTotal + 1 ))
        touch splitGenome_${countTotal}.fasta
    fi
    echo "this is the file size"
    echo "$FILESIZE"
    echo "this is the size limiter"
    echo "$size"
    cat firstSplitGenome${i}.fa >> splitGenome_${countTotal}.fasta
    rm firstSplitGenome${i}.fa
done
