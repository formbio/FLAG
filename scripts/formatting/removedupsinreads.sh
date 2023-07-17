#!/bin/bash
#Splitting the fasta file into multiple fasta files based off of proteins
usage() {
  echo "-h Help documentation for tasser.sh"
  echo "-i  --input fasta sequence"
  echo "Example: bash fasta_splitter.sh -i seq.fasta"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:h opt
do
    case $opt in
        i) fq=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $fq ]]
then
    usage
fi

prefix="${fq%.*}"
string=""
set line_num = 0

count=0
total=0
declare -i number=$(( $num + 0 ))
while IFS= read line; do
    line_num=$(( $line_num + 1 ))
    if [ $line_num == 1 ]; then
        newFile=${line:1}.fastq
    elif [ $line_num == 4 ]; then
        unset line_num
        set line_num = 0
    fi
        echo "$line" >> $newFile
done < "$fq"
rm $fq

find . -name "*.fastq" -type f -size -300 -delete
for f in *.fastq; do
    head -n 4 $f >> formatted_${prefix}.fastq
    rm $f
done
