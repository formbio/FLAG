#!/bin/bash
#Splitting the fasta file into multiple fasta files based off of proteins
usage() {
  echo "-h Help documentation for tasser.sh"
  echo "-i  --input tar.gz file"
  echo "-n  --genome name"
  echo "Example: bash fasta_splitter.sh -i seq.fasta"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:n:h opt
do
    case $opt in
        i) tarred=$OPTARG;;
        n) name=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $tarred ]]
then
    usage
fi
if [[ -z $name ]]
then
    name="genome"
fi

mkdir -p ref
tar -I pigz -xvf ${tarred} --no-same-owner --strip-components=1 -C ref

mv ref/genome.fa ${name}.fa
