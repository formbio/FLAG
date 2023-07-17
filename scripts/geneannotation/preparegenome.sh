#!/bin/bash
#preparegenome.sh

usage() {
  echo "-h Help documentation for liftoff.sh"
  echo "-i  --input genome file in fasta format"
  echo "Example: bash preparegenome.sh -i transcript.fa"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:h opt
do
    case $opt in
        i) genome=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${genome} ]]
then
    usage
fi

mv $genome genome.fa

samtools faidx genome.fa

java -jar /usr/local/bin/picard.jar CreateSequenceDictionary -R genome.fa -O genome.dict

mkdir genome

mv genome.fa genome.dict genome.fa.fai genome

tar cfz genomefa.tar.gz genome
