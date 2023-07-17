#!/bin/bash
#removedupsinfa.sh

usage() {
  echo "-h Help documentation for liftoff.sh"
  echo "-i  --input file in fasta format"
  echo "Example: bash removedupsinfa.sh -i rna.fa"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:h opt
do
    case $opt in
        i) input=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${input} ]]
then
    usage
fi

prefix="${input%.fa*}"
formatted=formatted_${prefix}.fa

seqkit rmdup -s $input | seqkit rename - > $formatted

rm $input
