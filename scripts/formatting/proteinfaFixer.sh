#!/bin/bash
#transdecoder.sh

usage() {
  echo "-h Help documentation for liftoff.sh"
  echo "-i  --input genome file in fasta format"
  echo "Example: bash transdecoder.sh -i transcript.fa -t threads"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:h opt
do
    case $opt in
        i) protein=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${protein} ]]
then
    usage
fi
#number=1000000
#>gi|6679997|ref|
declare -i v=1000000
awk 'BEGIN { v=1000000 } />/{sub(">", ">gi|"++v"|ref|")}1' $protein > proteinFinal.fa
