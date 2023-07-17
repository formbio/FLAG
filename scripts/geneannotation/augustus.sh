#!/bin/bash
#exonerate_p2g.sh

usage() {
  echo "-h Help documentation for tasser.sh"
  echo "-s  --input species"
  echo "-g  --input genome file"
  echo "-u  --UTR on or off"
  echo "Example: bash augustus.sh -s human -g GRCh38_latest_genomic.fna -u on"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :s:g:u:h opt
do
    case $opt in
        s) species=$OPTARG;;
        g) genome=$OPTARG;;
        u) utr=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${species} ]] || [[ -z ${genome} ]]
then
    usage
fi
if [[ -z ${utr} ]]
then
    utr='on'
fi

augustus --species=$species --UTR=$utr $genome
