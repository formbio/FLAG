#!/bin/bash
#augustus_full.sh

usage() {
  echo "-h Help documentation for tasser.sh"
  echo "-g  --input genome file"
  echo "-p  --input protein file"
  echo "-s  --input species"
  echo "-n  --number of cpus"
  echo "Example: bash augustus.sh -s human -g GRCh38_latest_genomic.fna -u on"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :g:p:s:n:h opt
do
    case $opt in
        g) genome=$OPTARG;;
        p) protein=$OPTARG;;
        s) species=$OPTARG;;
        n) cpu=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${protein} ]] || [[ -z ${genome} ]]
then
    usage
fi
if [[ -z ${species} ]]
then
    species='current'
fi
if [[ -z ${cpu} ]]
then
    cpu=16
fi

perl /root/augustus/scripts/autoAug.pl --genome=$genome --species=$species --trainingset=$protein --cpus=$cpu --workingdir=/data/ -v

#/data/autoAug/autoAugPred_abinitio/shells/aug1
rm /data/autoAug/autoAugPred_abinitio/shells/shellForAug

#for file in /data/autoAug/autoAugPred_abinitio/shells/*
#do
#  "$file"
#done
touch commands.txt
for FILE in aug*; do
    cat $FILE >> commands.txt
    echo -e "" >> commands.txt
done

parallel < commands.txt

perl /root/augustus/scripts/autoAug.pl --genome=$genome --species=$species --trainingset=$protein --cpus=$cpu --workingdir=/data/ -v --useexisting

cp /data/autoAug/autoAugPred_abinitio/predictions/* .

#convert gtf file to gff3 for evm
#gtf2gff.pl < augustus.gtf --out=augustus.gff3 --gff3

#if [ "$protein" == "referenceProtein.fa" ]; then
#    mv augustus.gff3 reference_augustus.gff3
#    sed -i -e 's/AUGUSTUS/HOMOLOGY_AUGUSTUS/g' reference_augustus.gff3
#fi

if [ "$protein" == "referenceProtein.fa" ]; then
    mv augustus.gtf reference_augustus.gtf
    sed -i -e 's/AUGUSTUS/HOMOLOGY_AUGUSTUS/g' reference_augustus.gtf
fi
