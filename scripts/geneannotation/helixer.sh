#!/bin/bash
#evidencemodeler.sh

usage() {
  echo "-h Help documentation for liftoff.sh"
  echo "-i  --input genome file in fasta format"
  echo "-n  --species name ie Sminthopsis_crassicaudata"
  echo "-m  --model"
  echo "-s  --size (small or normal)"
  echo "Example: bash helixer.sh -i genome.fa"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:n:m:s:h opt
do
    case $opt in
        i) genome=$OPTARG;;
        n) name=$OPTARG;;
        m) model=$OPTARG;;
        s) size=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${genome} ]]
then
    usage
fi
if [[ -z ${model} ]]
then
    model="vertebrate"
fi
if [[ -z ${size} ]]
then
    size="normal"
fi
mv $genome genome.fa

baseDir="`dirname \"$0\"`"


cp ${baseDir}/workflow_helixer_config.yaml .
#running it through a config. It has problems parsing args when run through nextflow
if [[ "$size" == "small" ]]; then
    sed -i "s/106920/10692/g" workflow_helixer_config.yaml
    sed -i "s/160380/16038/g" workflow_helixer_config.yaml
    sed -i "s/213840/21384/g" workflow_helixer_config.yaml
fi

sed -i "s/vertebrate/${model}/g" workflow_helixer_config.yaml
sed -i "s/Species_name/${name}/g" workflow_helixer_config.yaml

echo "My updated config.yaml"
cat workflow_helixer_config.yaml
cp workflow_helixer_config.yaml /data/workflow_helixer_config.yaml

Helixer.py --config-path /data/workflow_helixer_config.yaml --fasta-path genome.fa --gff-output-path Helixer.gff3
