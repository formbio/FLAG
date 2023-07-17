#!/bin/bash
#exonerate_p2g.sh

usage() {
  echo "-h Help documentation for tasser.sh"
  echo "-i  --input entap tsv"
  echo "-a  --input annotation gtf file"
  echo "-s  --species name"
  echo "Example: bash exonerate_p2g.sh -p GCF_000001905.1_Loxafr3.0_protein.faa -g GCF_000001905.1_Loxafr3.0_genomic.fna -t 24"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:a:s:h opt
do
    case $opt in
        i) entaptsv=$OPTARG;;
        a) annotation=$OPTARG;;
        s) speciesName=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${entaptsv} ]] || [[ -z ${annotation} ]]
then
    usage
fi
if [[ -z ${speciesName} ]]
then
    speciesName=Sample_species
fi

baseDir="`dirname \"$0\"`"

cat final_annotations_lvl0.tsv | awk -F '\t' '{ print $1"\t"$25"\t"$28"\t"$29"\t"$32"\t"$2 }' >> simplifiedEntap.tsv

sed -i 's/ /_/g' simplifiedEntap.tsv

sed -i '/^#/d' "$annotation"

awk '/\tgene\t/ || /\ttRNA\t/{n++}{print >"annotSplit"  n ".gtf" }' "$annotation"

rm $annotation

for file in ./annotSplit*.gtf; do
    echo -e "bash ${baseDir}/renameAnnots.sh -a ${file} -s ${speciesName}" >> renameParallel.txt
done

parallel < renameParallel.txt > renameParallel.log

count=0
for file in ./updatedAnnot_*.gtf; do
    count=$(( $count + 1 ))
done

for (( i=1; i<=$count; i++ )); do
    cat updatedAnnot_${i}.gtf >> finalAnnotation.gtf
    rm updatedAnnot_${i}.gtf
done

echo "adding comments to annotation"
currentTime=`date`
sed -i "1 i\#!Species: ${speciesName}" finalAnnotation.gtf
sed -i "1 i\#!Form Bio annotation workflow version 2.0" finalAnnotation.gtf
sed -i "1 i\#!Created with the Form Bio annotation workflow on: ${currentTime}" finalAnnotation.gtf
sed -i "1 i\##gtf-version 3" finalAnnotation.gtf

mv finalAnnotation.gtf final${speciesName}.gtf

