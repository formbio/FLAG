#!/bin/bash
#evidencemodeler.sh

usage() {
  echo "-h Help documentation for liftoff.sh"
  echo "-i  --input annotation file in gff3 format"
  echo "Example: bash helixer.sh -i genome.fa"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:h opt
do
    case $opt in
        i) annotation=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${annotation} ]]
then
    usage
fi

echo "fixing helixer annotations in gff3 file: ${annotation}"
while IFS= read line; do
    third=$(echo $line | awk '{print $3}')
    second=$(echo $line | awk '{print $2}')
    if [[ "$second" == "Helixer" ]] && [[ "$third" == "gene" ]]; then
        suffix="${line##*ID=}"
        echo -e "${line};Name=${suffix}" >> gff3togtf_tmp.gff3
    else
        echo -e "$line" >> gff3togtf_tmp.gff3
    fi
done < "${annotation}"

echo "adding comments to annotation"
currentTime=`date`
sed -i "1 i\#!Form Bio annotation workflow version 2.0" gff3togtf_tmp.gff3
sed -i "1 i\#!Created with the Form Bio annotation workflow on: ${currentTime}" gff3togtf_tmp.gff3
sed -i "1 i\##gff-version 3" gff3togtf_tmp.gff3

echo "tidying up annotation"
gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o gff3togtf_tmp2.gff3  gff3togtf_tmp.gff3
mv gff3togtf_tmp2.gff3 $annotation
rm gff3togtf_tmp.gff3

prefix="${annotation%.gff3}"

echo "converting annotation from gff3 to gtf"
agat_convert_sp_gff2gtf.pl --gff ${annotation} -o ${prefix}.gtf
