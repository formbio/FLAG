#!/bin/bash
#liftoff.sh

usage() {
  echo "-h Help documentation for liftoff.sh"
  echo "-q  --query fa file"
  echo "-t  --target fa file"
  echo "-r  --reference annotated gff file"
  echo "-g  --gaps (true or false)"
  echo "Example: bash liftoff.sh -q GCF_000001905.1_Loxafr3.0_genomic.fna -t GCF_000001905.1_Loxafr3.0_genomic.fna -r GCF_000001905.1_Loxafr3.0_genomic.gff"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :q:t:r:g:h opt
do
    case $opt in
        q) query=$OPTARG;;
        t) target=$OPTARG;;
        r) reference=$OPTARG;;
        g) gaps=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${query} ]] || [[ -z ${target} ]] || [[ -z ${reference} ]]
then
    usage
fi
if [[ -z ${gaps} ]]
then
    gaps="False"
fi

#here query is the target fasta genome to lift genes to
#target is reference fasta genome to lift genes from
referenceEnd="${reference: -3}"
echo "end of query file is: ${queryEnd}. This is to check if its a gtf file"
if [[ "${referenceEnd}" == "gtf" ]]; then
    agat_convert_sp_gxf2gxf.pl -g $reference -o liftoff_query.gff3
    rm $reference
else
    mv $reference liftoff_query.gff3
fi

if [[ "$gaps" == "False" ]] || [[ "$gaps" == "false" ]]; then
    liftoff $query $target -g liftoff_query.gff3 -o liftoff.gff3
else
    liftoff $query $target -d 4 -g liftoff_query.gff3 -o liftoff.gff3
fi
