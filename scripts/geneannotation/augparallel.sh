#!/bin/bash
#tblastn_prosplign.sh

usage() {
  echo "-h Help documentation for liftoff.sh"
  echo "-t  --threads"
  echo "Example: bash tblastn_prosplign.sh -g GCF_000001905.1_Loxafr3.0_genomic.fna -p GCF_000001905.1_Loxafr3.0_genomic.faa -n 16"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :t:h opt
do
    case $opt in
        t) cpu=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options

if [[ -z $cpu ]]
then
    cpu=`nproc`
fi

touch commands.txt

for query in ./*.fasta; do
    numQueries=$(( $numQueries + 1 ))
    back=${query%.fasta}
    count=${back:2}
    
    echo -e "augustus --species=tdevil ${query} > annot_${count}.gff" >> commands.txt
done

parallel --jobs $cpu --memfree 10G --retries 500  < commands.txt
