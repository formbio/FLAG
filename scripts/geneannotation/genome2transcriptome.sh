#!/bin/bash
#genome2transcriptomefa.sh

usage() {
  echo "-h Help documentation for genome2transcriptomefa.sh"
  echo "-i  --input dna file"
  echo "-d  --database"
  echo "-t  --number of threads"
  echo "Example: bash genome2transcriptomefa.sh -i GCF_000001905.1_Loxafr3.0_genome.fna -d refseq_rna.tar.gz -t 31"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:d:t:h opt
do
    case $opt in
        i) input=$OPTARG;;
        d) database=$OPTARG;;
        t) threads=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${input} ]]
then
    usage
fi
if [[ -s ${database} ]]
then
    mkdir -p db
    tar xvfz ${database} --strip-components=1 -C db
    dbname=$(ls db/*db | cut -f 1 -d '.'|sort -u)
else
    echo "Missing Database File"
    usage
fi
if [[ -z ${threads} ]]
then
    threads=31
fi

blastn -query $input -db ${dbname} -num_threads $threads >> blastn.out

newFile="referenceRNAAccessions.txt"
touch $newFile
set line_num = 0
while IFS= read line; do
    if [[ "$line" == ">"* ]]; then
        for word in $line; do
            if [ "${word:0:1}" == ">" ]; then
                echo -e "${word:1}" >> $newFile
            fi
        done
    fi
done < "blastn.out"

blastdbcmd -entry_batch referenceRNAAccessions.txt -db ${dbname} > referenceRNA.fa
