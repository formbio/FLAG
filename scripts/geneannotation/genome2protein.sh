#!/bin/bash
#genome2proteinfa.sh

usage() {
  echo "-h Help documentation for blastx.sh"
  echo "-i  --input dna file"
  echo "-d  --database"
  echo "-a  --algo"
  echo "-t  --number of threads"
  echo "Example: bash genome2proteinfa.sh -i GCF_000001905.1_Loxafr3.0_genome.fna -d swissprot"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:d:a:t:h opt
do
    case $opt in
        i) input=$OPTARG;;
        d) database=$OPTARG;;
        a) algo=$OPTARG;;
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

mv $input genome.fa

#first we split the genome up into separate reads to limit RAM usage
declare -i count=$(grep -o ">" genome.fa | wc -l)
awk '/>/{n++}{print >"___"  n ".fasta" }' "genome.fa"
rm genome.fa

touch blastx.out
for (( i=1; i<=$count; i++ )); do
    if [ "$algo" == "diamond" ]; then
        diamond blastx -d ${dbname} -q ___${i}.fasta -p${threads} --sensitive >> blastx${i}.out
    else
        blastx -query ___${i}.fasta -db ${dbname} -num_threads $threads >> blastx${i}.out
    fi
    cat blastx${i}.out >> blastx.out
done

newFile="referenceProteinAccessions.txt"
touch $newFile
if [ "$algo" == "diamond" ]; then
    while IFS= read line; do
        set word_num = 0
        word_num=$(( $word_num + 0 ))
        for word in $line; do
            word_num=$(( $word_num + 1 ))
            if [  "$word_num" = "2"  ]; then
                echo -e "${word}" >> $newFile
            fi
        done
        unset word_num
    done < "blastx.out"
else
    set line_num = 0
    while IFS= read line; do
        if [[ "$line" == ">"* ]]; then
            for word in $line; do
                if [ "${word:0:1}" == ">" ]; then
                    echo -e "${word:1}" >> $newFile
                fi
            done
        fi
    done < "blastx.out"
fi

blastdbcmd -entry_batch referenceProteinAccessions.txt -db ${dbname} > referenceProtein.fa

if [[ "$input" == *splitGenome* ]]; then
    backN=${input%.fasta}
    countN=${backN#*_}
    mv referenceProtein.fa referenceProtein${countN}.fa
fi
