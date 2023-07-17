#!/bin/bash
#windowmasker.sh

usage() {
  echo "-h Help documentation for tasser.sh"
  echo "-g  --Unmasked genome"
  echo "-t  --threads"
  echo "Example: bash windowmasker.sh -g genomeUnmasked.fna -t 30 -m 120G"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :g:t:h opt
do
    case $opt in
        g) genome=$OPTARG;;
        t) threads=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${genome} ]]
then
    usage
fi
if [[ -z $threads ]]
then
    threads=`nproc`
fi

mv $genome genome.fa

declare -i count=$(grep -o ">" genome.fa | wc -l)
if [ "$count" != "1" ]; then
    awk '/>/{n++}{print >"___"  n ".fasta" }' "genome.fa"
    rm genome.fa
fi

touch commands.txt
for FILE in ./___*.fasta; do
    back=${FILE%.fasta}
    count=${back:5}
    touch command_${count}.sh
    
    # First step of windowmasker
    echo -e "windowmasker -mk_counts -in ${FILE:2} -checkdup true -out wm_counts_${count}" >> command_${count}.sh
    #windowmasker -mk_counts -in $genome -out intermediate.txt

    # Seconds step of windowmasker
    echo -e "windowmasker -ustat wm_counts_${count} -in ${FILE:2} -outfmt fasta -out maskedGenome_${count}.fa -dust true" >> command_${count}.sh
    #windowmasker -ustat intermediate.txt
    
    echo -e "bash command_${count}.sh" >> commands.txt
done

parallel --jobs $threads < commands.txt

touch maskedGenome.fa
for (( i=1; i<=$count; i++ )); do
    echo $i
    cat maskedGenome_${i}.fa >> maskedGenome.fa
done


# First step of windowmasker
##windowmasker -mk_counts -in $genome -checkdup true -out wm_counts
#windowmasker -mk_counts -in $genome -out intermediate.txt

# Seconds step of windowmasker
##windowmasker -ustat wm_counts -in $genome -outfmt fasta -out maskedGenome.fa -dust true
#windowmasker -ustat intermediate.txt
