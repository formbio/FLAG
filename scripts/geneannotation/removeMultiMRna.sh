#!/bin/bash
#evidencemodeler.sh

usage() {
  echo "-h Help documentation for liftoff.sh"
  echo "-i  --input genome file in fasta format"
  echo "Example: bash removeMultiMRna.sh -i annot.gff3"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:h opt
do
    case $opt in
        i) input=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${input} ]]
then
    usage
fi

awk '/##/{n++}{print >"annotSplit"  n ".gff3" }' "$input"
rm $input

for FILE in ./annotSplit*.gff3; do
    count=$(( $count + 1 ))
done


for (( i=1; i<=$count; i++ )); do
    file="annotSplit${i}.gff3"
    echo "$file"
    set numMRNA=0
    while IFS= read line; do
        if [[ "$line" == *"mRNA"* ]]; then
            numMRNA=$(( $numMRNA + 1 ))
        fi
        if [[ $numMRNA -le 1 ]]; then
            echo -e "$line" >> noDupsAnnot.gff3
        else
            break
        fi
    done < "$file"
    unset numMRNA
    rm $file
done
