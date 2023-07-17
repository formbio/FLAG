#!/bin/bash
#Splitting the fasta file into multiple fasta files based off of proteins
usage() {
  echo "-h Help documentation for tasser.sh"
  echo "-g  --input fasta genome"
  echo "-a  --input annotation gff3"
  echo "Example: bash fasta_splitter.sh -i seq.fasta"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :g:a:h opt
do
    case $opt in
        g) genome=$OPTARG;;
        a) annotation=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $genome ]] || [[ -z $annotation ]]
then
    usage
fi

grep ">" $genome > genomeChrList.txt
rm $genome

while IFS= read line; do
    word=`echo "$line" | sed 's/>//g' | awk '{ print $1 }' | sed 's/ //g'`
    echo ${word} >> genomeChrListSimple.txt
done < "genomeChrList.txt"

awk '/\tgene\t/{n++}{print >"annotSplit"  n ".gtf" }' "$annotation"
rm $annotation

count=0
for FILE in ./annotSplit*.gtf; do
    count=$(( $count + 1 ))
done

for (( i=1; i<=$count; i++ )); do
    file="annotSplit${i}.gtf"
    #echo "${file}"
    currentHead=`head -n 1 ${file}`
    currentChr=`echo "${currentHead}" | awk '{ print $1 }'`
    currentProgram=`echo "${currentHead}" | awk '{ print $2 }'`
    currentStart=`echo "${currentHead}" | awk '{ print $4 }'`
    currentEnd=`echo "${currentHead}" | awk '{ print $5 }'`
    echo -e "${currentChr}\t${currentStart}\t${currentEnd}\t${currentProgram}\t${file}" >> intermediateUnsorted.txt
done

sort -V intermediateUnsorted.txt > intermediateSorted.txt

while IFS= read line; do
    word=`echo "$line"`
    grep "${line}" intermediateSorted.txt > currentFiles.txt
    while IFS= read line; do
        fileName=`echo "${line}" | awk '{ print $5 }'`
        cat $fileName >> sorted.gff3
        rm $fileName
    done < "currentFiles.txt"
done < "genomeChrListSimple.txt"
