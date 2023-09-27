#!/bin/bash
#tRNAScan.sh

usage() {
  echo "-h Help documentation for liftoff.sh"
  echo "-g  --input genome file in fasta format"
  echo "-t  --threads"
  echo "Example: bash tRNAScan.sh -g genome.fa -t 48"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :g:t:h opt
do
    case $opt in
        g) genome=$OPTARG;;
        t) cpu=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${genome} ]]
then
    usage
fi
if [[ -z $cpu ]]
then
    cpu=`nproc`
fi

mv $genome genomeIn.fa

awk '/>/{n++}{print >"SplitGenome"  n ".fa" }' "genomeIn.fa"

touch commands.txt
numQueries=0
for query in ./SplitGenome*.fa; do
    numQueries=$(( $numQueries + 1 ))
    front=${query#./SplitGenome}
    num=${front%.fa}
    echo -e "tRNAscan-SE -E --gff ${num}.gff3 ${query}" >> commands.txt
done

parallel --jobs $cpu --memfree 20G --retries 500  < commands.txt
prefixAssembly="${genome%.fa*}"

touch tRNAScan.gff3
for (( i=1; i<=$numQueries; i++ )); do
    cat ${i}.gff3 >> tRNAScan_${prefixAssembly}.gff3
done

cp tRNAScan_${prefixAssembly}.gff3 tRNAScan.gff3
