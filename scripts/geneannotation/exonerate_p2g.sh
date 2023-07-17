#!/bin/bash
#exonerate_p2g.sh

usage() {
  echo "-h Help documentation for tasser.sh"
  echo "-p  --input protein file"
  echo "-g  --input genome file"
  echo "-n  --number of chunks"
  echo "-t  --threads"
  echo "Example: bash exonerate_p2g.sh -p GCF_000001905.1_Loxafr3.0_protein.faa -g GCF_000001905.1_Loxafr3.0_genomic.fna -t 24"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :p:g:n:t:h opt
do
    case $opt in
        p) protein=$OPTARG;;
        g) genome=$OPTARG;;
        n) numchunks=$OPTARG;;
        t) threads=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${protein} ]] || [[ -z ${genome} ]]
then
    usage
fi
if [[ -z ${numchunks} ]]
then
    numchunks=30
fi
if [ "$protein" == "formatted_referenceProtein.fa" ]; then
    numchunks=1000
fi
cp $protein protein.fa
#numchunks=$(grep -o '>' protein.fa | wc -l)
numchunks=$threads
if [[ -z ${threads} ]]
then
    # here the protein is the query and the genome is the target. The protein is being mapped to the genome
    #exonerate --model protein2genome $protein $genome --showtargetgff TRUE --softmasktarget TRUE > exonerate.out
    exonerate --model protein2genome $protein $genome --showtargetgff TRUE --softmasktarget TRUE > exonerate.out
else
    touch commands.txt
    touch exonerate.out
    for (( i=1; i<=$numchunks; i++ )); do
        echo -e "exonerate --model protein2genome ${protein} ${genome} --showtargetgff TRUE --softmasktarget TRUE --querychunkid ${i} --querychunktotal ${numchunks} > exonerate${i}.out" >> commands.txt
    done
    parallel --jobs $threads --memfree 60G --retries 300 < commands.txt

    for (( i=1; i<=$numchunks; i++ )); do
        if [ $i -gt 1 ]; then
            sed -i '2d' exonerate${i}.out
        fi
        cat exonerate${i}.out >> exonerate.out
        rm exonerate${i}.out
    done
fi
if [[ "$genome" == *splitGenome* ]]; then
    back=${genome%.fasta}
    count=${back#*_}
    if [ "$protein" == "formatted_referenceProtein.fa" ]
    then
        mv exonerate.out reference_exonerate${count}.out
    else
        mv exonerate.out exonerate${count}.out
    fi
else
    if [ "$protein" == "formatted_referenceProtein.fa" ]
    then
        mv exonerate.out reference_exonerate.out
    fi
fi

