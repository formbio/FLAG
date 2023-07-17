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
    numchunks=100
fi

mv $protein protein.fa
mv $genome genome.fa
#numchunks=$(grep -o '>' protein.fa | wc -l)
numchunks=$threads
if [[ -z ${threads} ]]
then
    # here the protein is the query and the genome is the target. The protein is being mapped to the genome
    #exonerate --model protein2genome $protein $genome --showtargetgff TRUE --softmasktarget TRUE > exonerate.out
    gth --genomic genome.fa --protein protein.fa -gff3out -v -o genomethreader.gff3 -skipalignmentout -paralogs
else
    num=numseq
    string=""
    set line_num = 0

    count=0
    total=0
    declare -i number=$(( $num + 0 ))
    while IFS= read line; do
        line_num=$(( $line_num + 1 ))
        if [[ "$line" == ">"* ]]; then
            first=${line%% *}
            if [ $num == 0 ]; then
                oldFile=$newFile
                newFile="${prefix}_orig.fasta"
            else
                count=$(( $count + 1 ))
                if [ $count -le $number ]; then
                    oldFile=$newFile
                    newFile=${total}.fasta
                else
                    count=1
                    total=$(( $total + 1 ))
                    oldFile=$newFile
                    newFile=${total}.fasta
                fi
            fi
            if [ "$line_num" != 1 ]; then
                echo "$string" >> $oldFile
                echo "$first" >> $newFile
            else
                echo "$first" > $newFile
            fi
            unset string
            set string=""
        else
            if [[ "$line" != "" ]]; then
                string="$string$line"
            fi
        fi
    done < "protein.fa"
    echo "$string" >> $newFile
    rm protein.fa

    #make commands
    count=0
    for FILE in ./*.fasta; do
        count=$(( $count + 1 ))
    done
    
    touch commands.txt
    touch genomethreader.gff3
    for (( i=1; i<=$count; i++ )); do
        file="${i}.fasta"
        echo -e "gth --genomic genome.fa --protein ${file} -gff3out -v -o genomethreader_${i}.gff3 -skipalignmentout -paralogs" >> commands.txt
    done

    parallel --jobs $threads --memfree 60G --retries 300 < commands.txt

    for (( i=1; i<=$count; i++ )); do
        cat genomethreader_${i}.gff3 >> genomethreader.gff3
        rm genomethreader_${i}.gff3
    done
fi

#remove extra garbage
grep -v '^#' genomethreader.gff3 > tmp_genomethreader.gff3
mv tmp_genomethreader.gff3 genomethreader.gff3
sed -i "1 i\##gff-version 3" genomethreader.gff3

if [[ "$genome" == *splitGenome* ]]; then
    back=${genome%.fasta}
    count=${back#*_}
    if [ "$protein" == "formatted_referenceProtein.fa" ]
    then
        mv genomethreader.gff3 reference_genomethreader${count}.gff3
    else
        mv genomethreader.gff3 genomethreader${count}.gff3
    fi
else
    if [ "$protein" == "formatted_referenceProtein.fa" ]
    then
        mv genomethreader.gff3 reference_genomethreader.gff3
    fi
fi

