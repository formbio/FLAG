#!/bin/bash
#tblastn_prosplign.sh

usage() {
  echo "-h Help documentation for liftoff.sh"
  echo "-g  --input genome file in fasta format"
  echo "-p  --input protein file in fasta format"
  echo "-n  --number of sequences per protein"
  echo "-t  --threads"
  echo "Example: bash tblastn_prosplign.sh -g GCF_000001905.1_Loxafr3.0_genomic.fna -p GCF_000001905.1_Loxafr3.0_genomic.faa -n 16"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :g:p:n:t:h opt
do
    case $opt in
        g) genome=$OPTARG;;
        p) protein=$OPTARG;;
        n) numseq=$OPTARG;;
        t) cpu=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${genome} ]] || [[ -z ${protein} ]]
then
    usage
fi
if [[ -z $cpu ]]
then
    cpu=`nproc`
fi

if [[ -z $numseq ]]
then
    numseq=30
fi

#mv ${genome} subj.fa
#make the subject/genome 2 column format
declare -i w=30000000
awk 'BEGIN { w=30000000 } />/{sub(">", ">gi|"++w"|ref|")}1' $genome > subj.fa
#make the query/protein 2 column format
declare -i v=1000000
awk 'BEGIN { v=1000000 } />/{sub(">", ">gi|"++v"|ref|")}1' $protein > query.fa
#mv ${protein} query.fa
rm $protein
rm $genome

num=numseq
fa="query.fa"
prefix="${fa%.*}"
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
done < "$fa"
echo "$string" >> $newFile
rm $fa

makeblastdb -dbtype nucl -in subj.fa
    
touch commands.txt

numQueries=0
for query in ./*.fasta; do
    numQueries=$(( $numQueries + 1 ))
    back=${query%.fasta}
    count=${back:2}
    touch command_${count}.sh

    echo -e "tblastn -query ${query:2} -db subj.fa -num_threads 1 -outfmt 6  | sort -k 2,2 -k 1,1 > blast_${count}.hit" >> command_${count}.sh

    echo -e "cat blast_${count}.hit | procompart -t > comp${count}" >> command_${count}.sh

    echo -e "prosplign -i comp${count} -fasta subj.fa,${query:2} -nogenbank -o prosplign${count}.asn -eo pro${count}.txt" >> command_${count}.sh
    echo -e "bash command_${count}.sh" >> commands.txt
done

parallel --jobs $cpu --memfree 150G --retries 500  < commands.txt

touch prosplign.asn
for (( i=1; i<=$numQueries; i++ )); do
    echo $i
    cat prosplign${i}.asn >> prosplign.asn
done

if [ "$protein" == "formatted_referenceProtein.fa" ]
then
    annotwriter -i prosplign.asn -format gff3 -o reference_prosplign.gff3
    sed -i -e 's/RefSeq/HOMOLOGY_Prosplign/g' reference_prosplign.gff3
    if [[ "$genome" == *splitGenome* ]]; then
        backN=${genome%.fasta}
        countN=${backN#*_}
        mv reference_prosplign.gff3 reference_prosplign${countN}.gff3
    fi
else
    annotwriter -i prosplign.asn -format gff3 -o prosplign.gff3
    sed -i -e 's/RefSeq/Prosplign/g' prosplign.gff3
    if [[ "$genome" == *splitGenome* ]]; then
        backN=${genome%.fasta}
        countN=${backN#*_}
        mv prosplign.gff3 prosplign${countN}.gff3
    fi
fi
