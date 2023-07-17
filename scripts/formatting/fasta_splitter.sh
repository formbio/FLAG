#!/bin/bash
#Splitting the fasta file into multiple fasta files based off of proteins
usage() {
  echo "-h Help documentation for tasser.sh"
  echo "-i  --input fasta sequence"
  echo "-n  --number of sequences per fa"
  echo "-s  --size of files"
  echo "-c  --check if NT sequence, AA sequence, or dont care anything. (NT, AA, any)"
  echo "Example: bash fasta_splitter.sh -i seq.fasta"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:n:s:c:h opt
do
    case $opt in
        i) fa=$OPTARG;;
        n) num=$OPTARG;;
        s) size=$OPTARG;;
        c) check=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $fa ]]
then
    usage
fi
if [[ -z $num ]]
then
    num=1
fi
if [[ -z $size ]]
then
    size=0
fi
if [[ -z $check ]]
then
    check="any"
fi

writehtml () {
    echo "<html>
 <head>
   <title>
   The input fasta file contained a character that is not allowed or the file is empty.
   </title>
 </head>
 
 <body>
   ${1}
 </body>
 </html>" > results.html
 rm *.fa
 rm *.fasta
}

if [ ! -s $fa ]; then
    writehtml "No sequence found. The fasta file is empty."
    exit 0
fi

totalFilSizeMB=`ls -s --block-size=1048576 $fa | cut -d' ' -f1`
echo "The total file size in MB is ${totalFilSizeMB}"
if [[ $(( $totalFilSizeMB )) -ge 4000 ]]; then
    echo "file very large. Doing a simple split"
    mv $fa input.fa
    csplit input.fa "/^>/" "{*}"
    for file in xx*; do
        mv $file ${file:2}.fasta
    done
    rm input.fa
elif [[ "$size" == "0" ]]; then
    prefix="${fa%.*}"
    string=""
    set line_num = 0

    sed -i 's/\r//g' $fa
    sed -i 's/\n//g' $fa

    count=0
    total=0
    declare -i number=$(( $num + 0 ))
    OLDIFS=$IFS; IFS=$'\n';
    for line in $(<$fa); do
        line_num=$(( $line_num + 1 ))
        if [[ "$line" == ">"* ]]; then
            first=${line%% *}
            if [ $num == 0 ]; then
                oldFile=$newFile
                newFile="${prefix}_orig.fasta"
            elif [ $num == 1 ]; then
                oldFile=$newFile
                tempString=${first:1}
#                mod=${tempString//./_}
                newFile=${tempString}.fasta
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
                if [[ "$check" == "NT" ]]; then
                    if ! [[ $line =~ ^[ANTGCUantcgu]+$ ]]; then
                        echo "does not match nucleotides"
                        echo "$line"
                        writehtml "This sequence is not a nucleotide sequence: ${line}"
                        exit 0
                    fi
                elif [[ "$check" == "AA" ]]; then
                    if ! [[ $line =~ ^[ARNDBCEQZGHILKMFPSTVWYZarndbceqzghilkmfpstvwyz]+$ ]]; then
                        echo "does not match amino acids"
                        echo "$line"
                        writehtml "This sequence is not a amino acid sequence: ${line}"
                        exit 0
                    fi
                elif [[ "$check" == "NTAA" ]]; then
                    if ! [[ $line =~ ^[ARNDBCEQZGHILKMFPSTVWYZUarndbceqzghilkmfpstvwyzu]+$ ]]; then
                        echo "does not match amino acids or nucleotides"
                        echo "$line"
                        writehtml "This sequence is not a amino acid or nucleotide sequence: ${line}"
                        exit 0
                    fi
                fi
                string="$string$line"
            fi
        fi
    done
    echo "$string" >> $newFile
    numFiles=$(ls *.fa* 2> /dev/null | wc -l)
    if [[ "$numFiles" != "1" ]]; then
        rm $fa
    fi
else
    mv $fa genomeIn.fa
    declare -i count=$(grep -o ">" genomeIn.fa | wc -l)
    if [ "$count" != "1" ]; then
        awk '/>/{n++}{print >"firstSplitGenome"  n ".fa" }' "genomeIn.fa"
        rm genomeIn.fa
    fi

    countCurrent=0
    countTotal=0
    declare -i fsize=$(( $size + 0 ))
    for (( i=1; i<=$count; i++ )); do
        countCurrent=$(( $countCurrent + 1 ))
        FILENAME=firstSplitGenome${i}.fa
        declare -i FILESIZE=$(stat -c%s "$FILENAME")
        if [ $countCurrent == 1 ]; then
            countTotal=$(( $countTotal + 1 ))
            touch splitGenome_${countTotal}.fasta
        elif [ $countCurrent == ${num} ]; then
            countCurrent=0
        elif [ $FILESIZE -ge $fsize ]; then
            countCurrent=0
            countTotal=$(( $countTotal + 1 ))
            touch splitGenome_${countTotal}.fasta
        fi
        echo "this is the file size"
        echo "$FILESIZE"
        echo "this is the size limiter"
        echo "$size"
        cat firstSplitGenome${i}.fa >> splitGenome_${countTotal}.fasta
        rm firstSplitGenome${i}.fa
    done
fi

