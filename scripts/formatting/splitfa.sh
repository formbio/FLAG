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
while getopts :i:n:s:o:c:h opt
do
    case $opt in
        i) fa=$OPTARG;;
        n) num=$OPTARG;;
	o) outseqname=$OPTARG;;
        c) check=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"

# Check for mandatory options
if [[ -z $fa ]]
then
    usage
fi
if [[ -z $num ]]
then
    num=$(grep -c ">" $fa)
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
opt=''
if [[ $check == 'NT' ]]
then
    opt='-c NT'
fi
if [[ -n $outseqname ]]
then
    opt="$opt -o 1"
fi
perl $baseDir/format_fasta.pl -f $fa -n $num $opt
