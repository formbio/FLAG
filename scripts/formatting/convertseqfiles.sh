#!/bin/bash
#convert_files.sh

usage() {
  echo "-h Help documentation for convertseqfiles.sh"
  echo "-i  --input file"
  echo "-p  --prefix for output file name"
  echo "-t  --input type"
  echo "-o  --outfile"
  echo "-a  --algo"
  echo "Example: bash create_a3m_msa.sh -r seq.fa -i input.fa -p out -a jackhammer -t 6"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:t:a:c:o:p:h opt
do
    case $opt in
        i) infile=$OPTARG;;
        o) outfile=$OPTARG;;
	c) contents=$OPTARG;;
        p) prefix=$OPTARG;;
	a) algo=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"

# Check for mandatory options

if [[ -f $algo ]]
then
    temp=$algo
    algo=$(head -n 1 $temp)
fi

if [[ -z $threads ]]
then
    threads=`nproc`
fi

if [[ -f $prefix ]]
then
    echo "Prefix is a file"
    cat $prefix
    temp=$(head -n 1 $prefix)
    prefix=$temp
elif [[ $prefix == '/run-data/inputs/seqname' ]] || [[ -z $prefix ]]
then
    if [[ -f $infile ]]
    then
	prefix="${infile%.*}"
    else
	prefix='sequence'
    fi
else
    echo "Prefix is $prefix"
fi

if [[ $algo == 'splitfa' ]]
then
    opt=''
    if [[ -n $prefix ]]
    then
	opt='--by-size-prefix $prefix'
    fi
    seqkit split2 --by-size 1 -O fafiles $infile
elif [[ $algo == 'concatfx' ]]
then
    fxs=$@
    if [[ -z $outfile ]]
    then
	outfile="${prefix}.${ext}"
    fi	
    seqkit concat $fxs --full | gzip > $outfile
elif [[ $algo == 'makefa' ]]
then
    if [[ -f $infile ]]
    then
        contents=$(grep -v ">" $infile)
    fi
    if [[ -z $outfile ]]
    then
	outfile="${prefix}.fa"
    fi	
    echo "makingseq"
    echo ">$prefix" > $outfile
    echo "$contents" >> $outfile
elif [[ $algo == 'fq2fa' ]] 
then
    echo "seqkit fq2fa $infile -o $outfile"
    seqkit fq2fa $infile -o $outfile
elif [[ $algo == 'tbl2fx' ]] 
then
    echo "seqkit tab2fx $infile > $outfile"
    seqkit tab2fx $infile > $outfile
elif [[ $algo == 'interleavefq' ]]
then
    fqs=$@
    if [[ -z $outfile ]]
    then
	outfile="$prefix.fastq.gz"
    fi
    echo "seqtk mergepe $fqs |gzip > $outfile"
    seqtk mergepe $fqs | gzip > $outfile
    
else
    echo "no algo"
fi
