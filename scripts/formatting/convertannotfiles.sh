#!/bin/bash
#convert_files.sh

usage() {
  echo "-h Help documentation for convert_files.sh"
  echo "-i  --input file"
  echo "-p  --prefix for output file name"
  echo "-t  --input type"
  echo "-o  --outfile"
  echo "Example: bash convertannotfiles.sh -i infile -t intype -p outfile/seqname-prefix -o outfile"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:t:o:p:h opt
do
    case $opt in
        i) infile=$OPTARG;;
        t) intype=$OPTARG;;
        o) outfile=$OPTARG;;
        p) prefix=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $intype ]]
then
    usage
fi

if [[ -z $threads ]]
then
    threads=`nproc`
fi
if [[ -f $intype ]]
then
    algo=$(head -n 1 $intype)
else
    algo=$intype
fi

if [[ -f $prefix ]]
then
    pfile=$prefix
    prefix=$(head -n 1 $pfile)
elif [[ -z $prefix ]]
then
    prefix="${infile%.*}"
fi

baseDir="`dirname \"$0\"`"
echo $algo
if [[ $algo == 'snapgene' ]]
then
    python3 $baseDir/snpgene2fasta.py -d $infile -o $outfile -n $prefix
    bed="${outfile}.bed"
    fasta="${outfile}.fa"
    outfile=$bed
elif [[ $algo == 'gff' ]]
then
    if [[ -z $outfile ]]
    then
	outfile="${prefix}.gtf"
    fi
    agat_convert_sp_gff2gtf.pl --gff $infile -o $outfile
elif [[ $algo == 'gtf' ]] 
then
    if [[ -z $outfile ]]
    then
	outfile="${prefix}.gff"
    fi
    agat_convert_sp_gxf2gxf.pl -g $infile -o $outfile
fi

echo $outfile
cat $outfile
