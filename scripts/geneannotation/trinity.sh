#!/bin/bash
#trinity.sh

usage() {
  echo "-h Help documentation for tasser.sh"
  echo "-f  --input fasta sequence 1"
  echo "-r  --input fasta sequence 2 if its a paired read"
  echo "-t  --number of cpus"
  echo "-s  --seqtype fasta or fastaq"
  echo "-n  --name type"
  echo "-m  --max RAM"
  echo "-i  --max intron"
  echo "Example: bash trinity.sh -f LIB034841_TRA00126514_S11_L001_R1.fastq -r LIB034841_TRA00126514_S11_L001_R2.fastq -t 30 -m 120G"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :f:r:t:m:n:s:i:h opt
do
    case $opt in
        f) f1=$OPTARG;;
        r) f2=$OPTARG;;
        t) cpu=$OPTARG;;
        m) mem=$OPTARG;;
        n) nameType=$OPTARG;;
        s) seqtype=$OPTARG;;
        i) intron=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $f1 ]]
then
    usage
fi
if [[ -z ${cpu} ]]
then
    cpu=16
fi
if [[ -z ${mem} ]]
then
    mem=32G
fi
if [[ -z ${intron} ]]
then
    intron=10000
fi

if [ "$seqtype" == "bam" ]; then
    Trinity --genome_guided_bam $f1 --CPU $cpu --max_memory $mem --genome_guided_max_intron $intron --bflyHeapSpaceMax "80G"
    if [ "$f1" == "homology.star.bam" ]; then
        mv trinity_out_dir.Trinity.fasta reference_trinity.fasta
    else
        mv trinity_out_dir.Trinity.fasta trinity.fasta
    fi
else
    if [[ -z $f2 ]]
    then
        #perl /opt/trinityrnaseq-v2.13.2/util/insilico_read_normalization.pl --seqType fa --JM $mem --max_cov 30 --CPU $cpu --single $f1
        #rm $f1
        #Trinity --seqType $seqtype --single single.norm.fa --CPU $cpu --max_memory $mem --bflyHeapSpaceMax "80G"
        Trinity --seqType $seqtype --single $f1 --CPU $cpu --max_memory $mem --bflyHeapSpaceMax "180G"
        if [[ "$nameType" == "rename" ]]; then
            prefix="${f1%.fa*}"
            mv trinity_out_dir.Trinity.fasta ${prefix}_trinity.fasta
        elif [ "$f1" == "formatted_referenceRNA.fa" ]
        then
            mv trinity_out_dir.Trinity.fasta reference_trinity.fasta
        else
            mv trinity_out_dir.Trinity.fasta trinity.fasta
        fi
    else
        Trinity --seqType $seqtype --left $f1 --right $f2 --CPU $cpu --max_memory $mem --bflyHeapSpaceMax "80G"
    fi
fi
