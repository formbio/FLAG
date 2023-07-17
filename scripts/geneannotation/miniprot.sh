#!/bin/bash
#exonerate_p2g.sh

usage() {
  echo "-h Help documentation for tasser.sh"
  echo "-n  --name"
  echo "-p  --input protein file"
  echo "-g  --input genome file"
  echo "-t  --threads"
  echo "-o  --output type (gtf or paf)"
  echo "Example: bash exonerate_p2g.sh -p GCF_000001905.1_Loxafr3.0_protein.faa -g GCF_000001905.1_Loxafr3.0_genomic.fna -t 24"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :n:p:g:t:o:h opt
do
    case $opt in
        n) name=$OPTARG;;
        p) protein=$OPTARG;;
        g) genome=$OPTARG;;
        t) threads=$OPTARG;;
        o) outType=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${protein} ]] || [[ -z ${genome} ]]
then
    usage
fi
if [[ -z ${threads} ]]; then
    threads=90
fi
if [[ -z ${outType} ]]; then
    outType="gtf"
fi
if [[ -z ${name} ]]; then
    name="none"
fi

mv $protein protein.fa
genomeEnd="${genome: -7}"
echo "end of genome file is: ${genomeEnd}. This is to check if its a .tar.gz file"
if [[ "${genomeEnd}" != ".tar.gz" ]]; then
    mv $genome genome.fa
else
    mkdir -p db
    tar xvfz ${genome} --strip-components=1 -C db
    dbname=$(ls db/*.fa)
    if [[ "$dbname" == "" ]]; then
        dbname=$(ls db/*.fasta)
        if [[ "$dbname" == "" ]]; then
            echo "no fasta sequence for database found"
            usage
        fi
    fi
    mv $dbname genome.fa
fi
genomeName="${genome%%.*}"
#numchunks=$(grep -o '>' protein.fa | wc -l)
    # here the protein is the query and the genome is the target. The protein is being mapped to the genome
    #exonerate --model protein2genome $protein $genome --showtargetgff TRUE --softmasktarget TRUE > exonerate.out
    
#in miniprot the genome or nt sequence is the database
if [[ "$outType" == "gtf" ]]; then
    miniprot -ut${threads} --gtf genome.fa protein.fa > ${name}_miniprot.gtf
    if [ "$name" == "formatted_referenceProtein" ]
    then
        mv ${name}_miniprot.gtf reference_miniprot.gtf
    elif [ "$name" == "none" ]; then
        mv ${name}_miniprot.gtf miniprot.gtf
    fi
else
    miniprot -ut${threads} genome.fa protein.fa > ${name}_miniprot.paf
    cp ${name}_miniprot.paf ${name}_miniprot.txt
    sed -i "1s/^/QSeqName\tQSeqLength\tQSstartCoordinate\tQEndCoordinate\tQTStrandMatch\tTSeqName\tTSeqLength\tTSstartCoordinate\tTEndCoordinate\tNumMatchingBases\tNumBases\tMappingQuality\n/" ${name}_miniprot.txt
    sed -i "1s/^/#miniprot alignment of ${name} to ${genomeName}\n/" ${name}_miniprot.txt
fi


