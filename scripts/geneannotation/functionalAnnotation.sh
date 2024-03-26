#!/bin/bash
#entap.sh

usage() {
  echo "-h Help documentation for entap.sh"
  echo "-i  --input genome fasta"
  echo "-a  --input annotation"
  echo "-d  --databases"
  echo "-r  --reformat annotation (True or False)"
  echo "-p  --program (eggnog or entap)"
  echo "-t  --threads"
  echo "Example: bash evidencemodeler.sh -i genome.fa -l -t "
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:a:d:t:p:r:h opt
do
    case $opt in
        i) input=$OPTARG;;
        a) annotation=$OPTARG;;
        d) databases=$OPTARG;;
        r) reformat=$OPTARG;;
        p) program=$OPTARG;;
        t) threads=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${input} ]] || [[ -z ${annotation} ]]
then
    usage
fi
if [[ -s ${databases} ]]
then
    mkdir -p db
    tar xvfz ${databases} --strip-components=1 -C db
    mkdir /dbs
    mv db/* /dbs
else
    echo "Missing Database File"
    usage
fi
if [[ -z $threads ]]
then
    threads=`nproc`
fi
if [[ -z $program ]]
then
    program="eggnog"
fi

#checking if the fasta is a .gz
faEnd="${input: -3}"
if [[ "${faEnd}" == ".gz" ]]; then
    cp $input genome.fa.gz
    gunzip genome.fa.gz
else
    cp $input genome.fa
fi

agat_sp_extract_sequences.pl --clean_final_stop --gff $annotation -f genome.fa -p -o protein.fa

#EnTAP --config -d eggnog4.clustered_proteins.fa --out-dir /data/dbs/ --ini /opt/EnTAP/entap_config.ini

if [ "${program}" == "entap" ]; then

    EnTAP --runP -i protein.fa -d /dbs/uniprot_sprot.dmnd --ini /opt/EnTAP/entap_config.ini --threads $threads

    mv entap_outfiles/final_results/* .
else
    export EGGNOG_DATA_DIR=/dbs/

    emapper.py -i protein.fa -o eggnog --evalue 0.05 --cpu ${threads}
    cp eggnog.emapper.annotations eggnog_results.tsv
fi

#rm genome.fa

