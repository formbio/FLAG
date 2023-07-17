#!/bin/bash
#evidencemodeler.sh

usage() {
  echo "-h Help documentation for liftoff.sh"
  echo "-g  --input genome file in fasta format"
  echo "-a  --input EVM annoation"
  echo "-n  --name for output"
  echo "-f  --is the genome already formatted? (True or False)"
  echo "Example: bash evidencemodeler.sh -i genome.fa -l -t "
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :g:a:n:f:h opt
do
    case $opt in
        g) genome=$OPTARG;;
        a) annotation=$OPTARG;;
        n) name=$OPTARG;;
        f) formatted=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${genome} ]] || [[ -z ${annotation} ]]
then
    usage
fi

if [[ -z ${name} ]]
then
    name="gFacs"
fi
if [[ -z ${formatted} ]]
then
    formatted="False"
fi

if [[ "$formatted" == "False" ]]; then
    mkdir split
    mv $genome split/genome.fa

    cd split/

    awk '/>/{n++}{print >"___"  n ".fasta" }' "genome.fa"

    for file in *.fasta; do
        sed -i '1s/\s.*$//' $file
    done

    for file in *.fasta; do
        cat $file >> formattedgenome.fa
    done

    cp formattedgenome.fa ../genome.fa

    cd ..
fi

/opt/gFACs-master/gFACs.pl -f EVM_1.1.1_gff3 -p "$name" --statistics --statistics-at-every-step --rem-5prime-3prime-incompletes --rem-all-incompletes --min-exon-size 6 --min-intron-size 9 --min-CDS-size 50 --unique-genes-only --fasta genome.fa --splice-table --canonical-only --rem-genes-without-start-codon --allow-alternate-starts --rem-genes-without-stop-codon --rem-genes-without-start-and-stop-codon --allowed-inframe-stop-codons 0 --nt-content --get-fasta-with-introns --get-fasta-without-introns --get-protein-fasta --create-gtf --create-simple-gtf --create-gff3 --annotated-all-genes-only --distributions exon_lengths 10 ntron_lengths 15 CDS_lengths 20 gene_lengths 100 exon_position exon_position_data intron_position intron_position_data --compatibility SnpEff EVM_1.1.1_gene_prediction EVM_1.1.1_alignment -O . $annotation
