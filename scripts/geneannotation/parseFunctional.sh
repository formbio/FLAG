#!/bin/bash
#exonerate_p2g.sh

usage() {
  echo "-h Help documentation for tasser.sh"
  echo "-i  --input entap tsv"
  echo "-a  --input annotation gtf file"
  echo "-s  --species name"
  echo "-l  --busco lineage"
  echo "-g  --genome assembly file"
  echo "-p  --program (eggnog or entap)"
  echo "-d  --flag annotation checker database"
  echo "Example: bash exonerate_p2g.sh -p GCF_000001905.1_Loxafr3.0_protein.faa -g GCF_000001905.1_Loxafr3.0_genomic.fna -t 24"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:a:s:l:g:p:d:h opt
do
    case $opt in
        i) entaptsv=$OPTARG;;
        a) annotation=$OPTARG;;
        s) speciesName=$OPTARG;;
        l) lineage=$OPTARG;;
        g) genome=$OPTARG;;
        p) program=$OPTARG;;
        d) flagdb=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${entaptsv} ]] || [[ -z ${annotation} ]] || [[ -z ${genome} ]]
then
    usage
fi
if [[ -z ${speciesName} ]]
then
    speciesName=Sample_species
fi
if [[ -z ${threads} ]]
then
    threads=`nproc`
fi
if [[ -z ${program} ]]
then
    program="eggnog"
fi

#checking if the fasta is a .gz
faEnd="${genome: -3}"
if [[ "${faEnd}" == ".gz" ]]; then
    cp $genome genome.fa.gz
    gunzip genome.fa.gz
else
    cp $genome genome.fa
fi

############################################################
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
unset __conda_setup
# <<< conda initialize <<<
############################################################

baseDir="`dirname \"$0\"`"

if [[ "${program}" == "entap" ]]; then
    cat $entaptsv | awk -F '\t' '{ print $1"\t"$25"\t"$28"\t"$29"\t"$32"\t"$2 }' >> simplifiedFunctional.tsv
else
    cat $entaptsv  | awk -F '\t' '{ print $1"\t"$2"\t"$5"\t"$8"\t"$9 }' >> simplifiedFunctional.tsv
fi
sed -i 's/ /_/g' simplifiedFunctional.tsv

sed -i '/^#/d' "$annotation"

python3 /seqprg/scripts/geneannotation/renameAnnots.py --annotation_file $annotation --functional_simplified simplifiedFunctional.tsv --output_file finalAnnotation.gtf

#rm the intermediate file
rm CorrectedAnnotation.gtf

agat_convert_sp_gxf2gxf.pl -g finalAnnotation.gtf -o finalAnnotation.gff3

sed -i '/AGAT\tfive_prime_UTR/d' finalAnnotation.gff3
sed -i '/AGAT\tthree_prime_UTR/d' finalAnnotation.gff3

rm finalAnnotation.gtf

agat_convert_sp_gxf2gxf.pl --gff finalAnnotation.gff3 -o finalAnnotationSorted.gff3

mv finalAnnotationSorted.gff3 finalAnnotation.gff3

agat_convert_sp_gff2gtf.pl --gff finalAnnotation.gff3 -o finalAnnotation.gtf

echo "adding comments to annotation"
currentTime=`date`
sed -i "1 i\#!Species: ${speciesName}" finalAnnotation.gtf
sed -i "1 i\#!Form Bio annotation workflow version 2.0" finalAnnotation.gtf
sed -i "1 i\#!Created with the Form Bio annotation workflow on: ${currentTime}" finalAnnotation.gtf
sed -i "1 i\##gtf-version 3" finalAnnotation.gtf

mv finalAnnotation.gtf final${speciesName}.gtf
mv finalAnnotation.gff3 final${speciesName}.gff3
sed -i 's/""""/""/g' final${speciesName}.gtf
sed -i 's/""""/""/g' final${speciesName}.gff3

sed -i '/\tfive_prime_UTR\t/d' final${speciesName}.gtf
sed -i '/\tthree_prime_UTR\t/d' final${speciesName}.gtf

sed -i '/\tfive_prime_UTR\t/d' final${speciesName}.gff3
sed -i '/\tthree_prime_UTR\t/d' final${speciesName}.gff3

# make the protein fasta file
#echo "agat_sp_extract_sequences.pl --clean_final_stop --gff final${speciesName}.gff3 -f genome.fa -p -o proteins_${speciesName}.fa" >> parallel_agat_commands.txt
agat_sp_extract_sequences.pl --clean_final_stop --gff final${speciesName}.gff3 -f genome.fa -p -o proteins_${speciesName}.fa
# make the cdna fasta file
#echo "agat_sp_extract_sequences.pl --clean_final_stop --gff final${speciesName}.gff3 -f genome.fa --cdna -o cdna_${speciesName}.fa" >> parallel_agat_commands.txt
agat_sp_extract_sequences.pl --clean_final_stop --gff final${speciesName}.gff3 -f genome.fa --cdna -o cdna_${speciesName}.fa

# convert gtf to gff to supply both in the final output
#agat_convert_sp_gxf2gxf.pl -g final${speciesName}.gtf -o final${speciesName}.gff3

# Get agat stats
#echo "agat_sp_statistics.pl --gff final${speciesName}.gff3 --output final${speciesName}.AGAT.stats" >> parallel_agat_commands.txt
agat_sp_statistics.pl --gff final${speciesName}.gff3 --output final${speciesName}.AGAT.stats
#parallel < parallel_agat_commands.txt

# Get busco stats
conda activate BUSCO
busco -i proteins_${speciesName}.fa -l ${lineage} -o buscoout -m protein -c ${threads}
conda deactivate

mv buscoout/short_summary.*.buscoout.txt .
cp short_summary.*.buscoout.txt busco.txt

python3 /seqprg/scripts/geneannotation/finalReport.py --input_stats "final${speciesName}.AGAT.stats" --lineage "${lineage}" --scientific_name "${speciesName}"

conda activate biopython
python /seqprg/scripts/formatting/gff_to_genbank.py "final${speciesName}.gff3" "genome.fa"
conda deactivate


