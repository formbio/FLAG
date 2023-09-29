#!/bin/bash
#exonerate_p2g.sh

usage() {
  echo "-h Help documentation for tasser.sh"
  echo "-i  --input entap tsv"
  echo "-a  --input annotation gtf file"
  echo "-s  --species name"
  echo "-l  --busco lineage"
  echo "-g  --genome assembly file"
  echo "Example: bash exonerate_p2g.sh -p GCF_000001905.1_Loxafr3.0_protein.faa -g GCF_000001905.1_Loxafr3.0_genomic.fna -t 24"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:a:s:l:g:h opt
do
    case $opt in
        i) entaptsv=$OPTARG;;
        a) annotation=$OPTARG;;
        s) speciesName=$OPTARG;;
        l) lineage=$OPTARG;;
        g) genome=$OPTARG;;
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

############################################################
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
unset __conda_setup
# <<< conda initialize <<<
############################################################

baseDir="`dirname \"$0\"`"

cat final_annotations_lvl0.tsv | awk -F '\t' '{ print $1"\t"$25"\t"$28"\t"$29"\t"$32"\t"$2 }' >> simplifiedEntap.tsv

sed -i 's/ /_/g' simplifiedEntap.tsv

sed -i '/^#/d' "$annotation"

awk '/\tgene\t/ || /\ttRNA\t/{n++}{print >"annotSplit"  n ".gtf" }' "$annotation"

rm $annotation

for file in ./annotSplit*.gtf; do
    echo -e "bash ${baseDir}/renameAnnots.sh -a ${file} -s ${speciesName}" >> renameParallel.txt
done

parallel < renameParallel.txt > renameParallel.log

count=0
for file in ./updatedAnnot_*.gtf; do
    count=$(( $count + 1 ))
done

for (( i=1; i<=$count; i++ )); do
    cat updatedAnnot_${i}.gtf >> finalAnnotation.gtf
    rm updatedAnnot_${i}.gtf
done

echo "adding comments to annotation"
currentTime=`date`
sed -i "1 i\#!Species: ${speciesName}" finalAnnotation.gtf
sed -i "1 i\#!Form Bio annotation workflow version 2.0" finalAnnotation.gtf
sed -i "1 i\#!Created with the Form Bio annotation workflow on: ${currentTime}" finalAnnotation.gtf
sed -i "1 i\##gtf-version 3" finalAnnotation.gtf

mv finalAnnotation.gtf final${speciesName}.gtf

# make the protein fasta file
agat_sp_extract_sequences.pl --clean_final_stop --gff final${speciesName}.gtf -f genome.fasta -p -o proteins_${speciesName}.fa

# make the cdna fasta file
agat_sp_extract_sequences.pl --clean_final_stop --gff final${speciesName}.gtf -f genome.fasta --cdna -o cdna_${speciesName}.fa

# convert gtf to gff to supply both in the final output
agat_convert_sp_gxf2gxf.pl -g final${speciesName}.gtf -o final${speciesName}.gff3

# Input file
input_file="final${speciesName}.gff3"

# Temporary file for storing modified content
temp_file="temp_file.txt"

# Loop through each line in the input file
while IFS= read -r line; do
    if [[ $line == *"gene_biotype=protein_coding"* ]]; then
        # Replace "FLAG transcript" with "FLAG mrna"
        line=$(echo "$line" | sed 's/FormBio\ttranscript/FormBio\tmRNA/g')
    elif [[ $line == *"gene_biotype=tRNA"* ]]; then
        # Replace "FLAG transcript" with "FLAG mrna"
        line=$(echo "$line" | sed 's/FormBio\ttranscript/FormBio\ttRNA/g')
    fi
    # Append the modified line to the temporary file
    echo "$line" >> "$temp_file"
done < "$input_file"

mv "$temp_file" "final${speciesName}.gff3"

# Get agat stats
agat_sp_statistics.pl --gff final${speciesName}.gff3 --output final${speciesName}.AGAT.stats

# Get busco stats
conda activate BUSCO
busco -i proteins_${speciesName}.fa -l ${lineage} -o buscoout -m protein -c ${threads}
conda deactivate

mv buscoout/short_summary.*.buscoout.txt .



