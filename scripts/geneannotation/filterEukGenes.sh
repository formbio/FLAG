#!/bin/bash
#evidencemodeler.sh

usage() {
  echo "-h Help documentation for liftoff.sh"
  echo "-i  --input genome file in fasta format"
  echo "-l  --lineage"
  echo "-s  --size"
  echo "-z  --threads"
  echo "Example: bash evidencemodeler.sh -i genome.fa"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:l:s:z:h opt
do
    case $opt in
        i) genome=$OPTARG;;
        l) lineage=$OPTARG;;
        s) size=$OPTARG;;
        z) threads=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${genome} ]]
then
    usage
fi
if [[ -z ${lineage} ]]
then
    buscorun="False"
fi
if [[ -z $threads ]]
then
    threads=`nproc`
fi
if [[ -z ${size} ]]
then
    size="normal"
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

#process the inputs and format them for EVM and make the weights.txt file
mv ${genome} genome.fasta
touch weights.txt
touch gene_predictions.gff3
touch protein_alignments.gff3
touch transcript_alignments.gff3
touch other_predictions.gff3

samtools faidx genome.fasta

touch completeBuscos.gff3
augAdded="False"

if [ -f "reference_augustus.gtf" ]; then
    touch refAugCommands.txt
    echo -e "############################################################" > refAugCommands.txt
    echo -e "# >>> conda initialize >>>" >> refAugCommands.txt
    echo -e "# !! Contents within this block are managed by 'conda init' !!" >> refAugCommands.txt
    echo -e "__conda_setup="\""\$('conda' 'shell.bash' 'hook' 2> /dev/null)"\""" >> refAugCommands.txt
    echo -e 'eval "$__conda_setup"' >> refAugCommands.txt
    echo -e "unset __conda_setup" >> refAugCommands.txt
    echo -e "# <<< conda initialize <<<" >> refAugCommands.txt
    echo -e "############################################################" >> refAugCommands.txt
    
    
    echo -e "mkdir refAug" >> refAugCommands.txt
    echo -e "cp reference_augustus.gtf refAug/" >> refAugCommands.txt
    echo -e "cd refAug/" >> refAugCommands.txt
    echo -e "perl /opt/EVidenceModeler/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl reference_augustus.gtf > reference_augustus.gff3" >> refAugCommands.txt
    echo -e "gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o tmp.gff3  reference_augustus.gff3" >> refAugCommands.txt
    echo -e "mv tmp.gff3 reference_augustus.gff3" >> refAugCommands.txt
    echo -e "cp ../genome.fasta ." >> refAugCommands.txt
    echo -e "agat_sp_extract_sequences.pl --clean_final_stop --gff reference_augustus.gff3 -f genome.fasta -p -o protein.fa" >> refAugCommands.txt
    echo -e "conda activate BUSCO" >> refAugCommands.txt
    echo -e "busco -i protein.fa -l ${lineage} -o buscoout -m protein -c 8" >> refAugCommands.txt
    echo -e "mv buscoout/run_${lineage}/full_table.tsv ." >> refAugCommands.txt
    echo -e "conda deactivate" >> refAugCommands.txt
    echo -e "bash ${baseDir}/annotationFilter.sh -i reference_augustus.gff3 -o True -p Augustusgff3" >> refAugCommands.txt
    echo -e "mv AugBuscos.gff3 complete_Refaugustus.gff3" >> refAugCommands.txt
    echo -e "mv AugustusRemoved.gff3 removed_Refaugustus.gff3" >> refAugCommands.txt
    echo -e "sed -i -E 's/([[:digit:]])-g/\1-gR/g' complete_Refaugustus.gff3" >> refAugCommands.txt
    echo -e "sed -i -E 's/([[:digit:]])-g/\1-gR/g' removed_Refaugustus.gff3" >> refAugCommands.txt
    echo -e "sed -i -E 's/([[:digit:]])-g/\1-gR/g' reference_augustus.gff3" >> refAugCommands.txt
    echo -e "mv complete_Refaugustus.gff3 .." >> refAugCommands.txt
    echo -e "mv removed_Refaugustus.gff3 .." >> refAugCommands.txt
    echo -e "mv reference_augustus.gff3 .." >> refAugCommands.txt
    echo -e "cd .." >> refAugCommands.txt
    echo -e "cat complete_Refaugustus.gff3 >> completeBuscos.gff3" >> refAugCommands.txt
    echo -e "cat reference_augustus.gff3 >> gene_predictions.gff3" >> refAugCommands.txt

    echo -e "bash refAugCommands.txt" >> filterCommands.txt
    
    echo -e "ABINITIO_PREDICTION\tAugustus\t2" >> weights.txt
    augAdded="True"
fi

if [ -f "pretrained_augustus.gtf" ]; then
    touch preAugCommands.txt
    echo -e "############################################################" > preAugCommands.txt
    echo -e "# >>> conda initialize >>>" >> preAugCommands.txt
    echo -e "# !! Contents within this block are managed by 'conda init' !!" >> preAugCommands.txt
    echo -e "__conda_setup="\""\$('conda' 'shell.bash' 'hook' 2> /dev/null)"\""" >> preAugCommands.txt
    echo -e 'eval "$__conda_setup"' >> preAugCommands.txt
    echo -e "unset __conda_setup" >> preAugCommands.txt
    echo -e "# <<< conda initialize <<<" >> preAugCommands.txt
    echo -e "############################################################" >> preAugCommands.txt
    
    echo -e "mkdir preAug" >> preAugCommands.txt
    echo -e "cp pretrained_augustus.gtf preAug/" >> preAugCommands.txt
    echo -e "cd preAug/" >> preAugCommands.txt
    echo -e "perl /opt/EVidenceModeler/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl pretrained_augustus.gtf > pretrained_augustus.gff3" >> preAugCommands.txt
    echo -e "gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o tmp.gff3  pretrained_augustus.gff3" >> preAugCommands.txt
    echo -e "mv tmp.gff3 pretrained_augustus.gff3" >> preAugCommands.txt
    echo -e "cp ../genome.fasta ." >> preAugCommands.txt
    echo -e "agat_sp_extract_sequences.pl --clean_final_stop --gff pretrained_augustus.gff3 -f genome.fasta -p -o protein.fa" >> preAugCommands.txt
    echo -e "conda activate BUSCO" >> preAugCommands.txt
    echo -e "busco -i protein.fa -l ${lineage} -o buscoout -m protein -c 8" >> preAugCommands.txt
    echo -e "mv buscoout/run_${lineage}/full_table.tsv ." >> preAugCommands.txt
    echo -e "conda deactivate" >> preAugCommands.txt
    echo -e "bash ${baseDir}/annotationFilter.sh -i pretrained_augustus.gff3 -o True -p Augustusgff3" >> preAugCommands.txt
    echo -e "mv AugBuscos.gff3 complete_Preaugustus.gff3" >> preAugCommands.txt
    echo -e "mv AugustusRemoved.gff3 removed_Preaugustus.gff3" >> preAugCommands.txt
    echo -e "sed -i -E 's/([[:digit:]])-g/\1-gP/g' complete_Preaugustus.gff3" >> preAugCommands.txt
    echo -e "sed -i -E 's/([[:digit:]])-g/\1-gP/g' removed_Preaugustus.gff3" >> preAugCommands.txt
    echo -e "sed -i -E 's/([[:digit:]])-g/\1-gP/g' pretrained_augustus.gff3" >> preAugCommands.txt
    echo -e "mv complete_Preaugustus.gff3 .." >> preAugCommands.txt
    echo -e "mv removed_Preaugustus.gff3 .." >> preAugCommands.txt
    echo -e "mv pretrained_augustus.gff3 .." >> preAugCommands.txt
    echo -e "cd .." >> preAugCommands.txt
    echo -e "cat complete_Preaugustus.gff3 >> completeBuscos.gff3" >> preAugCommands.txt
    echo -e "cat pretrained_augustus.gff3 >> gene_predictions.gff3" >> preAugCommands.txt

    echo -e "bash preAugCommands.txt" >> filterCommands.txt
    
    if [[ "$augAdded" == "False" ]]; then
        echo -e "ABINITIO_PREDICTION\tAugustus\t2" >> weights.txt
        augAdded="True"
    fi
fi

if [ -f "liftoff_augustus.gtf" ]; then
    touch liftAugCommands.txt
    echo -e "############################################################" > liftAugCommands.txt
    echo -e "# >>> conda initialize >>>" >> liftAugCommands.txt
    echo -e "# !! Contents within this block are managed by 'conda init' !!" >> liftAugCommands.txt
    echo -e "__conda_setup="\""\$('conda' 'shell.bash' 'hook' 2> /dev/null)"\""" >> liftAugCommands.txt
    echo -e 'eval "$__conda_setup"' >> liftAugCommands.txt
    echo -e "unset __conda_setup" >> liftAugCommands.txt
    echo -e "# <<< conda initialize <<<" >> liftAugCommands.txt
    echo -e "############################################################" >> liftAugCommands.txt
    
    echo -e "mkdir liftAug" >> liftAugCommands.txt
    echo -e "cp liftoff_augustus.gtf liftAug/" >> liftAugCommands.txt
    echo -e "cd liftAug/" >> liftAugCommands.txt
    echo -e "perl /opt/EVidenceModeler/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl liftoff_augustus.gtf > liftoff_augustus.gff3" >> liftAugCommands.txt
    echo -e "gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o tmp.gff3  liftoff_augustus.gff3" >> liftAugCommands.txt
    echo -e "mv tmp.gff3 liftoff_augustus.gff3" >> liftAugCommands.txt
    echo -e "cp ../genome.fasta ." >> liftAugCommands.txt
    echo -e "agat_sp_extract_sequences.pl --clean_final_stop --gff liftoff_augustus.gff3 -f genome.fasta -p -o protein.fa" >> liftAugCommands.txt
    echo -e "conda activate BUSCO" >> liftAugCommands.txt
    echo -e "busco -i protein.fa -l ${lineage} -o buscoout -m protein -c 8" >> liftAugCommands.txt
    echo -e "mv buscoout/run_${lineage}/full_table.tsv ." >> liftAugCommands.txt
    echo -e "conda deactivate" >> liftAugCommands.txt
    echo -e "bash ${baseDir}/annotationFilter.sh -i liftoff_augustus.gff3 -o True -p Augustusgff3" >> liftAugCommands.txt
    echo -e "mv AugBuscos.gff3 complete_liftaugustus.gff3" >> liftAugCommands.txt
    echo -e "mv AugustusRemoved.gff3 removed_liftaugustus.gff3" >> liftAugCommands.txt
    echo -e "sed -i -E 's/([[:digit:]])-g/\1-gL/g' complete_liftaugustus.gff3" >> liftAugCommands.txt
    echo -e "sed -i -E 's/([[:digit:]])-g/\1-gL/g' removed_liftaugustus.gff3" >> liftAugCommands.txt
    echo -e "sed -i -E 's/([[:digit:]])-g/\1-gL/g' liftoff_augustus.gff3" >> liftAugCommands.txt
    echo -e "mv complete_liftaugustus.gff3 .." >> liftAugCommands.txt
    echo -e "mv removed_liftaugustus.gff3 .." >> liftAugCommands.txt
    echo -e "mv liftoff_augustus.gff3 .." >> liftAugCommands.txt
    echo -e "cd .." >> liftAugCommands.txt
    echo -e "cat complete_liftaugustus.gff3 >> completeBuscos.gff3" >> liftAugCommands.txt
    echo -e "cat liftoff_augustus.gff3 >> gene_predictions.gff3" >> liftAugCommands.txt

    echo -e "bash liftAugCommands.txt" >> filterCommands.txt
    
    if [[ "$augAdded" == "False" ]]; then
        echo -e "ABINITIO_PREDICTION\tAugustus\t2" >> weights.txt
        augAdded="True"
    fi
fi

if [ -f "Helixer_augustus.gtf" ]; then
    touch helixAugCommands.txt
    echo -e "############################################################" > helixAugCommands.txt
    echo -e "# >>> conda initialize >>>" >> helixAugCommands.txt
    echo -e "# !! Contents within this block are managed by 'conda init' !!" >> helixAugCommands.txt
    echo -e "__conda_setup="\""\$('conda' 'shell.bash' 'hook' 2> /dev/null)"\""" >> helixAugCommands.txt
    echo -e 'eval "$__conda_setup"' >> helixAugCommands.txt
    echo -e "unset __conda_setup" >> helixAugCommands.txt
    echo -e "# <<< conda initialize <<<" >> helixAugCommands.txt
    echo -e "############################################################" >> helixAugCommands.txt
    
    echo -e "mkdir helixAug" >> helixAugCommands.txt
    echo -e "cp Helixer_augustus.gtf helixAug/" >> helixAugCommands.txt
    echo -e "cd helixAug/" >> helixAugCommands.txt
    echo -e "perl /opt/EVidenceModeler/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl Helixer_augustus.gtf > Helixer_augustus.gff3" >> helixAugCommands.txt
    echo -e "gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o tmp.gff3  Helixer_augustus.gff3" >> helixAugCommands.txt
    echo -e "mv tmp.gff3 Helixer_augustus.gff3" >> helixAugCommands.txt
    echo -e "cp ../genome.fasta ." >> helixAugCommands.txt
    echo -e "agat_sp_extract_sequences.pl --clean_final_stop --gff Helixer_augustus.gff3 -f genome.fasta -p -o protein.fa" >> helixAugCommands.txt
    echo -e "conda activate BUSCO" >> helixAugCommands.txt
    echo -e "busco -i protein.fa -l ${lineage} -o buscoout -m protein -c 8" >> helixAugCommands.txt
    echo -e "mv buscoout/run_${lineage}/full_table.tsv ." >> helixAugCommands.txt
    echo -e "conda deactivate" >> helixAugCommands.txt
    echo -e "bash ${baseDir}/annotationFilter.sh -i Helixer_augustus.gff3 -o True -p Augustusgff3" >> helixAugCommands.txt
    echo -e "mv AugBuscos.gff3 complete_Helixaugustus.gff3" >> helixAugCommands.txt
    echo -e "mv AugustusRemoved.gff3 removed_Helixaugustus.gff3" >> helixAugCommands.txt
    echo -e "sed -i -E 's/([[:digit:]])-g/\1-gH/g' complete_Helixaugustus.gff3" >> helixAugCommands.txt
    echo -e "sed -i -E 's/([[:digit:]])-g/\1-gH/g' removed_Helixaugustus.gff3" >> helixAugCommands.txt
    echo -e "sed -i -E 's/([[:digit:]])-g/\1-gH/g' Helixer_augustus.gff3" >> helixAugCommands.txt
    echo -e "mv complete_Helixaugustus.gff3 .." >> helixAugCommands.txt
    echo -e "mv removed_Helixaugustus.gff3 .." >> helixAugCommands.txt
    echo -e "mv Helixer_augustus.gff3 .." >> helixAugCommands.txt
    echo -e "cd .." >> helixAugCommands.txt
    echo -e "cat complete_Helixaugustus.gff3 >> completeBuscos.gff3" >> helixAugCommands.txt
    echo -e "cat Helixer_augustus.gff3 >> gene_predictions.gff3" >> helixAugCommands.txt

    echo -e "bash helixAugCommands.txt" >> filterCommands.txt
    
    if [[ "$augAdded" == "False" ]]; then
        echo -e "ABINITIO_PREDICTION\tAugustus\t2" >> weights.txt
        augAdded="True"
    fi
fi

if [ -f "augustus.gtf" ]; then
    touch fullAugCommands.txt
    echo -e "############################################################" > fullAugCommands.txt
    echo -e "# >>> conda initialize >>>" >> fullAugCommands.txt
    echo -e "# !! Contents within this block are managed by 'conda init' !!" >> fullAugCommands.txt
    echo -e "__conda_setup="\""\$('conda' 'shell.bash' 'hook' 2> /dev/null)"\""" >> fullAugCommands.txt
    echo -e 'eval "$__conda_setup"' >> fullAugCommands.txt
    echo -e "unset __conda_setup" >> fullAugCommands.txt
    echo -e "# <<< conda initialize <<<" >> fullAugCommands.txt
    echo -e "############################################################" >> fullAugCommands.txt
    
    
    echo -e "mkdir Aug" >> fullAugCommands.txt
    echo -e "cp augustus.gtf Aug/" >> fullAugCommands.txt
    echo -e "cd Aug/" >> fullAugCommands.txt
    echo -e "perl /opt/EVidenceModeler/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl augustus.gtf > augustus.gff3" >> fullAugCommands.txt
    echo -e "gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o tmp.gff3  augustus.gff3" >> fullAugCommands.txt
    echo -e "mv tmp.gff3 augustus.gff3" >> fullAugCommands.txt
    echo -e "cp ../genome.fasta ." >> fullAugCommands.txt
    echo -e "agat_sp_extract_sequences.pl --clean_final_stop --gff augustus.gff3 -f genome.fasta -p -o protein.fa" >> fullAugCommands.txt
    echo -e "conda activate BUSCO" >> fullAugCommands.txt
    echo -e "busco -i protein.fa -l ${lineage} -o buscoout -m protein -c 8" >> fullAugCommands.txt
    echo -e "mv buscoout/run_${lineage}/full_table.tsv ." >> fullAugCommands.txt
    echo -e "conda deactivate" >> fullAugCommands.txt
    echo -e "bash ${baseDir}/annotationFilter.sh -i augustus.gff3 -o True -p Augustusgff3" >> fullAugCommands.txt
    echo -e "mv AugustusRemoved.gff3 fullAugustusRemoved.gff3" >> fullAugCommands.txt
    echo -e "mv AugBuscos.gff3 fullAugBuscos.gff3" >> fullAugCommands.txt
    echo -e "mv fullAugBuscos.gff3 .." >> fullAugCommands.txt
    echo -e "mv fullAugustusRemoved.gff3 .." >> fullAugCommands.txt
    echo -e "mv augustus.gff3 .." >> fullAugCommands.txt

    echo -e "cd .." >> fullAugCommands.txt
    echo -e "cat fullAugBuscos.gff3 >> completeBuscos.gff3" >> fullAugCommands.txt
    echo -e "cat augustus.gff3 >> gene_predictions.gff3" >> fullAugCommands.txt

    echo -e "bash fullAugCommands.txt" >> filterCommands.txt
    
    if [[ "$augAdded" == "False" ]]; then
        echo -e "ABINITIO_PREDICTION\tAugustus\t2" >> weights.txt
        augAdded="True"
    fi
fi

if [ -f "Helixer.gff3" ]; then
    touch helixerCommands.txt
    echo -e "############################################################" > helixerCommands.txt
    echo -e "# >>> conda initialize >>>" >> helixerCommands.txt
    echo -e "# !! Contents within this block are managed by 'conda init' !!" >> helixerCommands.txt
    echo -e "__conda_setup="\""\$('conda' 'shell.bash' 'hook' 2> /dev/null)"\""" >> helixerCommands.txt
    echo -e 'eval "$__conda_setup"' >> helixerCommands.txt
    echo -e "unset __conda_setup" >> helixerCommands.txt
    echo -e "# <<< conda initialize <<<" >> helixerCommands.txt
    echo -e "############################################################" >> helixerCommands.txt
    
    echo -e "mkdir helixer" >> helixerCommands.txt
    echo -e "cp Helixer.gff3 helixer/" >> helixerCommands.txt
    echo -e "cd helixer/" >> helixerCommands.txt
    echo -e "gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o tmp.gff3  Helixer.gff3" >> helixerCommands.txt
    echo -e "mv tmp.gff3 Helixer.gff3" >> helixerCommands.txt
    echo -e "cp ../genome.fasta ." >> helixerCommands.txt
    echo -e "agat_sp_extract_sequences.pl --clean_final_stop --gff Helixer.gff3 -f genome.fasta -p -o protein.fa" >> helixerCommands.txt
    echo -e "conda activate BUSCO" >> helixerCommands.txt
    echo -e "busco -i protein.fa -l ${lineage} -o buscoout -m protein -c 8" >> helixerCommands.txt
    echo -e "mv buscoout/run_${lineage}/full_table.tsv ." >> helixerCommands.txt
    echo -e "conda deactivate" >> helixerCommands.txt
    echo -e "bash ${baseDir}/annotationFilter.sh -i Helixer.gff3 -o True -p Helixer" >> helixerCommands.txt
    echo -e "mv HelixerBuscos.gff3 .." >> helixerCommands.txt
    echo -e "mv HelixerRemoved.gff3 .." >> helixerCommands.txt
    echo -e "cd .." >> helixerCommands.txt
    ##echo -e "cat HelixerRemoved.gff3 >> gene_predictions.gff3" >> helixerCommands.txt
    echo -e "cat HelixerBuscos.gff3 >> completeBuscos.gff3" >> helixerCommands.txt
    echo -e "cat Helixer.gff3 >> gene_predictions.gff3" >> helixerCommands.txt
    
    echo -e "bash helixerCommands.txt" >> filterCommands.txt
    
    echo -e "ABINITIO_PREDICTION\tHelixer\t3" >> weights.txt
fi

if [ -f "liftoff.gff3" ]; then
    touch liftoffCommands.txt
    echo -e "############################################################" > liftoffCommands.txt
    echo -e "# >>> conda initialize >>>" >> liftoffCommands.txt
    echo -e "# !! Contents within this block are managed by 'conda init' !!" >> liftoffCommands.txt
    echo -e "__conda_setup="\""\$('conda' 'shell.bash' 'hook' 2> /dev/null)"\""" >> liftoffCommands.txt
    echo -e 'eval "$__conda_setup"' >> liftoffCommands.txt
    echo -e "unset __conda_setup" >> liftoffCommands.txt
    echo -e "# <<< conda initialize <<<" >> liftoffCommands.txt
    echo -e "############################################################" >> liftoffCommands.txt
    
    echo -e "mkdir liftoff" >> liftoffCommands.txt
    echo -e "cp liftoff.gff3 liftoff/" >> liftoffCommands.txt
    echo -e "cd liftoff/" >> liftoffCommands.txt
    echo -e "gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o tmp.gff3  liftoff.gff3" >> liftoffCommands.txt
    echo -e "mv tmp.gff3 liftoff.gff3" >> liftoffCommands.txt
    echo -e "cp ../genome.fasta ." >> liftoffCommands.txt
    echo -e "agat_sp_extract_sequences.pl --clean_final_stop --gff liftoff.gff3 -f genome.fasta -p -o protein.fa" >> liftoffCommands.txt
    echo -e "conda activate BUSCO" >> liftoffCommands.txt
    echo -e "busco -i protein.fa -l ${lineage} -o buscoout -m protein -c 8" >> liftoffCommands.txt
    echo -e "mv buscoout/run_${lineage}/full_table.tsv ." >> liftoffCommands.txt
    echo -e "conda deactivate" >> liftoffCommands.txt
    echo -e "bash ${baseDir}/annotationFilter.sh -i liftoff.gff3 -o True -p Liftoff" >> liftoffCommands.txt
    echo -e "gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o tmp.gff3  LiftoffBuscos.gff3" >> liftoffCommands.txt
    echo -e "mv tmp.gff3 LiftoffBuscos.gff3" >> liftoffCommands.txt
    echo -e "sed -i '/^#/d' LiftoffBuscos.gff3" >> liftoffCommands.txt
    echo -e "gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o tmp.gff3  LiftoffRemoved.gff3" >> liftoffCommands.txt
    echo -e "mv tmp.gff3 LiftoffRemoved.gff3" >> liftoffCommands.txt
    echo -e "mv LiftoffBuscos.gff3 .." >> liftoffCommands.txt
    echo -e "mv LiftoffRemoved.gff3 .." >> liftoffCommands.txt
    echo -e "cd .." >> liftoffCommands.txt
    ##echo -e "cat LiftoffRemoved.gff3 >> gene_predictions.gff3" >> liftoffCommands.txt
    echo -e "cat LiftoffBuscos.gff3 >> completeBuscos.gff3" >> liftoffCommands.txt
    echo -e "cat liftoff.gff3 >> gene_predictions.gff3" >> liftoffCommands.txt
    
    echo -e "bash liftoffCommands.txt" >> filterCommands.txt
    
    echo -e "OTHER_PREDICTION\tLiftoff\t3" >> weights.txt
fi

#run filtering
parallel < filterCommands.txt > filter.log


if [ -f "exonerate.out" ]; then
    perl /opt/EVidenceModeler/EvmUtils/misc/Exonerate_to_evm_gff3.pl exonerate.out > exonerate.gff3
    cat exonerate.gff3 >> protein_alignments.gff3
    echo -e "PROTEIN\texonerate\t4" >> weights.txt
fi

if [ -f "reference_exonerate.out" ]; then
    perl /opt/EVidenceModeler/EvmUtils/misc/Exonerate_to_evm_gff3.pl reference_exonerate.out > rexonerate.gff3
    sed -i -e 's/exonerate/homology_exonerate/g' rexonerate.gff3
    cat rexonerate.gff3 >> protein_alignments.gff3
    echo -e "PROTEIN\thomology_exonerate\t4" >> weights.txt
fi

if [ -f "genomethreader.gff3" ]; then
    perl /opt/EVidenceModeler/EvmUtils/misc/genomeThreader_to_evm_gff3.pl genomethreader.gff3 > egenomethreader.gff3
    cat egenomethreader.gff3 >> gene_predictions.gff3
    echo -e "OTHER_PREDICTION\tgenomeThreader\t4" >> weights.txt
fi
if [ -f "reference_genomethreader.gff3" ]; then
    perl /opt/EVidenceModeler/EvmUtils/misc/genomeThreader_to_evm_gff3.pl reference_genomethreader.gff3 > ereference_genomethreader.gff3
    sed -i -e 's/genomeThreader/homology_genomeThreader/g' ereference_genomethreader.gff3
    cat ereference_genomethreader.gff3 >> gene_predictions.gff3
    echo -e "OTHER_PREDICTION\thomology_genomeThreader\t4" >> weights.txt
fi

if [ -f "miniprot.gtf" ]; then
    perl /opt/EVidenceModeler/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl miniprot.gtf > miniprot.gff3
    sed -i -e 's/Augustus/miniprot/g' miniprot.gff3
    cat miniprot.gff3 >> gene_predictions.gff3
    echo -e "OTHER_PREDICTION\tminiprot\t4" >> weights.txt
fi
if [ -f "reference_miniprot.gff3" ]; then
    perl /opt/EVidenceModeler/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl reference_miniprot.gtf > reference_miniprot.gff3
    sed -i -e 's/Augustus/homology_miniprot/g' reference_miniprot.gff3
    cat reference_miniprot.gff3 >> gene_predictions.gff3
    echo -e "OTHER_PREDICTION\thomology_miniprot\t4" >> weights.txt
fi

if [ -f "prosplign.gff3" ]; then
    cat prosplign.gff3 >> protein_alignments.gff3
    echo -e "PROTEIN\tProsplign\t4" >> weights.txt
fi
if [ -f "reference_prosplign.gff3" ]; then
    cat reference_prosplign.gff3 >> protein_alignments.gff3
    echo -e "PROTEIN\tHOMOLOGY_Prosplign\t4" >> weights.txt
fi
if [ -f "splign.gff3" ]; then
    cat splign.gff3 >> transcript_alignments.gff3
    echo -e "TRANSCRIPT\tRefSeq\t6" >> weights.txt
fi
if [ -f "reference_splign.gff3" ]; then
    cat reference_splign.gff3 >> transcript_alignments.gff3
    echo -e "TRANSCRIPT\tHOMOLOGY_RefSeq\t6" >> weights.txt
fi
if [ -f "pasa_assembly.gff3" ]; then
    cat pasa_assembly.gff3 >> transcript_alignments.gff3
    echo -e "TRANSCRIPT\tpasa\t6" >> weights.txt
fi
if [ -f "reference_pasa_assembly.gff3" ]; then
    cat reference_pasa_assembly.gff3 >> transcript_alignments.gff3
    echo -e "TRANSCRIPT\tHOMOLOGY_pasa\t6" >> weights.txt
fi
if [ -f "transdecoder.gff3" ]; then
    cat transdecoder.gff3 >> gene_predictions.gff3
    echo -e "OTHER_PREDICTION\ttransdecoder\t1" >> weights.txt
fi
if [ -f "reference_transdecoder.gff3" ]; then
    cat reference_transdecoder.gff3 >> gene_predictions.gff3
    echo -e "OTHER_PREDICTION\tHOMOLOGY_transdecoder\t1" >> weights.txt
fi

cp weights.txt /data/

#run EVM
bash ${baseDir}/runEVM.sh -s ${size} -z ${threads}
mkdir EVMRun
cp EVM.all.gff3 EVMRun/
cd EVMRun
agat_sp_extract_sequences.pl --clean_final_stop --gff EVM.all.gff3 -f ../genome.fasta -p -o protein.fa
conda activate BUSCO
busco -i protein.fa -l ${lineage} -o buscoout -m protein -c ${threads}
conda deactivate
mv buscoout/run_${lineage}/full_table.tsv .
bash ${baseDir}/annotationFilter.sh -i EVM.all.gff3 -o False -p EVM
mv EVMBuscos.gff3 ..
mv EVMRemoved.gff3 ..
cd ..
cat EVMBuscos.gff3 >> completeBuscos.gff3
    

#filter out duplicate BUSCOs and keep only the best ones
mkdir filter0
cp completeBuscos.gff3 filter0/
cd filter0/
agat_sp_extract_sequences.pl --clean_final_stop --gff completeBuscos.gff3 -f ../genome.fasta -p -o proteinOriginalCompleteBuscos.fa
conda activate BUSCO
busco -i proteinOriginalCompleteBuscos.fa -l ${lineage} -o buscoout -m protein -c ${threads}
conda deactivate
cp completeBuscos.gff3 ../completeBuscosUnfiltered.gff3
cp buscoout/run_${lineage}/full_table.tsv .

echo "filtering out duplicates round 1"
bash ${baseDir}/filterDuplicateBuscos.sh -f Duplicated
cp completeBuscos.gff3 ../completeBuscosNoDups.gff3

echo "filtering out duplicates round 2"

rm full_table*.tsv
agat_sp_extract_sequences.pl --clean_final_stop --gff completeBuscos.gff3 -f ../genome.fasta -p -o proteinCompleteBuscosNoDups.fa
conda activate BUSCO
busco -i proteinCompleteBuscosNoDups.fa -l ${lineage} -o buscooutnoDups -m protein -c ${threads}
conda deactivate
cp buscooutnoDups/run_${lineage}/full_table.tsv .

bash ${baseDir}/filterDuplicateBuscos.sh -f Missing

mv filteredCompleteBuscos.gff3 ../completeBuscos.gff3
cd ..

#run gfacs
mkdir split
cp genome.fasta split/genome.fa

cd split/

awk '/>/{n++}{print >"___"  n ".fasta" }' "genome.fa"

for file in *.fasta; do
    sed -i '1s/\s.*$//' $file
done

for file in *.fasta; do
    cat $file >> formattedgenome.fa
done

cp formattedgenome.fa ..

cd ..

/opt/gFACs-master/gFACs.pl -f EVM_1.1.1_gff3 -p gFacs_Filtered --statistics --statistics-at-every-step --min-CDS-size 250 --unique-genes-only --fasta formattedgenome.fa --splice-table --allow-alternate-starts --nt-content --create-gtf --create-simple-gtf --create-gff3 --annotated-all-genes-only  --compatibility SnpEff EVM_1.1.1_gene_prediction EVM_1.1.1_alignment -O . EVMRemoved.gff3


cat completeBuscos.gff3 > unsortedAnnotation.gff3

echo "sorting lenient"

cat unsortedAnnotation.gff3 >> FinalStructuralAnnotationLenientFilter.gff3
cat EVMRemoved.gff3 >> FinalStructuralAnnotationLenientFilter.gff3

cat unsortedAnnotation.gff3 >> FinalStructuralAnnotationStrictFilter.gff3
cat gFacs_Filtered_EVM_1.1.1_gene_prediction_format.gff3 >> FinalStructuralAnnotationStrictFilter.gff3

#echo -e "mkdir sortingLenient/" >> sortLenient.txt
#echo -e "cp unsortedAnnotation.gff3 sortingLenient/" >> sortLenient.txt
#echo -e "cp genome.fasta sortingLenient/genome.fa" >> sortLenient.txt
#echo -e "cd sortingLenient/" >> sortLenient.txt
#echo -e "cat ../EVMRemoved.gff3 >> unsortedAnnotation.gff3" >> sortLenient.txt
#echo -e "bash ${baseDir}/sortAnnotation.sh -g genome.fa -a unsortedAnnotation.gff3" >> sortLenient.txt
#echo -e "cp sorted.gff3 ../FinalStructuralAnnotationLenientFilter.gff3" >> sortLenient.txt
#echo -e "cd .." >> sortLenient.txt

#echo "sorting Strict"

#echo -e "mkdir sortingStrict/" >> sortStrict.txt
#echo -e "cp unsortedAnnotation.gff3 sortingStrict/" >> sortStrict.txt
#echo -e "cp genome.fasta sortingStrict/genome.fa" >> sortStrict.txt
#echo -e "cd sortingStrict/" >> sortStrict.txt
#echo -e "cat ../gFacs_Filtered_EVM_1.1.1_gene_prediction_format.gff3 >> unsortedAnnotation.gff3" >> sortStrict.txt
#echo -e "bash ${baseDir}/sortAnnotation.sh -g genome.fa -a unsortedAnnotation.gff3" >> sortStrict.txt
#echo -e "cp sorted.gff3 ../FinalStructuralAnnotationStrictFilter.gff3" >> sortStrict.txt
#echo -e "cd .." >> sortStrict.txt

#echo -e "bash sortStrict.txt" >> sortCommands.txt
#echo -e "bash sortLenient.txt" >> sortCommands.txt

#parallel < sortCommands.txt > sorting.log

cp FinalStructuralAnnotationStrictFilter.gff3 FinalStructuralAnnotationStrictFilterold.gff3
cp FinalStructuralAnnotationLenientFilter.gff3 FinalStructuralAnnotationLenientFilterold.gff3

gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o tmp_tRNAScan.gff3 tRNAScan.gff3
mv tmp_tRNAScan.gff3 tRNAScan.gff3

sed -i '/^#/d' tRNAScan.gff3

#tidying up and converting to gtf
echo "Tidying up the final gff3 annotations, converting to gtf, and running final agat stats"
cat tRNAScan.gff3 >> FinalStructuralAnnotationLenientFilter.gff3
echo -e "mkdir tidyLenient" >> tidyLenient.txt
echo -e "mv FinalStructuralAnnotationLenientFilter.gff3 tidyLenient/" >> tidyLenient.txt
echo -e "cd tidyLenient" >> tidyLenient.txt
echo -e "cp ../genome.fasta ." >> tidyLenient.txt
echo -e "bash ${baseDir}/gff3togtf.sh -i FinalStructuralAnnotationLenientFilter.gff3" >> tidyLenient.txt
echo -e "agat_sp_extract_sequences.pl --clean_final_stop --gff FinalStructuralAnnotationLenientFilter.gff3 -f genome.fasta -p -o proteinFinalLenient.fa" >> tidyLenient.txt
echo -e "agat_sp_statistics.pl --gff FinalStructuralAnnotationLenientFilter.gff3 --output FinalStructuralAnnotationLenientFilter.gff3.stats" >> tidyLenient.txt
echo -e "cp * .." >> tidyLenient.txt
echo -e "cd .." >> tidyLenient.txt

echo -e "mkdir tidyStrict" >> tidyStrict.txt
echo -e "mv FinalStructuralAnnotationStrictFilter.gff3 tidyStrict/" >> tidyStrict.txt
echo -e "cd tidyStrict" >> tidyStrict.txt
echo -e "cp ../genome.fasta ." >> tidyStrict.txt
#echo -e "bash ${baseDir}/gff3togtf.sh -i FinalStructuralAnnotationStrictFilter.gff3" >> tidyStrict.txt
echo -e "agat_sp_extract_sequences.pl --clean_final_stop --gff FinalStructuralAnnotationStrictFilter.gff3 -f genome.fasta -p -o proteinFinalStrict.fa" >> tidyStrict.txt
echo -e "agat_sp_statistics.pl --gff FinalStructuralAnnotationStrictFilter.gff3 --output FinalStructuralAnnotationStrictFilter.gff3.stats" >> tidyStrict.txt
echo -e "cp * .." >> tidyStrict.txt
echo -e "cd .." >> tidyStrict.txt

echo -e "bash tidyStrict.txt" >> tidyCommands.txt
echo -e "bash tidyLenient.txt" >> tidyCommands.txt

parallel < tidyCommands.txt > tidy.log

#running final busco stats
echo "running final busco stats"

conda activate BUSCO
busco -i proteinFinalLenient.fa -l ${lineage} -o buscooutLenient -m protein -c ${threads}
busco -i proteinFinalStrict.fa -l ${lineage} -o buscooutStrict -m protein -c ${threads}
mv buscooutLenient/*.txt .
mv buscooutStrict/*.txt .
conda deactivate
