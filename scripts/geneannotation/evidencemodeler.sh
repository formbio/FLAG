#!/bin/bash
#evidencemodeler.sh

usage() {
  echo "-h Help documentation for liftoff.sh"
  echo "-i  --input genome file in fasta format"
  echo "-b  --Run busco True or False"
  echo "-l  --lineage"
  echo "-z  --threads"
  echo "Example: bash evidencemodeler.sh -i genome.fa"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:b:l:z:h opt
do
    case $opt in
        i) genome=$OPTARG;;
        b) buscorun=$OPTARG;;
        l) lineage=$OPTARG;;
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
if [[ -z ${buscorun} ]]
then
    buscorun="True"
fi
if [[ -z ${lineage} ]]
then
    buscorun="False"
fi
if [[ -z $threads ]]
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

#process the inputs and format them for EVM and make the weights.txt file
mv ${genome} genome.fasta
touch weights.txt
touch gene_predictions.gff3
touch protein_alignments.gff3
touch transcript_alignments.gff3
touch other_predictions.gff3


if [ -f "filtered_augustus.gtf" ]; then
    if [[ "$buscorun" == "True" ]]; then
        gffread filtered_augustus.gtf -y filtered_augustus.fa -g genome.fasta
        conda activate BUSCO
        busco -i filtered_augustus.fa -l ${lineage} -o buscoout -m protein -c ${threads}
        mv buscoout/run_${lineage}/full_table.tsv .
        conda deactivate
        bash ${baseDir}/annotationFilter.sh -i filtered_augustus.gtf -s 101 -m 101 -e 101 -l 99999999999 -o True
        mv filteredAnnots.gtf complete_augustus.gtf
        mv removedAnnots.gtf removed_augustus.gtf
        perl /opt/EVidenceModeler-1.1.1/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl complete_augustus.gtf > complete_augustus.gff3
        perl /opt/EVidenceModeler-1.1.1/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl removed_augustus.gtf > removed_augustus.gff3
        sed -i -e 's/Augustus/Complete_Augustus/g' complete_augustus.gff3
        sed -i -e 's/Augustus/Removed_Augustus/g' removed_augustus.gff3
        cat complete_augustus.gff3 >> gene_predictions.gff3
        cat removed_augustus.gff3 >> gene_predictions.gff3
        echo -e "ABINITIO_PREDICTION\tComplete_Augustus\t10" >> weights.txt
        echo -e "ABINITIO_PREDICTION\tRemoved_Augustus\t4" >> weights.txt
        rm -rf buscoout
        rm full_table.tsv
    else
        perl /opt/EVidenceModeler-1.1.1/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl filtered_augustus.gtf > augustus.gff3
        cat augustus.gff3 >> gene_predictions.gff3
        echo -e "ABINITIO_PREDICTION\tAugustus\t5" >> weights.txt
    fi
fi

if [ -f "filtered_reference_augustus.gtf" ]; then
    if [[ "$buscorun" == "True" ]]; then
        gffread filtered_reference_augustus.gtf -y filtered_reference_augustus.fa -g genome.fasta
        conda activate BUSCO
        busco -i filtered_reference_augustus.fa -l ${lineage} -o buscoout -m protein -c ${threads}
        mv buscoout/run_${lineage}/full_table.tsv .
        conda deactivate
        bash ${baseDir}/annotationFilter.sh -i filtered_reference_augustus.gtf -s 101 -m 101 -e 101 -l 99999999999 -o True
        mv filteredAnnots.gtf complete_reference_augustus.gtf
        mv removedAnnots.gtf removed_reference_augustus.gtf
        perl /opt/EVidenceModeler-1.1.1/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl complete_reference_augustus.gtf > complete_reference_augustus.gff3
        perl /opt/EVidenceModeler-1.1.1/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl removed_reference_augustus.gtf > removed_reference_augustus.gff3
        sed -i -e 's/Augustus/Complete_RefAugustus/g' complete_reference_augustus.gff3
        sed -i -e 's/Augustus/Removed_RefAugustus/g' removed_reference_augustus.gff3
        cat complete_reference_augustus.gff3 >> gene_predictions.gff3
        cat removed_reference_augustus.gff3 >> gene_predictions.gff3
        echo -e "ABINITIO_PREDICTION\tComplete_RefAugustus\t10" >> weights.txt
        echo -e "ABINITIO_PREDICTION\tRemoved_RefAugustus\t4" >> weights.txt
        rm -rf buscoout
        rm full_table.tsv
    else
        perl /opt/EVidenceModeler-1.1.1/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl filtered_reference_augustus.gtf > reference_augustus.gff3
        sed -i -e 's/Augustus/HOMOLOGY_Augustus/g' reference_augustus.gff3
        cat reference_augustus.gff3 >> gene_predictions.gff3
        echo -e "ABINITIO_PREDICTION\tHOMOLOGY_Augustus\t2" >> weights.txt
    fi
fi

if [ -f "exonerate.out" ]; then
    perl /opt/EVidenceModeler-1.1.1/EvmUtils/misc/Exonerate_to_evm_gff3.pl exonerate.out > exonerate.gff3
    cat exonerate.gff3 >> protein_alignments.gff3
    echo -e "PROTEIN\texonerate\t4" >> weights.txt
fi

if [ -f "reference_exonerate.out" ]; then
    perl /opt/EVidenceModeler-1.1.1/EvmUtils/misc/Exonerate_to_evm_gff3.pl reference_exonerate.out > rexonerate.gff3
    sed -i -e 's/exonerate/homology_exonerate/g' rexonerate.gff3
    cat rexonerate.gff3 >> protein_alignments.gff3
    echo -e "PROTEIN\thomology_exonerate\t2" >> weights.txt
fi

if [ -f "prosplign.gff3" ]; then
    cat prosplign.gff3 >> protein_alignments.gff3
    echo -e "PROTEIN\tProsplign\t5" >> weights.txt
fi
if [ -f "reference_prosplign.gff3" ]; then
    cat reference_prosplign.gff3 >> protein_alignments.gff3
    echo -e "PROTEIN\tHOMOLOGY_Prosplign\t2" >> weights.txt
fi
if [ -f "splign.gff3" ]; then
    cat splign.gff3 >> transcript_alignments.gff3
    echo -e "TRANSCRIPT\tRefSeq\t10" >> weights.txt
fi
if [ -f "reference_splign.gff3" ]; then
    cat reference_splign.gff3 >> transcript_alignments.gff3
    echo -e "TRANSCRIPT\tHOMOLOGY_RefSeq\t3" >> weights.txt
fi
if [ -f "pasa_assembly.gff3" ]; then
    cat pasa_assembly.gff3 >> transcript_alignments.gff3
    echo -e "TRANSCRIPT\tpasa\t7" >> weights.txt
fi
if [ -f "reference_pasa_assembly.gff3" ]; then
    cat reference_pasa_assembly.gff3 >> transcript_alignments.gff3
    echo -e "TRANSCRIPT\tHOMOLOGY_pasa\t3" >> weights.txt
fi
if [ -f "transdecoder.gff3" ]; then
    cat transdecoder.gff3 >> gene_predictions.gff3
    echo -e "OTHER_PREDICTION\ttransdecoder\t2" >> weights.txt
fi
if [ -f "reference_transdecoder.gff3" ]; then
    cat reference_transdecoder.gff3 >> gene_predictions.gff3
    echo -e "OTHER_PREDICTION\tHOMOLOGY_transdecoder\t1" >> weights.txt
fi
if [ -f "liftoff.gff3" ]; then
    cat liftoff.gff3 >> gene_predictions.gff3
    echo -e "OTHER_PREDICTION\tLiftoff\t10" >> weights.txt
fi

cp weights.txt /data/

#the below line is for the simple case
#perl /opt/EVidenceModeler-1.1.1/evidence_modeler.pl --genome genome.fasta --weights weights.txt --gene_predictions gene_predictions.gff3 --protein_alignments protein_alignments.gff3 --transcript_alignments transcript_alignments.gff3 > evm.out
#the above command will output the EVM output file called evm.out


#the below is the advanced configuration
#partition the inputs
#perl /opt/EVidenceModeler-1.1.1/EvmUtils/partition_EVM_inputs.pl --genome genome.fasta \
#     --gene_predictions gene_predictions.gff3 --protein_alignments protein_alignments.gff3 \
#     --transcript_alignments transcript_alignments.gff3 \
#     --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out
     
perl /opt/EVidenceModeler-1.1.1/EvmUtils/partition_EVM_inputs.pl --genome genome.fasta --gene_predictions gene_predictions.gff3 --protein_alignments protein_alignments.gff3 --transcript_alignments transcript_alignments.gff3 --segmentSize 1000000 --overlapSize 200000 --partition_listing partitions_list.out

#generate evm command set
perl /opt/EVidenceModeler-1.1.1/EvmUtils/write_EVM_commands.pl --genome genome.fasta --weights /data/weights.txt \
      --gene_predictions gene_predictions.gff3 --protein_alignments protein_alignments.gff3 \
      --transcript_alignments transcript_alignments.gff3 \
      --output_file_name evm.out  --partitions partitions_list.out >  commands.list
      
#run the generated commands
parallel -j ${threads} < commands.list > evm.log

#combine the partitions
perl /opt/EVidenceModeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out

#convert the evm.out files in each contig directory to gff3
perl /opt/EVidenceModeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome genome.fasta

#combine each gff3 file in each contig directory to a single gff3 file
find . -regex ".*evm.out.gff3" -exec cat {} \; > EVM.all.gff3

cat tRNAScan.gff3 >> EVM.all.gff3


#perl /opt/EVidenceModeler-1.1.1/EvmUtils/partition_EVM_inputs.pl --genome genome.fasta  --gene_predictions gene_predictions.gff3 --protein_alignments protein_alignments.gff3 --transcript_alignments transcript_alignments.gff3 --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out


SZQC01000071.1-g13159
grep "SZQC01000071.1-g13159" test.gff3 | grep -P "Augustus\tgene"

7226.t1

116.t1

grep -E "\-116;Name" test.gff3 | grep -P "Augustus\tgene"
