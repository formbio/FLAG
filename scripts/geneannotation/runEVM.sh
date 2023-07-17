#!/bin/bash
#evidencemodeler.sh

usage() {
  echo "-h Help documentation for liftoff.sh"
  echo "-s  --genome size (small or normal)"
  echo "-z  --threads"
  echo "Example: bash evidencemodeler.sh -i genome.fa"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :s:z:h opt
do
    case $opt in
        s) size=$OPTARG;;
        z) threads=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $threads ]]
then
    threads=`nproc`
fi
if [[ -z $size ]]
then
    size="normal"
fi

#the below line is for the simple case
#perl /opt/EVidenceModeler/evidence_modeler.pl --genome genome.fasta --weights weights.txt --gene_predictions gene_predictions.gff3 --protein_alignments protein_alignments.gff3 --transcript_alignments transcript_alignments.gff3 > evm.out
#the above command will output the EVM output file called evm.out


#the below is the advanced configuration
#partition the inputs
#perl /opt/EVidenceModeler/EvmUtils/partition_EVM_inputs.pl --genome genome.fasta \
#     --gene_predictions gene_predictions.gff3 --protein_alignments protein_alignments.gff3 \
#     --transcript_alignments transcript_alignments.gff3 \
#     --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out
#7x     perl /opt/EVidenceModeler/EvmUtils/partition_EVM_inputs.pl --genome genome.fasta --gene_predictions gene_predictions.gff3 --protein_alignments protein_alignments.gff3 --transcript_alignments transcript_alignments.gff3 --segmentSize 125000 --overlapSize 25000 --partition_listing partitions_list.out
if [[ "$size" == "small" ]]; then
    perl /opt/EVidenceModeler/EvmUtils/partition_EVM_inputs.pl --genome genome.fasta --gene_predictions gene_predictions.gff3 --protein_alignments protein_alignments.gff3 --transcript_alignments transcript_alignments.gff3 --segmentSize 150000 --overlapSize 30000 --partition_listing partitions_list.out --partition_dir .

else
    perl /opt/EVidenceModeler/EvmUtils/partition_EVM_inputs.pl --genome genome.fasta --gene_predictions gene_predictions.gff3 --protein_alignments protein_alignments.gff3 --transcript_alignments transcript_alignments.gff3 --segmentSize 1000000 --overlapSize 200000 --partition_listing partitions_list.out --partition_dir .
fi

#generate evm command set
perl /opt/EVidenceModeler/EvmUtils/write_EVM_commands.pl --genome genome.fasta --weights /data/weights.txt \
      --gene_predictions gene_predictions.gff3 --protein_alignments protein_alignments.gff3 \
      --transcript_alignments transcript_alignments.gff3 \
      --output_file_name evm.out  --partitions partitions_list.out >  commands.list
      
#run the generated commands
parallel -j ${threads} < commands.list > evm.log

#combine the partitions
perl /opt/EVidenceModeler/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out

#convert the evm.out files in each contig directory to gff3
perl /opt/EVidenceModeler/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome genome.fasta

#combine each gff3 file in each contig directory to a single gff3 file
find . -regex ".*evm.out.gff3" -exec cat {} \; > EVM.all.gff3

