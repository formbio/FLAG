#!/bin/bash
#pasa.sh

usage() {
  echo "-h Help documentation for liftoff.sh"
  echo "-g  --input genome file in fasta format"
  echo "-r  --input transcript file in fasta format"
  echo "-t  --threads"
  echo "Example: bash pasa.sh -i transcript.fa -t threads"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :g:r:t:h opt
do
    case $opt in
        g) genome=$OPTARG;;
        r) transcript=$OPTARG;;
        t) threads=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${genome} ]] || [[ -z ${transcript} ]]
then
    usage
fi
if [[ -z ${threads} ]]
then
    threads=16
fi

Launch_PASA_pipeline.pl -c /usr/local/src/PASApipeline/sample_data/sqlite.confs/alignAssembly.config -C -R --ALIGNER blat -g $genome -t $transcript --CPU $threads

rm pblat_outdir/*.gff3
mv pblat_outdir/*.pslx.top_1 pasa_pblat.pslx

if [ "$transcript" == "reference_trinity.fasta" ] || [ "$transcript" == "formatted_referenceRNA.fa" ]
then
    mv sample_mydb_pasa.sqlite.pasa_assemblies.gff3 reference_pasa_assembly.gff3
    mv pasa_pblat.pslx reference_pasa_pblat.pslx
    sed -i -e 's/assembler-sample_mydb_pasa.sqlite/HOMOLOGY_pasa/g' reference_pasa_assembly.gff3
else
    mv sample_mydb_pasa.sqlite.pasa_assemblies.gff3 pasa_assembly.gff3
    sed -i -e 's/assembler-sample_mydb_pasa.sqlite/pasa/g' pasa_assembly.gff3
fi
