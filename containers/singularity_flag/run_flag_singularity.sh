#!/bin/bash
#run_flag_singularity.sh

usage() {
  echo "-h Help documentation for run_flag_singularity.sh"
  echo "-g  --input genome file in fasta format"
  echo "-r  --input rna or transcript file in fasta format"
  echo "-p  --input protein file in fasta format"
  echo "-m  --masker"
  echo "-t  --isTranscript"
  echo "-f  --fafile"
  echo "-a  --annotation file in gtf or gff3 format"
  echo "-l  --lineage"
  echo "-z  --annotationalgo such as Helixer,helixer_trained_augustus"
  echo "-q  --helixerModel"
  echo "-s  --size (normal or small)"
  echo "-n  --speciesScientificName"
  echo "-w  --proteinalgo"
  echo "-y  --runMode"
  echo "-u  --profile"
  echo "-o  --output"
  echo "Example: bash run_flag_singularity.sh -g Erynnis_tages-GCA_905147235.1-softmasked.fa -r curatedButterflyRNA.fa -p curatedButterflyProteins.fa -m skip -t true -l lepidoptera_odb10 -z Helixer,helixer_trained_augustus -q vertebrate -s small -n Eynnis_tages -w miniprot -y normal -u singularity -o outputdir"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :g:r:p:m:t:f:a:l:z:q:s:n:w:y:u:o:h opt
do
    case $opt in
        g) genome=$OPTARG;;
        r) rna=$OPTARG;;
        p) protein=$OPTARG;;
        m) masker=$OPTARG;;
        t) isTranscript=$OPTARG;;
        f) fafile=$OPTARG;;
        a) gtffile=$OPTARG;;
        l) lineage=$OPTARG;;
        z) annotationalgo=$OPTARG;;
        q) helixerModel=$OPTARG;;
        s) size=$OPTARG;;
        n) speciesScientificName=$OPTARG;;
        w) proteinalgo=$OPTARG;;
        y) runMode=$OPTARG;;
        u) profile=$OPTARG;;
        o) output=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${genome} ]]
then
    usage
fi
if [[ -z ${rna} ]]
then
    echo "No RNA or Transcript File Provided"
    rnaCommand=""
    transcriptInCommand=""
    externalalgoToUse=""
else
    rnaCommand="--rna /data/${rna}"
    externalalgoToUse="input_transcript"
    if [[ -z ${isTranscript} ]]
    then
        echo "isTranscript was not specified. Assuming it is a transcript"
        transcriptInCommand="--transcriptIn true"
    else
        transcriptInCommand="--transcriptIn ${isTranscript}"
    fi
fi
if [[ -z ${protein} ]]
then
    echo "No protein File Provided"
    proteinCommand=""
else
    proteinCommand="--proteins /data/${protein}"
    externalalgoToUse="${externalalgoToUse},input_proteins"
fi
if [[ -z ${fafile} ]]
then
    echo "No reference fafile File for liftoff Provided"
    fafileCommand=""
else
    fafileCommand="--fafile /data/${fafile}"
fi
if [[ -z ${gtffile} ]]
then
    echo "No reference annotation File for liftoff Provided"
    gtffileCommand=""
else
    gtffileCommand="--gtffile /data/${gtffile}"
fi
if [[ -z ${masker} ]]
then
    echo "Masking was not specified. Assuming the genome is masked"
    maskingCommand="--masker skip"
else
    maskingCommand="--masker ${masker}"
fi
if [[ -z ${lineage} ]]
then
    echo "Lineage was not specified. It is very much recommended to select a lineage that is specific as possible. Due to this though we are going to use eukaryota_odb10"
    lineageCommand="--lineage eukaryota_odb10"
else
    lineageCommand="--lineage ${lineage}"
fi
if [[ -z ${annotationalgo} ]]
then
    echo "annotationalgo was not specified. Defaulting to Helixer,helixer_trained_augustus"
    annotationalgoCommand="--annotationalgo Helixer,helixer_trained_augustus"
else
    annotationalgoCommand="--annotationalgo ${annotationalgo}"
fi
if [[ -z ${helixerModel} ]]
then
    echo "helixerModel was not specified. Defaulting to vertebrate"
    helixerModelCommand="--helixerModel vertebrate"
else
    helixerModelCommand="--helixerModel ${helixerModel}"
fi
if [[ -z ${size} ]]
then
    echo "size was not specified. Defaulting to normal"
    sizeCommand="--size normal"
else
    sizeCommand="--size ${size}"
fi
if [[ -z ${proteinalgo} ]]
then
    echo "proteinalgo was not specified. Defaulting to miniprot"
    proteinalgoCommand="--proteinalgo miniprot"
else
    proteinalgoCommand="--proteinalgo ${proteinalgo}"
fi
if [[ -z ${speciesScientificName} ]]
then
    echo "speciesScientificName was not specified. Defaulting to Replace_species"
    speciesScientificNameCommand="--speciesScientificName Replace_species"
else
    speciesScientificNameCommand="--speciesScientificName ${speciesScientificName}"
fi
if [[ -z ${runMode} ]]
then
    echo "runMode was not specified. Defaulting to normal"
    runModeCommand="--runMode normal"
else
    runModeCommand="--runMode ${runMode}"
fi
if [[ -z ${profile} ]]
then
    echo "profile was not specified. Defaulting to singularity"
    profileCommand="-profile singularity"
else
    profileCommand="-profile ${profile}"
fi
if [[ ${externalalgoToUse} != "" ]]; then
    externalalgoCommand="--externalalgo ${externalalgoToUse}"
else
    externalalgoCommand=""
fi
if [[ -z ${output} ]]; then
    echo "output is outputdir"
    output="outputdir"
else
    echo "outputdir is ${output}"
fi

echo "Formatting outputdir"
mkdir /data/$output
touch "/data/${output}/emptyProteinPlaceHolder.txt"
touch "/data/${output}/emptyTranscriptPlaceHolder.txt"
touch "/data/${output}/emptyPlaceHolder.txt"

cd /opt/FLAG

echo "RUNNING"
echo "nextflow run main.nf -w /data/workdir/ --output /data/${output}/ --genome /data/${genome} ${rnaCommand} ${proteinCommand} ${fafileCommand} ${gtffileCommand} ${maskingCommand} ${transcriptInCommand} ${lineageCommand} ${annotationalgoCommand} ${helixerModelCommand} ${externalalgoCommand} ${sizeCommand} ${proteinalgo} ${speciesScientificName} --funcAnnotProgram eggnog --eggnogDB /opt/FLAG/eggnogDB.tar.gz ${runModeCommand} ${profileCommand}"

nextflow run main.nf -w /data/workdir/ --output /data/${output}/ \
--genome /data/${genome} \
${rnaCommand} ${proteinCommand} \
${fafileCommand} ${gtffileCommand} \
${maskingCommand} ${transcriptInCommand} ${lineageCommand} \
${annotationalgoCommand} ${helixerModelCommand} \
${externalalgoCommand} ${sizeCommand} ${proteinalgoCommand} \
${speciesScientificName} --funcAnnotProgram eggnog --eggnogDB /opt/FLAG/eggnogDB.tar.gz \
${runModeCommand} ${profileCommand}
