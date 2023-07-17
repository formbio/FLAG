#!/bin/bash
#evidencemodeler.sh

usage() {
  echo "-h Help documentation for liftoff.sh"
  echo "-i  --input genome file in fasta format"
  echo "-a  --algo"
  echo "-l  --lib file"
  echo "-t  --threads"
  echo "Example: bash evidencemodeler.sh -i genome.fa -l -t "
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:a:l:t:h opt
do
    case $opt in
        i) genome=$OPTARG;;
        a) algo=$OPTARG;;
        l) library=$OPTARG;;
        t) threads=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${genome} ]]
then
    usage
fi
if [[ -z $threads ]]
then
    threads=`nproc`
fi

#process the inputs and format them for EVM and make the weights.txt file
mv ${genome} genome.fasta
prefixAssembly="${genome%.fa*}"
if [[ "$algo" = "RepeatMasker" ]]; then
    libCommand=""
    if [ $library ]; then
        libCommand="-lib ${library}"
    fi
    RepeatMasker ${libCommand} -no_is -gff -xsmall -pa ${threads} genome.fasta
    mv genome.fasta.out.gff ${prefixAssembly}_Repeatmaskermasked.fa.out.gff
    mv genome.fasta.masked ${prefixAssembly}_Repeatmaskermasked.fa
elif [[ "$algo" == "WindowMasker" ]]; then
    threads=`nproc`
    declare -i count=$(grep -o ">" genome.fa | wc -l)
    if [ "$count" != "1" ]; then
        awk '/>/{n++}{print >"___"  n ".fasta" }' "genome.fa"
        rm genome.fa
    fi

    touch commands.txt
    for FILE in ./___*.fasta; do
        back=${FILE%.fasta}
        count=${back:5}
        touch command_${count}.sh
        
        # First step of windowmasker
        echo -e "windowmasker -mk_counts -in ${FILE:2} -checkdup true -out wm_counts_${count}" >> command_${count}.sh
        #windowmasker -mk_counts -in $genome -out intermediate.txt

        # Seconds step of windowmasker
        echo -e "windowmasker -ustat wm_counts_${count} -in ${FILE:2} -outfmt fasta -out maskedGenome_${count}.fa -dust true" >> command_${count}.sh
        #windowmasker -ustat intermediate.txt
        
        echo -e "bash command_${count}.sh" >> commands.txt
    done

    parallel --jobs $threads < commands.txt

    touch maskedGenome.fa
    for (( i=1; i<=$count; i++ )); do
        echo $i
        cat maskedGenome_${i}.fa >> maskedGenome.fa
    done
elif [[ "$algo" == "RepeatModeler" ]]; then
    BuildDatabase -name repeatmodeldatabase genome.fasta
    RepeatModeler -database repeatmodeldatabase -pa ${threads}
    mv */consensi.fa.classified ${prefixAssembly}_consensirepeatmodelermasked.fa
elif [[ "$algo" == "RepeatMasker_with_RepeatModeler" ]]; then
    libCommand=""
    if [ $library ]; then
        libCommand="-lib ${library}"
    fi
    mkdir firstRound
    cd firstRound
    cp ../genome.fasta .
    RepeatMasker ${libCommand} -no_is -gff -xsmall -pa ${threads} genome.fasta
    mv genome.fasta.out.gff ../RepeatmaskermaskedFirstRound.fa.out.gff
    mv genome.fasta.masked ../RepeatmaskermaskedFirstRound.fa
    cd ..
    mkdir RepeatModel
    cd RepeatModel
    cp ../RepeatmaskermaskedFirstRound.fa .
    BuildDatabase -name repeatmodeldatabase Repeatmaskermasked.fa
    RepeatModeler -database repeatmodeldatabase -pa ${threads}
    mv */consensi.fa.classified ../consensirepeatmodelermasked.fa
    cd ..
    if [[ -s consensirepeatmodelermasked.fa ]]; then
        RepeatMasker --lib consensirepeatmodelermasked.fa -no_is -gff -xsmall -pa ${threads} genome.fasta
        mv genome.fasta.out.gff ${prefixAssembly}_Repeatmaskermasked.fa.out.gff
        mv genome.fasta.masked ${prefixAssembly}_Repeatmaskermasked.fa
    else
        mv RepeatmaskermaskedFirstRound.fa.out.gff ${prefixAssembly}_Repeatmaskermasked.fa.out.gff
        mv RepeatmaskermaskedFirstRound.fa ${prefixAssembly}_Repeatmaskermasked.fa
    fi
fi
