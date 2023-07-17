#!/bin/bash
#augustus_full.sh

usage() {
  echo "-h Help documentation for tasser.sh"
  echo "-g  --input genome file"
  echo "-r  --input transcript file"
  echo "-s  --input species"
  echo "-e  --input exonerate.out file"
  echo "-x  --input reference_exonerate.out file"
  echo "-p  --Is this a pretrained species? (True or False)"
  echo "-n  --number of cpus"
  echo "Example: bash augustus.sh -s human -g GRCh38_latest_genomic.fna -u on"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :g:r:s:e:p:n:h opt
do
    case $opt in
        g) genome=$OPTARG;;
        r) transcript=$OPTARG;;
        s) species=$OPTARG;;
        e) exoneratehints=$OPTARG;;
        p) pretrained=$OPTARG;;
        n) cpu=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${transcript} ]] || [[ -z ${genome} ]]
then
    usage
fi
if [[ -z ${species} ]]
then
    species='current'
fi
if [[ -z ${cpu} ]]
then
    cpu=`nproc`
fi
if [[ -z ${pretrained} ]]
then
    pretrained="False"
fi
mv $genome genome.fa

transcriptExtension="${transcript: -4}"
echo "Last 3 characters of transcript: ${transcriptExtension}"

if [[ "${transcriptExtension}" == ".gff" ]] || [[ "${transcriptExtension}" == ".gtf" ]] || [[ "${transcriptExtension}" == "gff3" ]]; then
    mv $transcript trainingannotations.gff
    mv $exoneratehints traininggenome.fa
    gff2gbSmallDNA.pl trainingannotations.gff traininggenome.fa 5000 genes.gb
    numLocus=`grep -c LOCUS genes.gb`
    echo "$numLocus"
    numTest=$(( numLocus - 5500 ))
    echo "$numTest"
    #make training set
    randomSplit.pl genes.gb $numTest
    
    grep -c LOCUS *
    
    #train for new species
    new_species.pl --species=$species
    etraining --species=$species genes.gb.train
    
    #split genome
    awk '/>/{n++}{print >"SplitGenome"  n ".fa" }' "genome.fa"
    
    #make commands
    for query in ./SplitGenome*.fa; do
        numQueries=$(( $numQueries + 1 ))
        back=${query%.fa}
        count=${back:2}
        
        echo -e "augustus --species=${species} ${query} > annot_${count}.gtf" >> commands.txt
    done
    
    #run parallel
    parallel --jobs ${cpu} --memfree 10G --retries 500 < commands.txt
    
    numQueries=0
    for query in ./annot_*; do
        numQueries=$(( $numQueries + 1 ))
        cat annot_SplitGenome${numQueries}.gtf >> reference_augustus.gtf
        rm annot_SplitGenome${numQueries}.gtf
    done
    rm traininggenome.fa
    rm trainingannotations.gff
    
    awk '/# start gene/{n++}{print >"annotSplit"  n ".gtf" }' "reference_augustus.gtf"
    
    rm reference_augustus.gtf
    count=0
    for FILE in ./annotSplit*.gtf; do
        count=$(( $count + 1 ))
    done

    touch reference_augustus.gtf
    for (( i=1; i<=$count; i++ )); do
        file="annotSplit${i}.gtf"
        currentName=`head -n 1 ${file} | awk '{ print $4 }'`
        sed -i "s/gene ${currentName}/gene g${i}/g" $file
        sed -i "s/\t${currentName}/\tg${i}/g" $file
        sed -i "s/\"${currentName}/\"g${i}/g" $file
        cat $file >> reference_augustus.gtf
        rm $file
    done
    
    gffread reference_augustus.gtf -y reference_augustus.aa -g genome.fa

else
    cp $transcript transcripts.fasta
    #cat $rtranscript >> transcripts.fasta

    cp genome.fa /data/genome.fa
    cp transcripts.fasta /data/transcripts.fasta

    autoAugpasa.pl --genome=genome.fa --species=$species --cdna=transcripts.fasta --cpus=$cpu --verbose=2 --pasa --workingdir=/data/

    #/data/autoAug/autoAugPred_abinitio/shells/aug1
    #rm /data/autoAug/autoAugPred_abinitio/shells/shellForAug

    #for file in /data/autoAug/autoAugPred_abinitio/shells/*
    #do
    #  "$file"
    #done
    touch commands.txt
    for FILE in /data/autoAug/autoAugPred_abinitio/shells/aug*; do
        cat $FILE >> commands.txt
        echo -e "" >> commands.txt
    done

    parallel < commands.txt

    autoAug.pl --species=$species --genome=/data/autoAug/seq/genome_clean.fa --useexisting --hints=/data/autoAug/hints/hints.E.gff  -v -v --pasa --index=1 --workingdir=/data/

    touch commands2.txt
    for FILE in /data/autoAug/autoAugPred_hints/shells/aug*; do
        cat $FILE >> commands2.txt
        echo -e "" >> commands2.txt
    done

    sed 'sK../../hints/hints.gffK/data/autoAug/hints/hints.gffKg' commands2.txt > commands3.txt

    parallel < commands3.txt

    cat reference_exonerate.out >> exoneratecombined.out
    cat exonerate.out >> exoneratecombined.out

    exonerate2hints.pl --in=exoneratecombined.out --source=P --out=exonerate.hints

    cat exonerate.hints >> /data/autoAug/hints/hints.E.gff

    autoAug.pl --species=$species --genome=/data/autoAug/seq/genome_clean.fa --useexisting --hints=/data/autoAug/hints/hints.E.gff --estali=/data/autoAug/cdna/cdna.f.psl -v -v -v --pasa --index=2 --workingdir=/data/

    cp /data/autoAug/autoAugPred_abinitio/predictions/* .
    
    mv /data/autoAug/hints/hints.E.gff .

    #convert gtf file to gff3 for evm
    #gtf2gff.pl < augustus.gtf --out=augustus.gff3 --gff3

    #if [ "$protein" == "referenceProtein.fa" ]; then
    #    mv augustus.gff3 reference_augustus.gff3
    #    sed -i -e 's/AUGUSTUS/HOMOLOGY_AUGUSTUS/g' reference_augustus.gff3
    #fi

    #if [ "$transcript" == "reference_trinity.fasta" ] || [ "$transcript" == "formatted_referenceRNA.fa" ]; then
    #    mv augustus.gtf reference_augustus.gtf
    #    sed -i -e 's/AUGUSTUS/HOMOLOGY_AUGUSTUS/g' reference_augustus.gtf
    #fi
fi
