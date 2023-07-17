#!/bin/bash
#augustus_full.sh

usage() {
  echo "-h Help documentation for tasser.sh"
  echo "-g  --input genome file"
  echo "-r  --input transcript file"
  echo "-s  --input species"
  echo "-e  --input exonerate.out file"
  echo "-p  --pretrained augustus (True or False)"
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
if [[ -z ${species} ]]
then
    species='current'
fi
if [[ -z ${pretrained} ]]
then
    pretrained='false'
fi
if [[ -z ${cpu} ]]
then
    cpu=`nproc`
fi

lineage="$species"

############################################################
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
unset __conda_setup
# <<< conda initialize <<<
############################################################

baseDir="`dirname \"$0\"`"

if [[ -z ${exoneratehints} ]]; then
    mv $genome genome.fa
    if [[ "$transcript" == "Helixer.gff3" ]]; then
        gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o tmp.gff3  Helixer.gff3
        mv tmp.gff3 Helixer.gff3
        agat_sp_extract_sequences.pl --clean_final_stop --gff Helixer.gff3 -f genome.fa -p -o protein.fa
        conda activate BUSCO
        busco -i protein.fa -l ${lineage} -o buscoout -m protein -c $cpu
        mv buscoout/run_${lineage}/full_table.tsv .
        conda deactivate
        bash ${baseDir}/annotationFilter.sh -i Helixer.gff3 -o True -p Helixer
        
        #make hints if they are present
        hints="False"
        if [[ -s pasa_pblat.pslx ]]; then
            blat2hints.pl --in=pasa_pblat.pslx --out=pasa_pblat.hints.gff --minintronlen=35 --trunkSS 1>blat2hints.stdout 2>blat2hints.stderr
            cat pasa_pblat.hints.gff >> hints.E.gff
            hints="True"
        fi
        if [[ -s reference_pasa_pblat.pslx ]]; then
            blat2hints.pl --in=reference_pasa_pblat.pslx --out=reference_pasa_pblat.hints.gff --minintronlen=35 --trunkSS 1>blat2hints.stdout 2>blat2hints.stderr
            cat reference_pasa_pblat.hints.gff >> hints.E.gff
            hints="True"
        fi
        if [[ -s reference_exonerate.out ]]; then
            cat reference_exonerate.out >> exoneratecombined.out
        fi
        if [[ -s exonerate.out ]]; then
            cat exonerate.out >> exoneratecombined.out
        fi
        if [[ -s exoneratecombined.out ]]; then
            exonerate2hints.pl --in=exoneratecombined.out --source=P --out=exonerate.hints
            cat exonerate.hints >> hints.E.gff
            hints="True"
        fi
        if [[ -s genomethreader.gff3 ]]; then
            perl /opt/GALBA/scripts/aln2hints.pl --in=genomethreader.gff3 --out=genomethreader.hints.gff --prg=gth
            cat genomethreader.hints.gff >> hints.E.gff
            hints="True"
        fi
        if [[ -s reference_genomethreader.gff3 ]]; then
            perl /opt/GALBA/scripts/aln2hints.pl --in=reference_genomethreader.gff3 --out=reference_genomethreader.hints.gff --prg=gth
            cat reference_genomethreader.hints.gff >> hints.E.gff
            hints="True"
        fi
        if [[ -s miniprot.gtf ]]; then
            perl /opt/GALBA/scripts/aln2hints.pl --in=miniprot.gtf --out=miniprot.hints.gff --prg=miniprot
            cat miniprot.hints.gff >> hints.E.gff
            hints="True"
        fi
        if [[ -s reference_miniprot.gff3 ]]; then
            perl /opt/GALBA/scripts/aln2hints.pl --in=reference_miniprot.gtf --out=reference_miniprot.hints.gff --prg=miniprot
            cat miniprot.hints.gff >> hints.E.gff
            hints="True"
        fi
        
        if [[ "$hints" == "True" ]]; then
            hintCommand="--hintsfile=hints.E.gff --extrinsicCfgFile=/root/augustus/config/extrinsic/extrinsic.M.RM.E.W.P.cfg "
        else
            hintCommand=""
        fi
        
        gff2gbSmallDNA.pl HelixerBuscos.gff3 genome.fa 5000 genes.gb

        #train for new species
        export AUGUSTUS_CONFIG_PATH=/root/augustus/config/
        new_species.pl --species=sample
        etraining --species=sample genes.gb
        
        #split genome
        awk '/>/{n++}{print >"SplitGenome"  n ".fa" }' "genome.fa"
        
        #make commands
        for query in ./SplitGenome*.fa; do
            numQueries=$(( $numQueries + 1 ))
            back=${query%.fa}
            count=${back:2}
            
            echo -e "augustus --species=sample ${hintCommand}${query} > annot_${count}.gtf" >> commands.txt
        done
        
        #run parallel
        parallel --jobs ${cpu} --memfree 10G --retries 500 < commands.txt
        
        numQueries=0
        for query in ./annot_*; do
            numQueries=$(( $numQueries + 1 ))
            cat annot_SplitGenome${numQueries}.gtf >> Helixer_augustus.gtf
            rm annot_SplitGenome${numQueries}.gtf
        done
        rm genome.fa
        rm trainingannotations.gff
        
        awk '/# start gene/{n++}{print >"annotSplit"  n ".gtf" }' "Helixer_augustus.gtf"
        
        rm Helixer_augustus.gtf
        count=0
        for FILE in ./annotSplit*.gtf; do
            count=$(( $count + 1 ))
        done

        touch Helixer_augustus.gtf
        for (( i=1; i<=$count; i++ )); do
            file="annotSplit${i}.gtf"
            currentName=`head -n 1 ${file} | awk '{ print $4 }'`
            sed -i "s/gene ${currentName}/gene g${i}/g" $file
            sed -i "s/\t${currentName}/\tg${i}/g" $file
            sed -i "s/\"${currentName}/\"g${i}/g" $file
            cat $file >> Helixer_augustus.gtf
            rm $file
        done
        rm annotSplit.gtf
        exit 0
    
    elif [[ "$transcript" = "liftoff.gff3" ]]; then
        gt  gff3  -force  -tidy  -sort  -retainids  -checkids  -o tmp.gff3  liftoff.gff3
        mv tmp.gff3 liftoff.gff3
        agat_sp_extract_sequences.pl --clean_final_stop --gff liftoff.gff3 -f genome.fa -p -o protein.fa
        conda activate BUSCO
        busco -i protein.fa -l ${lineage} -o buscoout -m protein -c $cpu
        mv buscoout/run_${lineage}/full_table.tsv .
        conda deactivate
        bash ${baseDir}/annotationFilter.sh -i liftoff.gff3 -o True -p Liftoff
        
        #make hints if they are present
        hints="False"
        if [[ -s pasa_pblat.pslx ]]; then
            blat2hints.pl --in=pasa_pblat.pslx --out=pasa_pblat.hints.gff --minintronlen=35 --trunkSS 1>blat2hints.stdout 2>blat2hints.stderr
            cat pasa_pblat.hints.gff >> hints.E.gff
            hints="True"
        fi
        if [[ -s reference_pasa_pblat.pslx ]]; then
            blat2hints.pl --in=reference_pasa_pblat.pslx --out=reference_pasa_pblat.hints.gff --minintronlen=35 --trunkSS 1>blat2hints.stdout 2>blat2hints.stderr
            cat reference_pasa_pblat.hints.gff >> hints.E.gff
            hints="True"
        fi
        if [[ -s reference_exonerate.out ]]; then
            cat reference_exonerate.out >> exoneratecombined.out
        fi
        if [[ -s exonerate.out ]]; then
            cat exonerate.out >> exoneratecombined.out
        fi
        if [[ -s exoneratecombined.out ]]; then
            exonerate2hints.pl --in=exoneratecombined.out --source=P --out=exonerate.hints
            cat exonerate.hints >> hints.E.gff
            hints="True"
        fi
        if [[ -s genomethreader.gff3 ]]; then
            perl /opt/GALBA/scripts/al2hints.pl --in=genomethreader.gff3 --out=genomethreader.hints.gff --prg=gth
            cat genomethreader.hints.gff >> hints.E.gff
            hints="True"
        fi
        if [[ -s reference_genomethreader.gff3 ]]; then
            perl /opt/GALBA/scripts/al2hints.pl --in=reference_genomethreader.gff3 --out=reference_genomethreader.hints.gff --prg=gth
            cat reference_genomethreader.hints.gff >> hints.E.gff
            hints="True"
        fi
        if [[ -s miniprot.gtf ]]; then
            perl /opt/GALBA/scripts/al2hints.pl --in=miniprot.gtf --out=miniprot.hints.gff --prg=miniprot
            cat miniprot.hints.gff >> hints.E.gff
            hints="True"
        fi
        if [[ -s reference_miniprot.gff3 ]]; then
            perl /opt/GALBA/scripts/al2hints.pl --in=reference_miniprot.gtf --out=reference_miniprot.hints.gff --prg=miniprot
            cat reference_miniprot.hints.gff >> hints.E.gff
            hints="True"
        fi
        
        if [[ "$hints" == "True" ]]; then
            hintCommand="--hintsfile=hints.E.gff --extrinsicCfgFile=/root/augustus/config/extrinsic/extrinsic.M.RM.E.W.P.cfg "
        else
            hintCommand=""
        fi
        
        gff2gbSmallDNA.pl LiftoffBuscos.gff3 genome.fa 5000 genes.gb

        #train for new species
        export AUGUSTUS_CONFIG_PATH=/root/augustus/config/
        new_species.pl --species=sample
        etraining --species=sample genes.gb
        
        #split genome
        awk '/>/{n++}{print >"SplitGenome"  n ".fa" }' "genome.fa"
        
        #make commands
        for query in ./SplitGenome*.fa; do
            numQueries=$(( $numQueries + 1 ))
            back=${query%.fa}
            count=${back:2}
            
            echo -e "augustus --species=sample ${hintCommand}${query} > annot_${count}.gtf" >> commands.txt
        done
        
        #run parallel
        parallel --jobs ${cpu} --memfree 10G --retries 500 < commands.txt
        
        numQueries=0
        for query in ./annot_*; do
            numQueries=$(( $numQueries + 1 ))
            cat annot_SplitGenome${numQueries}.gtf >> liftoff_augustus.gtf
            rm annot_SplitGenome${numQueries}.gtf
        done
        rm genome.fa
        rm trainingannotations.gff
        
        awk '/# start gene/{n++}{print >"annotSplit"  n ".gtf" }' "liftoff_augustus.gtf"
        
        rm liftoff_augustus.gtf
        count=0
        for FILE in ./annotSplit*.gtf; do
            count=$(( $count + 1 ))
        done

        touch liftoff_augustus.gtf
        for (( i=1; i<=$count; i++ )); do
            file="annotSplit${i}.gtf"
            currentName=`head -n 1 ${file} | awk '{ print $4 }'`
            sed -i "s/gene ${currentName}/gene g${i}/g" $file
            sed -i "s/\t${currentName}/\tg${i}/g" $file
            sed -i "s/\"${currentName}/\"g${i}/g" $file
            cat $file >> liftoff_augustus.gtf
            rm $file
        done
        rm annotSplit.gtf
        exit 0
        
    else
        usage
    fi

else
    mv $genome genome.fa
    #all old augustus run methods

    transcriptExtension="${transcript: -4}"
    echo "Last 3 characters of transcript: ${transcriptExtension}"
    if [[ "$pretrained" == "True" ]] || [[ "$pretrained" == "true" ]]; then
        #make hints if they are present
        hints="False"
        if [[ -s pasa_pblat.pslx ]]; then
            blat2hints.pl --in=pasa_pblat.pslx --out=pasa_pblat.hints.gff --minintronlen=35 --trunkSS 1>blat2hints.stdout 2>blat2hints.stderr
            cat pasa_pblat.hints.gff >> hints.E.gff
            hints="True"
        fi
        if [[ -s reference_pasa_pblat.pslx ]]; then
            blat2hints.pl --in=reference_pasa_pblat.pslx --out=reference_pasa_pblat.hints.gff --minintronlen=35 --trunkSS 1>blat2hints.stdout 2>blat2hints.stderr
            cat reference_pasa_pblat.hints.gff >> hints.E.gff
            hints="True"
        fi
        if [[ -s reference_exonerate.out ]]; then
            cat reference_exonerate.out >> exoneratecombined.out
        fi
        if [[ -s exonerate.out ]]; then
            cat exonerate.out >> exoneratecombined.out
        fi
        if [[ -s exoneratecombined.out ]]; then
            exonerate2hints.pl --in=exoneratecombined.out --source=P --out=exonerate.hints
            cat exonerate.hints >> hints.E.gff
            hints="True"
        fi
        if [[ -s genomethreader.gff3 ]]; then
            perl /opt/GALBA/scripts/al2hints.pl --in=genomethreader.gff3 --out=genomethreader.hints.gff --prg=gth
            cat genomethreader.hints.gff >> hints.E.gff
            hints="True"
        fi
        if [[ -s reference_genomethreader.gff3 ]]; then
            perl /opt/GALBA/scripts/al2hints.pl --in=reference_genomethreader.gff3 --out=reference_genomethreader.hints.gff --prg=gth
            cat reference_genomethreader.hints.gff >> hints.E.gff
            hints="True"
        fi
        if [[ -s miniprot.gtf ]]; then
            perl /opt/GALBA/scripts/al2hints.pl --in=miniprot.gtf --out=miniprot.hints.gff --prg=miniprot
            cat miniprot.hints.gff >> hints.E.gff
            hints="True"
        fi
        if [[ -s reference_miniprot.gff3 ]]; then
            perl /opt/GALBA/scripts/al2hints.pl --in=reference_miniprot.gtf --out=reference_miniprot.hints.gff --prg=miniprot
            cat miniprot.hints.gff >> hints.E.gff
            hints="True"
        fi
        
        if [[ "$hints" == "True" ]]; then
            hintCommand="--hintsfile=hints.E.gff --extrinsicCfgFile=/root/augustus/config/extrinsic/extrinsic.M.RM.E.W.P.cfg "
        else
            hintCommand=""
        fi
        
        #split genome
        awk '/>/{n++}{print >"SplitGenome"  n ".fa" }' "genome.fa"
        
        #make commands
        for query in ./SplitGenome*.fa; do
            numQueries=$(( $numQueries + 1 ))
            back=${query%.fa}
            count=${back:2}
            
            echo -e "augustus --species=${species} ${hintCommand}${query} > annot_${count}.gtf" >> commands.txt
        done
        
        #run parallel
        parallel --jobs ${cpu} --memfree 10G --retries 500 < commands.txt
        
        numQueries=0
        for query in ./annot_*; do
            numQueries=$(( $numQueries + 1 ))
            cat annot_SplitGenome${numQueries}.gtf >> pretrained_augustus.gtf
            rm annot_SplitGenome${numQueries}.gtf
        done
        
        awk '/# start gene/{n++}{print >"annotSplit"  n ".gtf" }' "pretrained_augustus.gtf"
        
        rm pretrained_augustus.gtf
        count=0
        for FILE in ./annotSplit*.gtf; do
            count=$(( $count + 1 ))
        done

        touch pretrained_augustus.gtf
        for (( i=1; i<=$count; i++ )); do
            file="annotSplit${i}.gtf"
            currentName=`head -n 1 ${file} | awk '{ print $4 }'`
            sed -i "s/gene ${currentName}/gene g${i}/g" $file
            sed -i "s/\t${currentName}/\tg${i}/g" $file
            sed -i "s/\"${currentName}/\"g${i}/g" $file
            cat $file >> pretrained_augustus.gtf
            rm $file
        done
        rm annotSplit.gtf
        exit 0
    elif [[ "${transcriptExtension}" == ".gff" ]] || [[ "${transcriptExtension}" == ".gtf" ]] || [[ "${transcriptExtension}" == "gff3" ]]; then

        #make hints if they are present
        hints="False"
        if [[ -s pasa_pblat.pslx ]]; then
            blat2hints.pl --in=pasa_pblat.pslx --out=pasa_pblat.hints.gff --minintronlen=35 --trunkSS 1>blat2hints.stdout 2>blat2hints.stderr
            cat pasa_pblat.hints.gff >> hints.E.gff
            hints="True"
        fi
        if [[ -s reference_pasa_pblat.pslx ]]; then
            blat2hints.pl --in=reference_pasa_pblat.pslx --out=reference_pasa_pblat.hints.gff --minintronlen=35 --trunkSS 1>blat2hints.stdout 2>blat2hints.stderr
            cat reference_pasa_pblat.hints.gff >> hints.E.gff
            hints="True"
        fi
        if [[ -s reference_exonerate.out ]]; then
            cat reference_exonerate.out >> exoneratecombined.out
        fi
        if [[ -s exonerate.out ]]; then
            cat exonerate.out >> exoneratecombined.out
        fi
        if [[ -s exoneratecombined.out ]]; then
            exonerate2hints.pl --in=exoneratecombined.out --source=P --out=exonerate.hints
            cat exonerate.hints >> hints.E.gff
            hints="True"
        fi
        if [[ -s genomethreader.gff3 ]]; then
            perl /opt/GALBA/scripts/al2hints.pl --in=genomethreader.gff3 --out=genomethreader.hints.gff --prg=gth
            cat genomethreader.hints.gff >> hints.E.gff
            hints="True"
        fi
        if [[ -s reference_genomethreader.gff3 ]]; then
            perl /opt/GALBA/scripts/al2hints.pl --in=reference_genomethreader.gff3 --out=reference_genomethreader.hints.gff --prg=gth
            cat reference_genomethreader.hints.gff >> hints.E.gff
            hints="True"
        fi
        if [[ -s miniprot.gtf ]]; then
            perl /opt/GALBA/scripts/al2hints.pl --in=miniprot.gtf --out=miniprot.hints.gff --prg=miniprot
            cat miniprot.hints.gff >> hints.E.gff
            hints="True"
        fi
        if [[ -s reference_miniprot.gff3 ]]; then
            perl /opt/GALBA/scripts/al2hints.pl --in=reference_miniprot.gtf --out=reference_miniprot.hints.gff --prg=miniprot
            cat miniprot.hints.gff >> hints.E.gff
            hints="True"
        fi
        
        if [[ "$hints" == "True" ]]; then
            hintCommand="--hintsfile=hints.E.gff --extrinsicCfgFile=/root/augustus/config/extrinsic/extrinsic.M.RM.E.W.P.cfg "
        else
            hintCommand=""
        fi
        
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
        export AUGUSTUS_CONFIG_PATH=/root/augustus/config/
        new_species.pl --species=$species
        etraining --species=$species genes.gb.train
        
        #split genome
        awk '/>/{n++}{print >"SplitGenome"  n ".fa" }' "genome.fa"
        
        #make commands
        for query in ./SplitGenome*.fa; do
            numQueries=$(( $numQueries + 1 ))
            back=${query%.fa}
            count=${back:2}
            
            echo -e "augustus --species=${species} ${hintCommand}${query} > annot_${count}.gtf" >> commands.txt
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
        
        #gffread reference_augustus.gtf -y reference_augustus.aa -g genome.fa
        rm annotSplit.gtf
        exit 0
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
fi
