#!/bin/bash
#evidencemodeler.sh

usage() {
  echo "-h Help documentation for liftoff.sh"
  echo "-f  --filter out (Duplicated or Missing)"
  echo "Example: bash evidencemodeler.sh -i genome.fa -l -t "
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :f:h opt
do
    case $opt in
        f) filter=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${filter} ]]
then
    filter="Duplicated"
fi

grep -v '^#' full_table.tsv > full_tableEdited1.tsv
if [[ "$filter" == "Duplicated" ]]; then
    grep 'Duplicated' full_tableEdited1.tsv > full_tableEdited2.tsv
    
    oldHmm=""
    currentHmm=""
    dominantProgram=""
    currentProgram=""
    while IFS= read line; do
        oldHmm="$currentHmm"
        currentHmm=`echo "${line}" | awk '{ print $1 }' | awk -Fg '{print $NF}'`
        currentBusco=`echo "${line}" | awk '{ print $3 }' | awk -Fg '{print $NF}'`
        currentBusco="${currentBusco}"
        echo "new"
        echo "$oldHmm"
        echo "$currentHmm"
        echo "$currentBusco"
        if [[ "$currentHmm" == "$oldHmm" ]]; then
            mRNALineLift=`grep "${currentBusco}" completeBuscos.gff3 | grep -P "Liftoff\tmRNA"`
            geneLineAug=`grep "\-g${currentBusco::-3};Name" completeBuscos.gff3 | grep -P "Augustus\tgene"`
            mRNALineEVM=`grep "${currentBusco};" completeBuscos.gff3 | grep -P "EVM\tmRNA"`
            echo "signs"
            echo "$mRNALineLift"
            echo "$geneLineAug"
            if [[ "$geneLineAug" != "" ]]; then
                echo "Aug"
                grep -v -P "\-g${currentBusco::-3};Name" completeBuscos.gff3 > tmpfile && mv tmpfile completeBuscos.gff3
                grep -v -P "\-g${currentBusco}$" completeBuscos.gff3 > tmpfile && mv tmpfile completeBuscos.gff3
            elif [[ "$mRNALineLift" != "" ]]; then
                echo "liftoff"
                if [[ "$dominantProgram" != "Liftoff" ]]; then
                    first=`echo "$mRNALineLift" | cut -d$'\t' -f 4`
                    second=`echo "$mRNALineLift" | cut -d$'\t' -f 5`
                    liftoffGene=`grep -P "Liftoff\tgene\t${first}\t${second}" completeBuscos.gff3`
                    if [[ "$liftoffGene" != "" ]]; then
                        grep -v -P "Liftoff\tgene\t${first}\t${second}" completeBuscos.gff3 > tmpfile && mv tmpfile completeBuscos.gff3
                    else
                        tmp="${mRNALineLift##*Parent=}"
                        tmp2="${tmp%%;*}"
                        geneLine=`grep "${tmp2}" completeBuscos.gff3 | grep -P "Liftoff\ttranscript"`
                        echo "$tmp"
                        echo "$tmp2"
                        first=`echo "$geneLine" | cut -d$'\t' -f 4`
                        second=`echo "$geneLine" | cut -d$'\t' -f 5`
                        liftoffGene=`grep -P "Liftoff\tgene\t${first}\t${second}" completeBuscos.gff3`
                        mRNALineLift=`grep "${tmp2}" completeBuscos.gff3 | grep -P "Liftoff\ttranscript"`
                        first=`echo "$mRNALineLift" | cut -d$'\t' -f 4`
                        second=`echo "$mRNALineLift" | cut -d$'\t' -f 5`
                        grep -v -P "Liftoff\tgene\t${first}\t${second}" completeBuscos.gff3 > tmpfile && mv tmpfile completeBuscos.gff3
                        tmp="${geneLine##*ID=}"
                        tmp2="${tmp%%;*}"
                        grep -v "${tmp2}" completeBuscos.gff3 > tmpfile && mv tmpfile completeBuscos.gff3
                    fi
                fi
                grep -v "${currentBusco}" completeBuscos.gff3 > tmpfile && mv tmpfile completeBuscos.gff3
            elif [[ "$mRNALineEVM" != "" ]]; then
                echo "EVM"
                grep -v "ID=evm.TU${currentBusco:9};" completeBuscos.gff3 > tmpfile && mv tmpfile completeBuscos.gff3
                grep -v "Parent=evm.TU${currentBusco:9};" completeBuscos.gff3 > tmpfile && mv tmpfile completeBuscos.gff3
                grep -v "Parent=${currentBusco}$" completeBuscos.gff3 > tmpfile && mv tmpfile completeBuscos.gff3
                echo "evm.TU${currentBusco:9}"
            else
                echo "helix"
                grep -v -P "${currentBusco::-2}" completeBuscos.gff3 > tmpfile && mv tmpfile completeBuscos.gff3
            fi
        else
            echo "NEW"
            mRNALineLift=`grep "${currentBusco}" completeBuscos.gff3 | grep -P "Liftoff\tmRNA"`
            geneLineAug=`grep "\-g${currentBusco::-3};Name" completeBuscos.gff3 | grep -P "Augustus\tgene"`
            echo "testing"
            echo "old dominant ${dominantProgram}"
            if [[ "$geneLineAug" == "" ]] && [[ "$mRNALineLift" != "" ]]; then
                echo "dominant Lift"
                dominantProgram="Liftoff"
            else
                echo "dominant other"
                dominantProgram="Other"
            fi
        fi
    done < "full_tableEdited2.tsv"
else
    grep -v 'Missing' full_tableEdited1.tsv > full_tableEdited2.tsv
        
    while IFS= read line; do
        currentBusco=`echo "${line}" | awk '{ print $3 }' | awk -Fg '{print $NF}'`
        currentBusco="${currentBusco}"
        echo "new"

        echo "$currentBusco"
        mRNALineLift=`grep "${currentBusco::-2}" completeBuscos.gff3 | grep -P "Liftoff\tmRNA"`
        geneLineAug=`grep "\-g${currentBusco::-3};Name" completeBuscos.gff3 | grep -P "Augustus\tgene"`
        mRNALineEVM=`grep "${currentBusco};" completeBuscos.gff3 | grep -P "EVM\tmRNA"`
        echo "signs"
        echo "$mRNALineLift"
        echo "$geneLineAug"
        if [[ "$geneLineAug" != "" ]]; then
            echo "Aug"
            echo  -e "$geneLineAug" >> filteredCompleteBuscos.gff3
            grep -P "\-g${currentBusco}" completeBuscos.gff3 >> filteredCompleteBuscos.gff3
        elif [[ "$mRNALineLift" != "" ]]; then
            echo "liftoff"
            first=`echo "$mRNALineLift" | cut -d$'\t' -f 4`
            second=`echo "$mRNALineLift" | cut -d$'\t' -f 5`
            liftoffGene=`grep -P "Liftoff\tgene\t${first}\t${second}" completeBuscos.gff3`
            if [[ "$liftoffGene" != "" ]]; then
                echo "good"
                echo "$liftoffGene"
                grep -P "Liftoff\tgene\t${first}\t${second}" completeBuscos.gff3 >> filteredCompleteBuscos.gff3
            else
                echo "bad"
                tmp="${mRNALineLift##*Parent=}"
                tmp2"${tmp%%;*}"
                echo "$tmp"
                echo "$tmp2"
                geneLine=`grep "${tmp2}" completeBuscos.gff3 | grep -P "Liftoff\ttranscript"`
                first=`echo "$geneLine" | cut -d$'\t' -f 4`
                second=`echo "$geneLine" | cut -d$'\t' -f 5`
                liftoffGene=`grep -P "Liftoff\tgene\t${first}\t${second}" completeBuscos.gff3`
                tmp="${geneLine##*ID=}"
                tmp2="${tmp%%;*}"
                grep -P "Liftoff\tgene\t${first}\t${second}" completeBuscos.gff3 >> filteredCompleteBuscos.gff3
                grep "=${tmp2};"
            fi
            grep "${currentBusco}" completeBuscos.gff3 >> filteredCompleteBuscos.gff3
        elif [[ "$mRNALineEVM" != "" ]]; then
            echo "EVM"
            grep "ID=evm.TU${currentBusco:9};" completeBuscos.gff3 >> filteredCompleteBuscos.gff3
            grep "Parent=evm.TU${currentBusco:9};" completeBuscos.gff3 >> filteredCompleteBuscos.gff3
            grep "Parent=${currentBusco}$" completeBuscos.gff3 >> filteredCompleteBuscos.gff3
            echo "evm.TU${currentBusco:9}"
        else
            echo "helix"
            grep -P "${currentBusco::-2}" completeBuscos.gff3 >> filteredCompleteBuscos.gff3
        fi
    done < "full_tableEdited2.tsv"
fi


