#!/bin/bash
#evidencemodeler.sh

usage() {
  echo "-h Help documentation for liftoff.sh"
  echo "-i  --input genome file in fasta format"
  echo "-s  --cutoff for single exon genes percent certainty"
  echo "-m  --cutoff for multi exon genes percent certainty"
  echo "-l  --max gene length"
  echo "-e  --minimum percent of hint evidence needed to be considered a good gene"
  echo "-o  --only complete buscos"
  echo "-p  --program (Augustus or other)"
  echo "Example: bash evidencemodeler.sh -i genome.fa -l -t "
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:s:m:e:l:o:p:h opt
do
    case $opt in
        i) input=$OPTARG;;
        s) minPercentSingle=$OPTARG;;
        m) minPercentMulti=$OPTARG;;
        e) minPercentEvidence=$OPTARG;;
        l) maxGeneLength=$OPTARG;;
        o) onlycomplete=$OPTARG;;
        p) predictionProgram=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${input} ]]
then
    usage
fi
if [[ -z ${minPercentSingle} ]]
then
    minPercentSingle="99"
fi
if [[ -z ${minPercentMulti} ]]
then
    minPercentMulti="96"
fi
if [[ -z ${minPercentEvidence} ]]
then
    minPercentEvidence="90"
fi
if [[ -z ${maxGeneLength} ]]
then
    maxGeneLength="2300000"
fi
if [[ -z ${onlycomplete} ]]
then
    onlycomplete="False"
fi
if [[ -z ${predictionProgram} ]]
then
    predictionProgram="Augustus"
fi


if [[ "$predictionProgram" == "Augustus" ]]; then

    underTen=0
    underTwenty=0
    underThirty=0
    underForty=0
    underFifty=0
    underSixty=0
    underSeventy=0
    underEighty=0
    underNinty=0
    underOnehundred=0

    underTenSingle=0
    underTwentySingle=0
    underThirtySingle=0
    underFortySingle=0
    underFiftySingle=0
    underSixtySingle=0
    underSeventySingle=0
    underEightySingle=0
    underNintySingle=0
    underOnehundredSingle=0

    removeSingle=0
    removeMulti=0

    grep -v '^#' full_table.tsv > full_tableEdited1.tsv
    grep -v 'Missing' full_tableEdited1.tsv > full_tableEdited2.tsv
    if [[ "$onlycomplete" == "True" ]]; then
        mv full_tableEdited2.tsv tmpfull_tableEdited2.tsv
        grep -v 'Fragmented' tmpfull_tableEdited2.tsv > full_tableEdited2.tsv
    fi

    touch full_tableEdited3.tsv
    while IFS= read line; do
            currentBusco=`echo "${line}" | awk '{ print $3 }' | awk -F. '{print $1}' | awk -Fg '{print $NF}'`
            echo "$currentBusco" >> full_tableEdited3.tsv
    done < "full_tableEdited2.tsv"

    sort -h full_tableEdited3.tsv > full_tableEdited3Sorted.tsv

    awk '/# start gene/{n++}{print >"annotSplit"  n ".gtf" }' "$input"
    rm $input

    count=0
    for FILE in ./annotSplit*.gtf; do
        count=$(( $count + 1 ))
    done

    for (( i=1; i<=$count; i++ )); do
        file="annotSplit${i}.gtf"
        lineNum=0
        percentSupport=0
        numExons=0
        percentSupportFinal=0
        percentSupportFinal2=0
        numExonsFinal=0
        percentSupportFinalInt=0
        scoretmp=0
        score=0
        lowerLimit=0
        upperLimit=0
        geneLenth=0

        while IFS= read line; do
            if [[ $line == *"% of transcript supported by hints"* ]]; then
                percentSupport="${line: -3}"
                percentSupportFinal=${percentSupport#*/}
                percentSupportFinal2=${percentSupport#* }
                percentSupportFinalInt=${percentSupportFinal2%.*}
                #echo "$line"
            elif [[ $line == *"CDS exons:"* ]]; then
                numExons="${line: -3}"
                numExonsFinal=${numExons#*/}
                #echo "$line"
                break
            elif [[ $line == *"AUGUSTUS    gene"* ]]; then
                #echo "$line"
                wordNum=0
                for word in $line; do
                    wordNum=$(( $wordNum + 1 ))
                    if [ "$wordNum" == "6" ]; then
                        #echo "$word"
                        if [ "$word" == "1" ]; then
                            score="100"
                        else
                            scoretmp=${word#*.}
                            if [ "${#scoretmp}" == "1" ]; then
                                score="${scoretmp}0"
                            else
                                score=${scoretmp#*0}
                            fi
                        fi
                        #echo "$score"
                    elif [ "$wordNum" == "5" ]; then
                        upperLimit="$word"
                        geneLenth=$(( $upperLimit - $lowerLimit ))
                        if [ $((geneLenth)) -ge 2300000 ]; then
                            echo "$geneLenth"
                            echo "$line"
                        fi
                    elif [ "$wordNum" == "4" ]; then
                        lowerLimit="$word"
                    fi
                done
            fi
        done < "$file"
        

        #echo "new"
        #echo "$percentSupportFinalInt"
        #echo "$numExonsFinal"
        currentBusco=`head -n 1 full_tableEdited3Sorted.tsv | awk '{ print $1 }'`
        if [ $((i)) -eq $((currentBusco)) ]; then
            cat $file >> filteredAnnots.gtf
            sed -i '1d' full_tableEdited3Sorted.tsv
            if  [ $((numExonsFinal)) -le 1 ]; then
                buscoSingle=$(( $buscoSingle + 1 ))
            else
                buscoMulti=$(( $buscoMulti + 1 ))
            fi
        elif [ $((percentSupportFinalInt)) -lt $((minPercentEvidence)) ] && [ $((numExonsFinal)) -le 1 ] && [ $((score)) -lt $((minPercentSingle)) ]; then
            #rm $file
            #echo "$file"
            cat $file >> removedAnnots.gtf
            removeSingle=$(( $removeSingle + 1 ))
        elif [ $((percentSupportFinalInt)) -lt $((minPercentEvidence)) ] && [ $((score)) -lt $((minPercentMulti)) ]; then
            cat $file >> removedAnnots.gtf
            removeMulti=$(( $removeMulti + 1 ))
        elif [ $((geneLenth)) -ge $((maxGeneLength)) ]; then
            cat $file >> removedAnnots.gtf
            if  [ $((numExonsFinal)) -le 1 ]; then
                removeSingle=$(( $removeSingle + 1 ))
            else
                removeMulti=$(( $removeMulti + 1 ))
            fi
        else
            cat $file >> filteredAnnots.gtf
        fi
        if [ $((score)) -le 10 ]; then
            underTen=$(( $underTen + 1 ))
            if  [ $((numExonsFinal)) -le 1 ]; then
                underTenSingle=$(( $underTenSingle + 1 ))
            fi
        elif [ $((score)) -le 20 ]; then
            underTwenty=$(( $underTwenty + 1 ))
            if  [ $((numExonsFinal)) -le 1 ]; then
                underTwentySingle=$(( $underTwentySingle + 1 ))
            fi
        elif [ $((score)) -le 30 ]; then
            underThirty=$(( $underThirty + 1 ))
            if  [ $((numExonsFinal)) -le 1 ]; then
                underThirtySingle=$(( $underThirtySingle + 1 ))
            fi
        elif [ $((score)) -le 40 ]; then
            underForty=$(( $underForty + 1 ))
            if  [ $((numExonsFinal)) -le 1 ]; then
                underFortySingle=$(( $underFortySingle + 1 ))
            fi
        elif [ $((score)) -le 50 ]; then
            underFifty=$(( $underFifty + 1 ))
            if  [ $((numExonsFinal)) -le 1 ]; then
                underFiftySingle=$(( $underFiftySingle + 1 ))
            fi
        elif [ $((score)) -le 60 ]; then
            underSixty=$(( $underSixty + 1 ))
            if  [ $((numExonsFinal)) -le 1 ]; then
                underSixtySingle=$(( $underSixtySingle + 1 ))
            fi
        elif [ $((score)) -le 70 ]; then
            underSeventy=$(( $underSeventy + 1 ))
            if  [ $((numExonsFinal)) -le 1 ]; then
                underSeventySingle=$(( $underSeventySingle + 1 ))
            fi
        elif [ $((score)) -le 80 ]; then
            underEighty=$(( $underEighty + 1 ))
            if  [ $((numExonsFinal)) -le 1 ]; then
                underEightySingle=$(( $underEightySingle + 1 ))
            fi
        elif [ $((score)) -le 90 ]; then
            underNinty=$(( $underNinty + 1 ))
            if  [ $((numExonsFinal)) -le 1 ]; then
                underNintySingle=$(( $underNintySingle + 1 ))
            fi
        elif [ $((score)) -le 100 ]; then
            underOnehundred=$(( $underOnehundred + 1 ))
            if  [ $((numExonsFinal)) -le 1 ]; then
                underOnehundredSingle=$(( $underOnehundredSingle + 1 ))
            fi
        fi
        rm $file
    done

    total=$(( $underOnehundred + $underNinty + $underEighty + $underSeventy + $underSixty + $underFifty + $underForty + $underThirty + $underTwenty + $underTen ))
    totalSingle=$(( $underOnehundredSingle + $underNintySingle + $underEightySingle + $underSeventySingle + $underSixtySingle + $underFiftySingle + $underFortySingle + $underThirtySingle + $underTwentySingle + $underTenSingle ))
    totalBusco=$(( $buscoMulti + $buscoSingle ))
    underTenMulti=$(( $underTen - $underTenSingle ))
    underTwentyMulti=$(( $underTwenty - $underTwentySingle ))
    underThirtyMulti=$(( $underThirty - $underThirtySingle ))
    underFortyMulti=$(( $underForty - $underFortySingle ))
    underFiftyMulti=$(( $underFifty - $underFiftySingle ))
    underSixtyMulti=$(( $underSixty - $underSixtySingle ))
    underSeventyMulti=$(( $underSeventy - $underSeventySingle ))
    underEightyMulti=$(( $underEighty - $underEightySingle ))
    underNintyMulti=$(( $underNinty - $underNintySingle ))
    underOnehundredMulti=$(( $underOnehundred - $underOnehundredSingle ))
    totalMulti=$(( $total - $totalSingle ))
    removeTotal=$(( $removeMulti + $removeSingle ))
    buscoTotal=$(( $buscoMulti + $buscoSingle ))

    touch stats.txt
    echo -e "Percent\tAll\tMulti Exon\tSingle Exon" > stats.txt
    echo -e "10%:\t${underTen}\t${underTenMulti}\t${underTenSingle}" >> stats.txt
    echo -e "20%:\t${underTwenty}\t${underTwentyMulti}\t${underTwentySingle}" >> stats.txt
    echo -e "30%:\t${underThirty}\t${underThirtyMulti}\t${underThirtySingle}" >> stats.txt
    echo -e "40%:\t${underForty}\t${underFortyMulti}\t${underFortySingle}" >> stats.txt
    echo -e "50%:\t${underFifty}\t${underFiftyMulti}\t${underFiftySingle}" >> stats.txt
    echo -e "60%:\t${underSixty}\t${underSixtyMulti}\t${underSixtySingle}" >> stats.txt
    echo -e "70%:\t${underSeventy}\t${underSeventyMulti}\t${underSeventySingle}" >> stats.txt
    echo -e "80%:\t${underEighty}\t${underEightyMulti}\t${underEightySingle}" >> stats.txt
    echo -e "90%:\t${underNinty}\t${underNintyMulti}\t${underNintySingle}" >> stats.txt
    echo -e "100%:\t${underOnehundred}\t${underOnehundredMulti}\t${underOnehundredSingle}" >> stats.txt
    echo -e "Total Genes Found:\t${total}\t${totalMulti}\t${totalSingle}" >> stats.txt
    echo -e "Removed:\t${removeTotal}\t${removeMulti}\t${removeSingle}" >> stats.txt
    echo -e "Buscos Found:\t${buscoTotal}\t${buscoMulti}\t${buscoSingle}" >> stats.txt
    mv stats.txt ${input}.stats
elif [[ "$predictionProgram" == "Augustusgff3" ]]; then
    cp $input AugustusRemoved.gff3
    
    grep -v '^#' full_table.tsv > full_tableEdited1.tsv
    grep -v 'Missing' full_tableEdited1.tsv > full_tableEdited2.tsv
    if [[ "$onlycomplete" == "True" ]]; then
        mv full_tableEdited2.tsv tmpfull_tableEdited2.tsv
        grep -v 'Fragmented' tmpfull_tableEdited2.tsv > full_tableEdited2.tsv
    fi

    touch full_tableEdited3.tsv
    while IFS= read line; do
            currentBusco=`echo "${line}" | awk '{ print $3 }' | awk -Fg '{print $NF}'`
            echo "$currentBusco" >> full_tableEdited3.tsv
    done < "full_tableEdited2.tsv"
    
    buscoFound="False"
    while IFS= read line; do
        currentBusco="${line}"
        echo "$currentBusco"
        geneLineAug=`grep "\-g${currentBusco::-3};Name" AugustusRemoved.gff3 | grep -P "Augustus\tgene"`
        echo "$geneLineAug"
        echo "Aug"
        grep -P "\-g${currentBusco::-3};Name" AugustusRemoved.gff3 >> AugBuscos.gff3
        grep -P "\-g${currentBusco}$" AugustusRemoved.gff3 >> AugBuscos.gff3
        grep -v -P "\-g${currentBusco::-3};Name" AugustusRemoved.gff3 > tmpfile && mv tmpfile AugustusRemoved.gff3
        grep -v -P "\-g${currentBusco}$" AugustusRemoved.gff3 > tmpfile && mv tmpfile AugustusRemoved.gff3
    done < "full_tableEdited3.tsv"
    
elif [[ "$predictionProgram" == "Liftoff" ]]; then
    cp $input LiftoffRemoved.gff3
    
    grep -v '^#' full_table.tsv > full_tableEdited1.tsv
    grep -v 'Missing' full_tableEdited1.tsv > full_tableEdited2.tsv
    if [[ "$onlycomplete" == "True" ]]; then
        mv full_tableEdited2.tsv tmpfull_tableEdited2.tsv
        grep -v 'Fragmented' tmpfull_tableEdited2.tsv > full_tableEdited2.tsv
    fi

    touch full_tableEdited3.tsv
    while IFS= read line; do
            currentBusco=`echo "${line}" | awk '{ print $3 }' | awk -Fg '{print $NF}'`
            echo "$currentBusco" >> full_tableEdited3.tsv
    done < "full_tableEdited2.tsv"
    
    buscoFound="False"
    while IFS= read line; do
        currentBusco="${line}"
        echo "$currentBusco"
        mRNALine=`grep "${currentBusco}" LiftoffRemoved.gff3 | grep -P "Liftoff\tmRNA"`
        first=`echo "$mRNALine" | cut -d$'\t' -f 4`
        second=`echo "$mRNALine" | cut -d$'\t' -f 5`
        #geneLine="#duplicate"
        geneLine=`grep -P "Liftoff\tgene\t${first}\t${second}" LiftoffRemoved.gff3`
        echo "$geneLine"
        if [[ "$geneLine" != "" ]]; then
           echo "$geneLine" >> LiftoffBuscos.gff3
           grep -v -P "Liftoff\tgene\t${first}\t${second}" LiftoffRemoved.gff3 > tmpfile && mv tmpfile LiftoffRemoved.gff3
        elif [[ "$mRNALine" != "" ]]; then
            tmp="${mRNALine##*Parent=}"
            tmp2="${tmp%%;*}"
            echo "$tmp"
            echo "$tmp2"
            echo "liftoff gene is empty"
            #geneLine=`grep "${tmp2}" liftoff.gff3 | grep -P "Liftoff\tgene"`
            #echo "$geneLine" >> LiftoffBuscos.gff3
            #first=`echo "$geneLine" | cut -d$'\t' -f 4`
            #second=`echo "$geneLine" | cut -d$'\t' -f 5`
            #grep -v -P "Liftoff\tgene\t${first}\t${second}" LiftoffRemoved.gff3 > tmpfile && mv tmpfile LiftoffRemoved.gff3
            first=`echo "$mRNALine" | cut -d$'\t' -f 4`
            second=`echo "$mRNALine" | cut -d$'\t' -f 5`
            chrName=`echo "$mRNALine" | cut -d$'\t' -f 1`
            frStrand=`echo "$mRNALine" | cut -d$'\t' -f 7`
            newGeneLine="${chrName}\tLiftoff\tgene\t${first}\t${second}\t.\t${frStrand}\t.\tID=${tmp2};Name=${tmp2##*-};"
            testLiftLine=`grep -P "\tID=${tmp2};" LiftoffBuscos.gff3 | grep -P "${chrName}\tLiftoff\tgene\t"`
            if [[ "$testLiftLine" != "" ]]; then
                continue
            else
                firstTestLift=`echo "$testLiftLine" | cut -d$'\t' -f 4`
                secondTestLift=`echo "$testLiftLine" | cut -d$'\t' -f 5`
                firstTestLiftInt=$(( $firstTestLift + 0 ))
                secondTestLiftInt=$(( $secondTestLiftLift + 0 ))
                firstInt=$(( $first + 0 ))
                secondInt=$(( $second + 0 ))
                echo "$firstTestLiftInt vs $firstInt"
                echo "$secondTestLiftInt vs $secondInt"
                echo "${testLiftLine}"
                echo "$tmp2"
                echo "$newGeneLine"
                if [[ $firstTestLift -gt 0 ]] && [[ $secondTestLiftInt -gt 0 ]] && [[ $firstTestLift -lt $firstInt || $secondTestLiftInt -gt $secondInt ]]; then
                    grep -v -P "${newGeneLine}" LiftoffBuscos.gff3 > tmpfile && mv tmpfile LiftoffBuscos.gff3
                else
                    echo "this will make a brand new gene line"
                    echo -e "${newGeneLine}" >> LiftoffBuscos.gff3
                fi
            fi
        fi
        grep "=${currentBusco};" LiftoffRemoved.gff3 >> LiftoffBuscos.gff3
        grep -v "=${currentBusco};" LiftoffRemoved.gff3 > tmpfile && mv tmpfile LiftoffRemoved.gff3
    done < "full_tableEdited3.tsv"
    
elif [[ "$predictionProgram" == "Helixer" ]]; then
    cp $input HelixerRemoved.gff3

    grep -v '^#' full_table.tsv > full_tableEdited1.tsv
    grep -v 'Missing' full_tableEdited1.tsv > full_tableEdited2.tsv
    if [[ "$onlycomplete" == "True" ]]; then
        mv full_tableEdited2.tsv tmpfull_tableEdited2.tsv
        grep -v 'Fragmented' tmpfull_tableEdited2.tsv > full_tableEdited2.tsv
    fi

    touch full_tableEdited3.tsv
    while IFS= read line; do
            currentBusco=`echo "${line}" | awk '{ print $3 }'`
            echo "$currentBusco" >> full_tableEdited3.tsv
    done < "full_tableEdited2.tsv"

    sort -h full_tableEdited3.tsv > full_tableEdited3Sorted.tsv
    
    while IFS= read line; do
        currentBusco="${line}"
        echo "$currentBusco"
        #grep "${currentBusco}" annot.gff3 >> testfile.gff3
        #grep -v -P "${currentBusco}" annot.gff3 > tmpfile && mv tmpfile annot.gff3
        grep "ID=${currentBusco::-2}$" HelixerRemoved.gff3 >> HelixerBuscos.gff3
        grep "Parent=${currentBusco::-2}$" HelixerRemoved.gff3 >> HelixerBuscos.gff3
        grep "Parent=${currentBusco}$" HelixerRemoved.gff3 >> HelixerBuscos.gff3
        grep -v "ID=${currentBusco::-2}$" HelixerRemoved.gff3 > tmpfile && mv tmpfile HelixerRemoved.gff3
        grep -v "Parent=${currentBusco::-2}$" HelixerRemoved.gff3 > tmpfile && mv tmpfile HelixerRemoved.gff3
        grep -v "Parent=${currentBusco}$" HelixerRemoved.gff3 > tmpfile && mv tmpfile HelixerRemoved.gff3
    done < "full_tableEdited3Sorted.tsv"
    
elif [[ "$predictionProgram" == "EVM" ]]; then
    cp $input EVMRemoved.gff3

    grep -v '^#' full_table.tsv > full_tableEdited1.tsv
    grep -v 'Missing' full_tableEdited1.tsv > full_tableEdited2.tsv
    if [[ "$onlycomplete" == "True" ]]; then
        mv full_tableEdited2.tsv tmpfull_tableEdited2.tsv
        grep -v 'Fragmented' tmpfull_tableEdited2.tsv > full_tableEdited2.tsv
    fi

    touch full_tableEdited3.tsv
    while IFS= read line; do
            currentBusco=`echo "${line}" | awk '{ print $3 }' | awk -Fg '{print $NF}'`
            echo "${currentBusco}" >> full_tableEdited3.tsv
    done < "full_tableEdited2.tsv"

    sort -h full_tableEdited3.tsv > full_tableEdited3Sorted.tsv
    
    while IFS= read line; do
        currentBusco="${line}"
        echo "$currentBusco"
        echo "yes"
        grep "ID=evm.TU${currentBusco:9};" EVMRemoved.gff3 >> EVMBuscos.gff3
        testvar=`grep "ID=evm.TU${currentBusco:9};" EVMRemoved.gff3`
        echo "$testvar"
        echo `grep "Parent=evm.TU${currentBusco:9};" EVMRemoved.gff3`
        grep "Parent=evm.TU${currentBusco:9};" EVMRemoved.gff3 >> EVMBuscos.gff3
        grep "Parent=${currentBusco}$" EVMRemoved.gff3 >> EVMBuscos.gff3
        grep -v "ID=evm.TU${currentBusco:9};" EVMRemoved.gff3 > tmpfile && mv tmpfile EVMRemoved.gff3
        grep -v "Parent=evm.TU${currentBusco:9};" EVMRemoved.gff3 > tmpfile && mv tmpfile EVMRemoved.gff3
        grep -v "Parent=${currentBusco}$" EVMRemoved.gff3 > tmpfile && mv tmpfile EVMRemoved.gff3
    done < "full_tableEdited3Sorted.tsv"

elif [[ "$predictionProgram" == "OtherAug" ]]; then

    grep -v '^#' full_table.tsv > full_tableEdited1.tsv
    grep -v 'Missing' full_tableEdited1.tsv > full_tableEdited2.tsv
    if [[ "$onlycomplete" == "True" ]]; then
        mv full_tableEdited2.tsv tmpfull_tableEdited2.tsv
        grep -v 'Fragmented' tmpfull_tableEdited2.tsv > full_tableEdited2.tsv
    fi

    touch full_tableEdited3.tsv
    while IFS= read line; do
            currentBusco=`echo "${line}" | awk '{ print $3 }' | awk -F. '{print $1}' | awk -Fg '{print $NF}'`
            echo "$currentBusco" >> full_tableEdited3.tsv
    done < "full_tableEdited2.tsv"

    sort -h full_tableEdited3.tsv > full_tableEdited3Sorted.tsv

    awk '/# start gene/{n++}{print >"annotSplit"  n ".gtf" }' "$input"
    
    cp $input LiftoffRemoved.gff3
    
    grep -v '^#' full_table.tsv > full_tableEdited1.tsv
    grep -v 'Missing' full_tableEdited1.tsv > full_tableEdited2.tsv
    if [[ "$onlycomplete" == "True" ]]; then
        mv full_tableEdited2.tsv tmpfull_tableEdited2.tsv
        grep -v 'Fragmented' tmpfull_tableEdited2.tsv > full_tableEdited2.tsv
    fi

    touch full_tableEdited3.tsv
    while IFS= read line; do
            currentBusco=`echo "${line}" | awk '{ print $3 }' | awk -Fg '{print $NF}'`
            echo "$currentBusco" >> full_tableEdited3.tsv
    done < "full_tableEdited2.tsv"
    sort -h full_tableEdited3.tsv > full_tableEdited3Sorted.tsv

    awk '/gene/{n++}{print >"annotSplit"  n ".gtf" }' "$input"
    rm $input

    count=0
    for FILE in ./annotSplit*.gtf; do
        count=$(( $count + 1 ))
    done

    for (( i=1; i<=$count; i++ )); do
        file="annotSplit${i}.gtf"
        buscoFound="False"
        while IFS= read line; do
            currentBusco="${line:1}"
            if grep -q "$currentBusco" "$file"; then
                buscoFound="True"
                break
            fi
        done < "full_tableEdited3Sorted.tsv"

        
        if [ "$buscoFound" == "True" ]; then
            cat $file >> filteredAnnots.gtf
            sed -i '1d' full_tableEdited3Sorted.tsv
            grep -v "currentBusco" full_tableEdited3Sorted.tsv > tmpfile
            mv tmpfile full_tableEdited3Sorted.tsv
        else
            cat $file >> removedAnnots.gtf
        fi
    done

else

    cp $input LiftoffRemoved.gff3
    
    grep -v '^#' full_table.tsv > full_tableEdited1.tsv
    grep -v 'Missing' full_tableEdited1.tsv > full_tableEdited2.tsv
    if [[ "$onlycomplete" == "True" ]]; then
        mv full_tableEdited2.tsv tmpfull_tableEdited2.tsv
        grep -v 'Fragmented' tmpfull_tableEdited2.tsv > full_tableEdited2.tsv
    fi

    touch full_tableEdited3.tsv
    while IFS= read line; do
            currentBusco=`echo "${line}" | awk '{ print $3 }' | awk -Fg '{print $NF}'`
            echo "$currentBusco" >> full_tableEdited3.tsv
    done < "full_tableEdited2.tsv"
    sort -h full_tableEdited3.tsv > full_tableEdited3Sorted.tsv

    awk '/gene/{n++}{print >"annotSplit"  n ".gtf" }' "$input"
    rm $input

    count=0
    for FILE in ./annotSplit*.gtf; do
        count=$(( $count + 1 ))
    done

    for (( i=1; i<=$count; i++ )); do
        file="annotSplit${i}.gtf"
        buscoFound="False"
        while IFS= read line; do
            currentBusco="${line:1}"
            if grep -q "$currentBusco" "$file"; then
                buscoFound="True"
                break
            fi
        done < "full_tableEdited3Sorted.tsv"

        
        if [ "$buscoFound" == "True" ]; then
            cat $file >> filteredAnnots.gtf
            sed -i '1d' full_tableEdited3Sorted.tsv
            grep -v "currentBusco" full_tableEdited3Sorted.tsv > tmpfile
            mv tmpfile full_tableEdited3Sorted.tsv
        else
            cat $file >> removedAnnots.gtf
        fi
    done
fi


rm full_tableEdited*.tsv

