#!/bin/bash
#evidencemodeler.sh

usage() {
  echo "-h Help documentation for liftoff.sh"
  echo "-a  --annotation gtf file named annotSplit*.gtf"
  echo "-s  --species name"
  echo "Example: bash evidencemodeler.sh -i genome.fa -l -t "
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :a:s:h opt
do
    case $opt in
        a) annotation=$OPTARG;;
        s) speciesName=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${annotation} ]]
then
    usage
fi
if [[ -z ${speciesName} ]]
then
    speciesName=Species_name
fi

firstSpeciesLetter=${speciesName:0:1}
speciesSuffix=${speciesName##*_}
speciesShortName="${firstSpeciesLetter}${speciesSuffix:0:3}"
fileSuffix=${annotation##*annotSplit}
fileNum=${fileSuffix%%.gtf}
fileNumIntZeros=$(printf "%06d" $fileNum)

firstLine=`head -2 ${annotation} | tail -1`
echo "$firstLine"
#firstLine=`head -n 1 ${annotation}`
echo "$firstLine"
prefix=${firstLine##*ID \"}
suffix=${prefix%%\"*}

geneID="FBG${speciesShortName}${fileNumIntZeros}"
entapAnnot=`grep -P "${suffix}\t" simplifiedEntap.tsv`
echo "${entapAnnot}"
EggNOGSeedOrtholog=`echo "$entapAnnot" | awk -F '\t' '{ print $2 }'`
EggNOGPredictedGene=`echo "$entapAnnot" | awk -F '\t' '{ print $3 }'`
EggNOGTaxScope=`echo "$entapAnnot" | awk -F '\t' '{ print $4 }'`
EggNOGDescription=`echo "$entapAnnot" | awk -F '\t' '{ print $5 }'`
if [[ "$EggNOGPredictedGene" == "" ]]; then
    EggNOGPredictedGeneLine2=`echo "$entapAnnot" | awk -F '\t' '{ print $6 }'`
    if [[ "$EggNOGPredictedGeneLine2" == "" ]]; then
        EggNOGPredictedGene="$geneID"
    else
        EggNOGPredictedGeneSuffix=${EggNOGPredictedGeneLine2##*|}
        EggNOGPredictedGene=${EggNOGPredictedGeneSuffix%%_*}
    fi
fi

annotType=`echo "$firstLine" | awk -F '\t' '{ print $2 }'`
firstLineType=`echo "$firstLine" | awk -F '\t' '{ print $3 }'`
OLDIFS=$IFS; IFS=$'\n';
exonIncrement="None"
echo "first type: ${firstLineType}"
if [[ "$annotType" != "tRNAscan-SE" ]]; then
    #it is a transcript or mRNA
    for line in $(<$annotation); do
        chr=`echo "$line" | awk -F '\t' '{ print $1 }'`
        lineType=`echo "$line" | awk -F '\t' '{ print $3 }'`
        middle=`echo "$line" | awk -F '\t' '{ print $4"\t"$5"\t"$6"\t"$7"\t"$8 }'`
        if [[ "$lineType" == "mRNA" ]] && [[ "$firstLineType" == "mRNA" ]]; then
            echo -e "${chr}\tFormBio\tgene\t${middle}\tgene_id \"${geneID}\"; gene_version \"1\"; gene_source \"FormBio\"; gene_biotype \"protein_coding\"; eggnog_ortholog \"${EggNOGSeedOrtholog}\"; gene_name \"${EggNOGPredictedGene}\"; eggnog_taxscope \"${EggNOGTaxScope}\"; gene_desription \"${EggNOGDescription}\";" >> updatedAnnot_${fileNum}.gtf
        elif [[ "$lineType" == "gene" ]]; then
            echo -e "${chr}\tFormBio\t${lineType}\t${middle}\tgene_id \"${geneID}\"; gene_version \"1\"; gene_source \"FormBio\"; gene_biotype \"protein_coding\"; eggnog_ortholog \"${EggNOGSeedOrtholog}\"; gene_name \"${EggNOGPredictedGene}\"; eggnog_taxscope \"${EggNOGTaxScope}\"; gene_desription \"${EggNOGDescription}\";" >> updatedAnnot_${fileNum}.gtf
        elif [[ "$lineType" == "transcript" ]]; then
            echo -e "${chr}\tFormBio\t${lineType}\t${middle}\tgene_id \"${geneID}\"; gene_version \"1\"; transcript_id \"${geneID}.t1\"; transcript_version \"1\"; gene_source \"FormBio\"; gene_biotype \"protein_coding\"; eggnog_ortholog \"${EggNOGSeedOrtholog}\"; gene_name \"${EggNOGPredictedGene}\"; eggnog_taxscope \"${EggNOGTaxScope}\"; gene_desription \"${EggNOGDescription}\";" >> updatedAnnot_${fileNum}.gtf
            
       elif [[ "$lineType" == "exon" ]]; then
            if [[ "$exonIncrement" == "None" ]]; then
                #this is the first exon so we find if its starting low or high
                if [[ "$annotType" == "Helixer" ]]; then
                    echo "Helixer"
                    linePrefix=${line##*.exon.}
                    lineSuffix=${linePrefix%%\";*}
                    echo "suffix"
                    echo "$linePrefix"
                    echo "$lineSuffix"
                elif [[ "$annotType" == "Liftoff" ]]; then
                    echo "Liftoff"
                    echo "line: $line"
            
                    linePrefix=${line##*exon-}
                    lineIntermediateSuffix=${linePrefix%%\"*}
                    echo "prefix: ${linePrefix}"
                    echo "suffix: ${lineIntermediateSuffix}"
                    lineSuffix=${lineIntermediateSuffix##*-}
                    echo "final: ${lineSuffix}"
                else
                    echo "EVM or Augustus"
                    linePrefix=${line##*.exon}
                    lineSuffix=${linePrefix%%\";*}
                                        echo "suffix"
                    echo "$linePrefix"
                    echo "$lineSuffix"
                fi
        echo "Exon num: ${exonNum}"
                exonNum=$(( $lineSuffix + 0 ))
                if [[ $exonNum -eq 1 ]]; then
                    exonIncrement="up"
                else
                    exonIncrement="down"
                fi
            elif [[ "$exonIncrement" == "up" ]]; then
                exonNum=$(( $exonNum + 1 ))
            else
                exonNum=$(( $exonNum - 1 ))
            fi
            echo -e "${chr}\tFormBio\t${lineType}\t${middle}\tgene_id \"${geneID}\"; gene_version \"1\"; transcript_id \"${geneID}.t1\"; transcript_version \"1\"; exon_number \"${exonNum}\"; gene_source \"FormBio\"; gene_biotype \"protein_coding\"; eggnog_ortholog \"${EggNOGSeedOrtholog}\"; gene_name \"${EggNOGPredictedGene}\"; eggnog_taxscope \"${EggNOGTaxScope}\"; gene_desription \"${EggNOGDescription}\";" >> updatedAnnot_${fileNum}.gtf
            
       elif [[ "$lineType" == "cds" ]] || [[ "$lineType" == "CDS" ]]; then
           echo -e "${chr}\tFormBio\t${lineType}\t${middle}\tgene_id \"${geneID}\"; gene_version \"1\"; transcript_id \"${geneID}.t1\"; transcript_version \"1\"; gene_source \"FormBio\"; gene_biotype \"protein_coding\"; eggnog_ortholog \"${EggNOGSeedOrtholog}\"; gene_name \"${EggNOGPredictedGene}\"; eggnog_taxscope \"${EggNOGTaxScope}\"; gene_desription \"${EggNOGDescription}\";" >> updatedAnnot_${fileNum}.gtf
           
       elif [[ "$lineType" == "three_prime_utr" ]]; then
           echo -e "${chr}\tFormBio\t${lineType}\t${middle}\tgene_id \"${geneID}\"; gene_version \"1\"; transcript_id \"${geneID}.t1\"; transcript_version \"1\"; gene_source \"FormBio\"; gene_biotype \"protein_coding\"; eggnog_ortholog \"${EggNOGSeedOrtholog}\"; gene_name \"${EggNOGPredictedGene}\"; eggnog_taxscope \"${EggNOGTaxScope}\"; gene_desription \"${EggNOGDescription}\";" >> updatedAnnot_${fileNum}.gtf
           
       elif [[ "$lineType" == "five_prime_utr" ]]; then
           echo -e "${chr}\tFormBio\t${lineType}\t${middle}\tgene_id \"${geneID}\"; gene_version \"1\"; transcript_id \"${geneID}.t1\"; transcript_version \"1\"; gene_source \"FormBio\"; gene_biotype \"protein_coding\"; eggnog_ortholog \"${EggNOGSeedOrtholog}\"; gene_name \"${EggNOGPredictedGene}\"; eggnog_taxscope \"${EggNOGTaxScope}\"; gene_desription \"${EggNOGDescription}\";" >> updatedAnnot_${fileNum}.gtf
        fi
    done
    
else
    #it is tRNA
    for line in $(<$annotation); do
        chr=`echo "$line" | awk -F '\t' '{ print $1 }'`
        lineType=`echo "$line" | awk -F '\t' '{ print $3 }'`
        middle=`echo "$line" | awk -F '\t' '{ print $4"\t"$5"\t"$6"\t"$7"\t"$8 }'`
        if [[ "$lineType" == "gene" ]]; then
            anticodonPrefix=${line##*anticodon \"}
            anticodonSuffix=${anticodonPrefix%%\"*}
            isotypePrefix=${line##*isotype \"}
            isotypeSuffix=${isotypePrefix%%\"*}
            echo -e "${chr}\tFormBio\t${lineType}\t${middle}\tgene_id \"${geneID}\"; gene_version \"1\"; gene_source \"FormBio\"; anticodon \"anticodonSuffix\"; gene_biotype \"tRNA\"; isotype \"isotypeSuffix\";" >> updatedAnnot_${fileNum}.gtf
        elif [[ "$lineType" == "transcript" ]]; then
            echo -e "${chr}\tFormBio\t${lineType}\t${middle}\tgene_id \"${geneID}\"; gene_version \"1\"; transcript_id \"${geneID}.t1\"; transcript_version \"1\"; gene_source \"FormBio\"; anticodon \"anticodonSuffix\"; gene_biotype \"tRNA\"; isotype \"isotypeSuffix\";" >> updatedAnnot_${fileNum}.gtf
        elif [[ "$lineType" == "exon" ]]; then
            echo -e "${chr}\tFormBio\t${lineType}\t${middle}\tgene_id \"${geneID}\"; gene_version \"1\"; transcript_id \"${geneID}.t1\"; transcript_version \"1\"; exon_number \"${exonNum}\"; gene_source \"FormBio\"; anticodon \"anticodonSuffix\"; gene_biotype \"tRNA\"; isotype \"isotypeSuffix\";" >> updatedAnnot_${fileNum}.gtf
        fi
    done
fi

rm $annotation
