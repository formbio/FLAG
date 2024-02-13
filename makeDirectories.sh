# If you are running FLAG as in the example runs then your output directory will be named as outputdir
# In this directory for some runs fake files are needed to substitute for other inputs that arent provided
# This script simply makes those fake files in the directory defined as outDir
outDir='outputdir'
mkdir "$outDir"
touch "${outDir}/emptyProteinPlaceHolder.txt"
touch "${outDir}/emptyTranscriptPlaceHolder.txt"
touch "${outDir}/emptyPlaceHolder.txt"
