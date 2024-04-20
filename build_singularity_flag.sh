#!/bin/bash

# cd to the directory
echo "Entering singularity_flag container directory"
cd containers/singularity_flag/

#build the singularity image
echo "Building the singularity_flag singularity image"
singularity build singularity_flag.image singularity_flag.def

# cd to examples
echo "Entering the examples directory"
cd ../../examples

# move the image
echo "Moving the singularity_flag singularity image to the examples directory"
mv ../containers/singularity_flag/singularity_flag.image .

# Setup the run directory
echo "Creating initial files/directories needed to run flag from the singularity image"
mkdir nxf_temp
mkdir nxf_home; cd nxf_home; mkdir framework; cd framework; mkdir 23.10.0; cd 23.10.0
wget https://www.nextflow.io/releases/v23.10.0/nextflow-23.10.0-one.jar
cd ../../../
touch pipeline_trace.txt
mkdir .nextflow
mkdir tempdir

echo "Singularity FLAG image built and initial files setup in the examples directory."

