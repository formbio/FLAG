#!/bin/bash

# cd to the directory
cd containers/singularity_flag/

#build the singularity image
singularity build singularity_flag.image singularity_flag.def

# cd to examples
cd ../../examples

# move the image
mv ../containers/singularity_flag/singularity_flag.image .
