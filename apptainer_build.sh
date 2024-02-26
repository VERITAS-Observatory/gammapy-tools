#!/bin/bash

[[ "$1" ]] && IMAGE_NAME=$1  || IMAGE_NAME="gammapy-tools.sif"


if command -v apptainer &> /dev/null
then
    apptainer build ./$IMAGE_NAME gammapy-tools.def
else
    singularity build ./$IMAGE_NAME gammapy-tools.def
fi
