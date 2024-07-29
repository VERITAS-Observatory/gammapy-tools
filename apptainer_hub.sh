#!/bin/bash

[[ "$1" ]] && DIR=$1  || DIR=$(pwd)
if command -v apptainer &> /dev/null
then
    apptainer run -B $DIR:/local_data gammapy-tools.sif
else
    singularity run -B $DIR:/local_data gammapy-tools.sif
fi

