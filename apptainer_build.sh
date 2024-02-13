#!/bin/bash

if command -v apptainer &> /dev/null
then
    apptainer build ./gammapy-tools.sif gammapy-tools.def
else
    singularity build ./gammapy-tools.sif gammapy-tools.def
fi
