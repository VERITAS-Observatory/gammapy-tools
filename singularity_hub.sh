#!/bin/bash

[[ "$1" ]] && DIR=$1  || DIR=$(pwd)
apptainer run -B $DIR:/local_data gammapy-tools.sif
#bash $(pwd)/launch_jupyter.sh
