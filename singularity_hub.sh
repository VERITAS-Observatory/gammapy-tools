#!/bin/bash

[[ "$1" ]] && DIR=$1  || DIR=$(pwd)
apptainer run -B $DIR:/local_data gamma-tools.sif
#bash $(pwd)/launch_jupyter.sh
