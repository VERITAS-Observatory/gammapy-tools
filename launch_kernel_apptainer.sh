#!/bin/bash

# Get the filename of the connection file
connection_file=$(basename $1)
container_name="gammapy_server"
image_name=$GAMMAPY_KERNEL_IMAGE
# Check if the server instance is currently running
if [ $(apptainer instance list $container_name | wc | awk '{print $1}') -ge 2 ]; then
    echo "Server is running"
else
    # if it isn't start an instance
    echo "Starting the server"
    # Bind the runtime-dir (where the connection files are) to /connections within the container
    apptainer instance start --bind $GAMMAPY_WORK_DIR:$GAMMAPY_WORK_DIR --bind `jupyter --runtime-dir`:/connections $image_name $container_name
fi

# Attach and run ipykernel_laucher
apptainer exec instance://$container_name python -m ipykernel_launcher -f /connections/$connection_file
