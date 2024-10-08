#!/bin/bash

# Get the filename of the connection file
connection_file=$(basename $1)
container_name="gammapy_server"
image_name=$GAMMAPY_KERNEL_IMAGE

# Check if the server is currently runnong
if [ "$( docker container inspect -f '{{.State.Running}}' $container_name )" = "true" ]; then
    echo "Server is running"
else
    # if not then start the server
    echo "Starting the server"
    # Get the path of where the connection files will be stored
    connection_path="$(jupyter --runtime-dir)"
    # Stop and remove the container if it already exists
    docker stop $container_name
    docker rm $container_name
    # Create a new container
    # Add it to the host network and mounting the connection_path
    docker create --name=$container_name -v  $connection_path:/connections --network=host  $image_name
    # Start this server instance
    docker start  $container_name
fi

# Launch a ipykernel using the connection file that was passed as arg 1
docker exec  $container_name python -m ipykernel_launcher -f /connections/$connection_file
