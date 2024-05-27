#!/bin/bash

[[ "$1" ]] && DOCKER_NAME=$1  || DOCKER_NAME=$(whoami)

docker build . -t $DOCKER_NAME/gammapy_tools:latest
