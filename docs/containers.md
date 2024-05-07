

# Working with containers

At some stage you might reach one of these situations:

1. You have a new/old machine and installing all the various dependencies will be a pain.
2. You have exhausted the power of your personal machine, and you'd like to move to cluster. Installing all the various dependencies will be a pain.

The solution to both of these problems are "[Containers](https://www.techtarget.com/searchitoperations/definition/container-containerization-or-container-based-virtualization)". If you are working on your own personal machine you might decide to use "[Docker](https://docs.docker.com/get-started/overview/)". If you are working on a cluster you might be surprised to find out that HPC centres/clusters don't really like people using docker, instead they prefer people use something like [Singularity/Apptainer](https://hsf-training.github.io/hsf-training-singularity-webpage/). Why?

1. HPC clusters use queue systems (e.g. slurm) to manage resources and schedule jobs. A docker command is a call to the docker-daemon which manages resources, this adds an extra layer of complexity for scheduling and resource management. Singularity/Apptainer doesn't require a daemon, instead it spawns a new process. This allows for easy management of the resources available to that process.

2. Permissions. By default, Docker will give the user superuser privileges. This is because Docker is designed to run services on machines that aren't "multi-user". In a multi-user environment this could have very destructive. Singularity/Apptainer is a process launched by a user which maintains the permissions of that user.


## Should you use Docker or Singularity/Apptainer?

Docker is widely used in industry. Mastering Docker will go a long way outside an academic setting. Apptainer is gaining wider acceptance in the HPC community. However, they know that most people will still work with Docker images. One can convert from a Docker image to a singularity image. An example of how to do so can be found in the `create_singularity.sh` image.

Mastering Apptainer will help you create images with smaller footprints. At the moment it will mainly benefit scientific computing rather than in-industry computing.

## Building and running the Docker image

To create the docker image run:
```
docker build -t local/gammapy-tools:latest .
```

Alternatively if you want to start a jupyter lab instance:
```
docker compose up --detach
```
Will automatically build the image, launch it and start a jupyter lab on port 8888 with the password `letmein`. This is all handled by the `docker-compose.yaml` file and `Dockerfile`.

To run an interactive shell:
```
docker run --rm -it local/gammapy-tools:latest bash
```


## Building and running the Singularity image

To build the image using singularity/apptainer run:
```
singularity/apptainer build gammapy-tools.sif singularity.def
```

Which creates the `gammapy-tools.sif` image file, at the time of writing the file is ~700MB.


To launch a jupyter lab instance with this image run:
```
singularity/apptainer run -B $(pwd):/local_data  gammapy-tools.sif
```

We can now open a jupyter lab instance in our browser at port 1234.

You'll notice that you cannot access /local_data from the jupyter notebook. The solution is to `cd /local_data` before launching the jupyter lab:
```
cd /local_data
jupyter-notebook
```

To simplify all of this, there are two scripts `singularity_hub.sh`:
```
#!/bin/bash

[[ "$1" ]] && DIR=$1  || DIR=$(pwd)
singularity exec -B $DIR:/local_data gammapy-tools.sif bash $(pwd)/launch_jupyter.sh
```
and `launch_jupyer.sh`:
```
cd /local_data
jupyter-notebook
```

Using these two scripts we can launch a jupyter lab instance with a directory of our choice:
```
singularity_hub.sh /path/to/my/data
```
Will start the singularity image, mounting `/path/to/my/data` to `/local_data`, `cd` into `/local_data` and setup up a jupyter lab.

We can also run standalone scripts using the singularity image:
```
singularity/apptainer run -B /path/to/my/data:/local_data  gammapy-tools.sif python my_super_awesome_script.py
```
This will have access to data mounted in `/local_data` and run the python version within the container.

Finally, you can start an interactive shell using:
```
singularity/apptainer shell -B /path/to/my/data:/local_data gammapy-tools.sif
```