---
layout: post
title: Docker on your local machine
---

Docker is a platform for running applications, it enables us to ensure that we are all running on exactly the version of python, ipython and associated libraries. It can be thought of as a virtual machine, although with many handy additions. 

Firstly, you will need to install docker, the instructions can be found here:
[http://docs.docker.com/mac/started/]
(You’ll need to be root)

The next thing you will need to do is to download the docker container we use for ag1k analysis, we have named this *biipy* (bioinformatics ipython). This can be downloaded by:

    docker pull cggh/biipy:v0.8

v0.8 is the most recent version. Although check [here](https://github.com/cggh/biipy/releases)

Next we should think about what data and code the container has access to when it runs. When we run the container we can map a directory on the host machine to the docker container. I typically map 2 directories (aka volumes).

The first is `/data`, which contains all the ag1k data, and is also where I will write any outputs from my analyses. The directory structure on my host machine mimics the data structure on the cluster. However, it’s not a complete copy as there is obviously insufficient disk on my machine to mirror the cluster filesystem entirely. Instead I just have a subset with some careful use of the very powerful rsync tool. 

Firstly create a “sync” directory somewhere on your system. Home directory is fine.

    mkdir -p ~/sync/data
Then symlink the root /data to your new directory. The reason I do it this way is we run into less issues with permissions, and being in a working directory directly below data allows us to use the `-R` function of `Rsync` neatly.

    ln -s /data ~/sync/data

Now cd to your ~/sync directory and run an rsync command to the part of the filesystem you want to synchronise (the . at the end (meaning current directory) is important!):

    rsync --progress -vrRun \
      clusterpath:/data/anopheles/example/example/example \
      .

The flags are important here:

- `-v` verbose
- `-r` recursive
- `-R` use relative filenames, will ensure the directory structures match
- `-u` only copy newer files, this saves time by not recopying unnecessarily
- `-n` dry run! Remove this flag when you are satisfied with the file list.

The second thing you probably want to map to the container is your code repository. I have all my analysis code (mainly ipython notebooks) checked out in `~/git`. I manage my git pull/push outside of the docker container as normal.

so when it comes to running the container I run something like:

    biipy_version=v0.8
    docker run -d \
      -v ${HOME}/git:/home \
      -v /data:/data \
      -p 31778:8888 \
      --name biipy_$biipy_version \
      -e "docker_image=$biipy_version" \
      cggh/biipy:${biipy_version}

Now fire up a browser and navigate to localhost:31778 (or whatever IP address you chose to set on your machine)

### Caution
Only the mapped volumes will have data persist after you havve finished working, *anything* you write to *anywhere* else on disk will not be kept when the container is exited. 
