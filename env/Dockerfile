FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gfortran \
    make \
        cmake \
    pkg-config \
    openmpi-bin \
    libopenmpi-dev \
    libpetsc-real3.19-dev \
    python3 \
    python3-tk \
    python3-numpy \
    python3-matplotlib \
    python3-pandas \
    python3-scipy \
    ca-certificates \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /work
CMD ["/bin/bash"]
