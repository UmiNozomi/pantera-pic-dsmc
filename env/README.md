# Fast Environment Setup

This project builds most reliably on Ubuntu 24.04 / WSL2 with MPI and PETSc
provided by apt packages.

## Option A: WSL / Ubuntu Native

From the project root:

```bash
bash env/setup_ubuntu24.sh
```

This installs the build dependencies, checks `mpifort`, `pkg-config PETSc`, and
Python plotting packages, then builds `pantera.exe`.

Run a case:

```bash
cp pantera.exe 2026/H2_low_recycling/
cd 2026/H2_low_recycling
mpirun -np 10 ./pantera.exe
```

If you only want to check an existing machine:

```bash
bash env/check_env.sh
make -j$(nproc)
```

## Option B: Docker

Build a reusable image:

```bash
docker build -f env/Dockerfile -t pantera-env .
```

Open a shell inside the image with the project mounted:

```bash
docker run --rm -it -v "$PWD:/work" pantera-env
```

Inside the container:

```bash
make -j$(nproc)
cp pantera.exe 2026/H2_low_recycling/
cd 2026/H2_low_recycling
mpirun -np 4 ./pantera.exe
```

Docker is good for reproducible builds. For long production runs on WSL or HPC,
native Ubuntu/WSL or Apptainer/Singularity is usually more convenient.
