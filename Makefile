############################################################
#                        CLEAN MAKEFILE                    #
############################################################

PETSC_CFLAGS := $(shell pkg-config --cflags PETSc)
PETSC_LIBS   := $(shell pkg-config --libs PETSc)

FC      = mpifort
SRC_DIR = src
OBJ_DIR = src
TARGET  = pantera.exe

OBJS = src/mpi_common.o src/particle.o src/global.o src/screen.o \
       src/mt19937.o src/tools.o src/secondary_electron_emission.o \
       src/grid_and_partition.o src/collisions.o src/fields.o \
       src/washboard.o src/postprocess.o src/initialization.o \
       src/timecycle.o src/pantera.o

FCFLAGS = -cpp -fimplicit-none -O3 $(PETSC_CFLAGS)
LDFLAGS = $(PETSC_LIBS)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) -o $@ $^ $(LDFLAGS) -lm

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FCFLAGS) -c -J$(OBJ_DIR) -o $@ $<

# Module dependencies
src/global.o: src/global.f90 src/particle.o src/mpi_common.o
src/collisions.o: src/collisions.f90 src/global.o src/mpi_common.o src/screen.o src/tools.o
src/fields.o: src/fields.f90 src/global.o src/mpi_common.o src/screen.o
src/postprocess.o: src/postprocess.f90 src/global.o src/mpi_common.o src/screen.o
src/initialization.o: src/initialization.f90 src/global.o src/mpi_common.o src/screen.o src/tools.o src/secondary_electron_emission.o src/grid_and_partition.o src/collisions.o src/fields.o src/washboard.o
src/timecycle.o: src/timecycle.f90 src/global.o src/mpi_common.o src/screen.o src/tools.o src/secondary_electron_emission.o src/grid_and_partition.o src/collisions.o src/fields.o src/washboard.o src/postprocess.o
src/pantera.o: src/pantera.f90 src/global.o src/mpi_common.o src/screen.o src/initialization.o src/timecycle.o
src/tools.o: src/tools.f90 src/global.o src/mpi_common.o src/mt19937.o
src/secondary_electron_emission.o: src/secondary_electron_emission.f90 src/global.o src/mpi_common.o src/tools.o
src/grid_and_partition.o: src/grid_and_partition.f90 src/global.o src/mpi_common.o src/screen.o src/tools.o
src/washboard.o: src/washboard.f90 src/global.o src/mpi_common.o src/tools.o
src/screen.o: src/screen.f90 src/mpi_common.o

clean:
	rm -rf $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(TARGET)
