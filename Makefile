############################################################
#                        CLEAN MAKEFILE                    #
############################################################

PETSC_CFLAGS := $(shell pkg-config --cflags PETSc)
PETSC_LIBS   := $(shell pkg-config --libs PETSc)

FC      = mpifort
SRC_DIR = src
OBJ_DIR = src
TARGET  = pantera.exe

OBJS = src/mpi_common.o src/particle.o src/global.o src/boundary_control_core.o src/screen.o \
       src/mt19937.o src/tools.o src/secondary_electron_emission.o \
       src/grid_and_partition.o src/relativistic_collision_core.o src/collisions.o src/fields.o \
       src/washboard.o src/postprocess.o src/initialization.o \
       src/timecycle.o src/pantera.o

FCFLAGS = -cpp -fimplicit-none -O3 $(PETSC_CFLAGS)
LDFLAGS = $(PETSC_LIBS)

.PHONY: all clean test test-boundary-core test-relativistic-collision-core test-static-p0

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) -o $@ $^ $(LDFLAGS) -lm

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FCFLAGS) -c -J$(OBJ_DIR) -o $@ $<

# Module dependencies
src/global.o: src/global.f90 src/particle.o src/mpi_common.o
src/boundary_control_core.o: src/boundary_control_core.f90
src/relativistic_collision_core.o: src/relativistic_collision_core.f90
src/collisions.o: src/collisions.f90 src/global.o src/mpi_common.o src/screen.o src/tools.o src/relativistic_collision_core.o
src/fields.o: src/fields.f90 src/global.o src/mpi_common.o src/screen.o src/tools.o src/grid_and_partition.o src/secondary_electron_emission.o src/boundary_control_core.o
src/postprocess.o: src/postprocess.f90 src/global.o src/mpi_common.o src/screen.o src/fields.o
src/initialization.o: src/initialization.f90 src/global.o src/mpi_common.o src/screen.o src/tools.o src/secondary_electron_emission.o src/grid_and_partition.o src/collisions.o src/fields.o src/washboard.o
src/timecycle.o: src/timecycle.f90 src/global.o src/mpi_common.o src/screen.o src/tools.o src/secondary_electron_emission.o src/grid_and_partition.o src/collisions.o src/fields.o src/washboard.o src/postprocess.o
src/pantera.o: src/pantera.f90 src/global.o src/mpi_common.o src/screen.o src/initialization.o src/timecycle.o
src/tools.o: src/tools.f90 src/global.o src/mpi_common.o src/mt19937.o
src/secondary_electron_emission.o: src/secondary_electron_emission.f90 src/global.o src/mpi_common.o src/tools.o
src/grid_and_partition.o: src/grid_and_partition.f90 src/global.o src/mpi_common.o src/screen.o src/tools.o
src/washboard.o: src/washboard.f90 src/global.o src/mpi_common.o src/tools.o
src/screen.o: src/screen.f90 src/mpi_common.o

test: test-boundary-core test-relativistic-collision-core test-static-p0

test-boundary-core:
	mkdir -p /tmp/pantera-unit
	gfortran -std=f2008 -O0 -g -Wall -Wextra -Werror -fcheck=all \
	  -ffpe-trap=invalid,zero,overflow -J/tmp/pantera-unit \
	  src/boundary_control_core.f90 tests/test_boundary_control_core.f90 \
	  -o /tmp/pantera-unit/test_boundary_control_core
	/tmp/pantera-unit/test_boundary_control_core

test-relativistic-collision-core:
	mkdir -p /tmp/pantera-unit
	gfortran -std=f2008 -O0 -g -Wall -Wextra -Werror -fcheck=all \
	  -ffpe-trap=invalid,zero,overflow -J/tmp/pantera-unit \
	  src/relativistic_collision_core.f90 tests/test_relativistic_collision_core.f90 \
	  -o /tmp/pantera-unit/test_relativistic_collision_core
	/tmp/pantera-unit/test_relativistic_collision_core

test-static-p0:
	python3 -m unittest tests/test_static_p0.py

clean:
	rm -rf $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(TARGET)
