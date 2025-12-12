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

clean:
	rm -rf $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(TARGET)
