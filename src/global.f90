! Copyright (C) 2025 von Karman Institute for Fluid Dynamics (VKI)
!
! This file is part of PANTERA PIC-DSMC, a software for the simulation
! of rarefied gases and plasmas using particles.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.PANTERA PIC-DSMC

! This module holds global variables

MODULE global

   USE particle
   USE mpi_common


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Constants !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   REAL(KIND=8) :: PI   = 3.1415926535897932d0           ! https://dlmf.nist.gov/3.12
   REAL(KIND=8) :: EPS0 = 8.8541878128d-12               ! https://physics.nist.gov/cgi-bin/cuu/Value?ep0
   REAL(KIND=8) :: MU0  = 1.25663706212d-6               ! https://physics.nist.gov/cgi-bin/cuu/Value?mu0
   REAL(KIND=8) :: KB   = 1.380649d-23                   ! https://physics.nist.gov/cgi-bin/cuu/Value?k
   REAL(KIND=8) :: QE   = 1.602176634d-19                ! https://physics.nist.gov/cgi-bin/cuu/Value?e
   REAL(KIND=8) :: NA   = 6.02214076e23                  ! https://physics.nist.gov/cgi-bin/cuu/Value?na
   REAL(KIND=8) :: ME   = 9.1093837139d-31               ! https://physics.nist.gov/cgi-bin/cuu/Value?me

   REAL(KIND=8) :: EPS_SCALING = 1.d0
      
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Particle variables and arrays !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   INTEGER :: NP_PROC = 0 ! Number of particles in local MPI process
   INTEGER :: MPI_PARTICLE_DATA_STRUCTURE ! We need this for MPI
   TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: particles, part_dump, part_inject
   INTEGER :: NP_DUMP_PROC = 0
   INTEGER :: NP_INJECT_PROC = 0
   INTEGER :: DUMP_TRAJECTORY_START = -1
   INTEGER :: DUMP_TRAJECTORY_NUMBER = 0
   CHARACTER*256 :: TRAJDUMP_SAVE_PATH
   CHARACTER*256 :: PARTDUMP_SAVE_PATH
   CHARACTER*256 :: CHECKS_SAVE_PATH
   CHARACTER*256 :: RESIDUAL_SAVE_PATH
   CHARACTER*256 :: INJECT_FILENAME
   REAL(KIND=8) :: INJECT_PROBABILITY = 1
   LOGICAL :: BOOL_INJECT_FROM_FILE = .FALSE.

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!! Fluid electrons !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   REAL(KIND=8) :: BOLTZ_N0 = 0., BOLTZ_PHI0 = 0., BOLTZ_TE = 0.
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: BOLTZ_NRHOE
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: BOLTZ_SOLID_NODES

   ! Option for Kappa distribution
   LOGICAL :: BOOL_KAPPA_FLUID    = .FALSE.
   REAL(KIND=8) :: KAPPA_FLUID_C  = 4.d0

   LOGICAL :: BOOL_CONDUCTIVE_BC = .FALSE.
   REAL(KIND=8) :: WALL_METAL_POTENTIAL

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Geometry, domain and grid !!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   LOGICAL      :: BOOL_AXI = .FALSE. ! Assign default value!
   INTEGER      :: NX, NY, NZ
   INTEGER      :: DIMS = 2
   REAL(KIND=8) :: XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
   REAL(KIND=8) :: CELL_VOL
   LOGICAL      :: AXI = .FALSE.
   LOGICAL      :: BOOL_X_PERIODIC = .FALSE. ! Assign default value!
   LOGICAL      :: BOOL_Y_PERIODIC = .FALSE.
   LOGICAL      :: BOOL_Z_PERIODIC = .FALSE.
   LOGICAL, DIMENSION(4) :: BOOL_PERIODIC = .FALSE.
   LOGICAL, DIMENSION(4) :: BOOL_SPECULAR = .FALSE.
   LOGICAL, DIMENSION(4) :: BOOL_REACT    = .FALSE.
   LOGICAL, DIMENSION(4) :: BOOL_DIFFUSE  = .FALSE.
   INTEGER, DIMENSION(:), ALLOCATABLE :: BOUNDARY_COLL_COUNT, WALL_COLL_COUNT, LINE_EMIT_COUNT
   REAL(KIND=8) :: BOUNDTEMP
   CHARACTER(LEN=256) :: GRID_FILENAME
   LOGICAL :: BOOL_BINARY_OUTPUT = .TRUE.
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: XCOORD, YCOORD
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: XSIZE, YSIZE
   INTEGER, DIMENSION(:), ALLOCATABLE :: CELL_PROCS
   INTEGER :: NCELLS, NNODES, NBOUNDCELLS, NBOUNDNODES
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: CELL_VOLUMES


   ENUM, BIND(C)
      ENUMERATOR RECTILINEAR_UNIFORM, RECTILINEAR_NONUNIFORM, QUADTREE, UNSTRUCTURED
   END ENUM
   INTEGER(KIND(RECTILINEAR_UNIFORM)) :: GRID_TYPE = RECTILINEAR_UNIFORM


   ENUM, BIND(C)
      ENUMERATOR STRIPSX, STRIPSY, STRIPSZ, SLICESX, SLICESZ
   END ENUM
   INTEGER(KIND(STRIPSX)) :: PARTITION_STYLE = STRIPSX

   LOGICAL :: LOAD_BALANCE = .FALSE.
   INTEGER :: LOAD_BALANCE_EVERY = 0

   TYPE UNSTRUCTURED_0D_GRID_DATA_STRUCTURE
      INTEGER :: NUM_NODES, NUM_POINTS
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: NODE_COORDS
      INTEGER, DIMENSION(:), ALLOCATABLE        :: POINT_PG
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: PG_NODES
      INTEGER, DIMENSION(:), ALLOCATABLE      :: POINT_NODES
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: VERTEX_AREAS
   END TYPE UNSTRUCTURED_0D_GRID_DATA_STRUCTURE

   TYPE(UNSTRUCTURED_0D_GRID_DATA_STRUCTURE) :: U0D_GRID


   TYPE UNSTRUCTURED_1D_GRID_DATA_STRUCTURE
      INTEGER :: NUM_NODES, NUM_CELLS
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: NODE_COORDS
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: CELL_NODES
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: CELL_NEIGHBORS
      REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: EDGE_NORMAL
      INTEGER, DIMENSION(:), ALLOCATABLE        :: CELL_PG
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: CELL_EDGES_PG
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: PG_NODES
      REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: BASIS_COEFFS
      INTEGER, DIMENSION(:), ALLOCATABLE        :: PERIODIC_RELATED_NODE
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: SEGMENT_NODES_BOUNDARY_INDEX
      INTEGER, DIMENSION(:), ALLOCATABLE        :: NODES_BOUNDARY_INDEX
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: CELL_VOLUMES
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: SEGMENT_AREAS
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: SEGMENT_LENGTHS
   END TYPE UNSTRUCTURED_1D_GRID_DATA_STRUCTURE

   TYPE(UNSTRUCTURED_1D_GRID_DATA_STRUCTURE) :: U1D_GRID


   TYPE UNSTRUCTURED_2D_GRID_DATA_STRUCTURE
      INTEGER :: NUM_NODES, NUM_CELLS, NUM_LINES
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: NODE_COORDS
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: CELL_NODES
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: LINE_NODES
      INTEGER, DIMENSION(:), ALLOCATABLE        :: POINT_NODES
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: CELL_NEIGHBORS
      REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: EDGE_NORMAL
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: CELL_EDGES_LEN
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: CELL_EDGES_PG
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: CELL_EDGES_BOUNDARY_INDEX
      INTEGER, DIMENSION(:), ALLOCATABLE        :: NODES_BOUNDARY_INDEX
      INTEGER, DIMENSION(:), ALLOCATABLE        :: CELL_PG
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: PG_NODES
      REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: BASIS_COEFFS
      INTEGER, DIMENSION(:), ALLOCATABLE        :: PERIODIC_RELATED_NODE
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: CELL_VOLUMES
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: CELL_AREAS
   END TYPE UNSTRUCTURED_2D_GRID_DATA_STRUCTURE

   TYPE(UNSTRUCTURED_2D_GRID_DATA_STRUCTURE) :: U2D_GRID


   TYPE UNSTRUCTURED_3D_GRID_DATA_STRUCTURE
      INTEGER :: NUM_NODES, NUM_CELLS, NUM_FACES
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: NODE_COORDS
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: CELL_NODES
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: CELL_FACES
      REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: CELL_FACES_COEFFS
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE    :: FACE_NODES
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: CELL_NEIGHBORS
      REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: FACE_NORMAL
      REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: FACE_TANG1
      REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: FACE_TANG2
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: FACE_AREA
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: CELL_FACES_PG
      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: CELL_FACES_BOUNDARY_INDEX
      INTEGER, DIMENSION(:), ALLOCATABLE        :: NODES_BOUNDARY_INDEX
      INTEGER, DIMENSION(:), ALLOCATABLE        :: CELL_PG
      REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: BASIS_COEFFS
      INTEGER, DIMENSION(:), ALLOCATABLE        :: PERIODIC_RELATED_NODE
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: CELL_VOLUMES
   END TYPE UNSTRUCTURED_3D_GRID_DATA_STRUCTURE

   TYPE(UNSTRUCTURED_3D_GRID_DATA_STRUCTURE) :: U3D_GRID

   !CHARACTER*256 :: MESH_FILENAME = 'meshlofthousecoarse.su2'
   !CHARACTER*256 :: MESH_FILENAME = 'meshlofthousenew.su2'
   !CHARACTER*256 :: MESH_FILENAME = 'meshtestfine.su2'
   !CHARACTER*256 :: MESH_FILENAME = 'meshrectfine.su2'
   !CHARACTER*256 :: MESH_FILENAME = 'meshtestbound.su2'
   !CHARACTER*256 :: MESH_FILENAME = 'meshlh.su2'
   CHARACTER*256 :: MESH_FILENAME

   INTEGER         :: N_GRID_BC = 0

   ENUM, BIND(C)
      ENUMERATOR VACUUM, SPECULAR, DIFFUSE, CLL, REACT, AXIS, PERIODIC_MASTER, PERIODIC_SLAVE, EMIT, WB_BC
   END ENUM

   ENUM, BIND(C)
      ENUMERATOR DIRICHLET_BC, NEUMANN_BC, DIELECTRIC_BC, CONDUCTIVE_BC, ROBIN_BC, PERIODIC_MASTER_BC, PERIODIC_SLAVE_BC, &
                 RF_VOLTAGE_BC, DECOUPLED_RF_VOLTAGE_BC, SPICE_NODE_BC, NO_BC
   END ENUM

   ENUM, BIND(C)
      ENUMERATOR FLUID, SOLID
   END ENUM

   TYPE BOUNDARY_CONDITION_DATA_STRUCTURE
      CHARACTER(LEN=256)       :: PHYSICAL_GROUP_NAME
      INTEGER(KIND(VACUUM))    :: PARTICLE_BC = VACUUM
      INTEGER(KIND(DIRICHLET_BC)) :: FIELD_BC = NO_BC
      INTEGER(KIND(FLUID)) :: VOLUME_BC = FLUID

      REAL(KIND=8) :: WALL_TEMP
      REAL(KIND=8) :: WALL_POTENTIAL
      REAL(KIND=8) :: WALL_EFIELD
      REAL(KIND=8) :: ACC_N
      REAL(KIND=8) :: ACC_T

      REAL(KIND=8) :: EPS_REL

      REAL(KIND=8) :: WALL_RF_POTENTIAL
      REAL(KIND=8) :: RF_FREQUENCY
      REAL(KIND=8) :: CAPACITANCE

      REAL(KIND=8) :: SPICE_NODE_POTENTIAL = 0.d0
      REAL(KIND=8) :: SPICE_NODE_CURRENT

      REAL(KIND=8), DIMENSION(2) :: TRANSLATEVEC

      LOGICAL :: REACT = .FALSE.
      LOGICAL :: DUMP_FLUXES = .FALSE.

      ! Washboard model
      REAL(KIND=8) :: A
      REAL(KIND=8) :: B
      REAL(KIND=8) :: W
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: MAX_P_DN
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: MAX_P_UP
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: P_COLL_UP

      ! Constant current control parameters
      LOGICAL :: IS_CONSTANT_CURRENT = .FALSE.
      REAL(KIND=8) :: TARGET_CURRENT            ! Target current [A]
      REAL(KIND=8) :: INITIAL_VOLTAGE           ! Initial voltage for CV mode [V]
      REAL(KIND=8) :: PID_KP, PID_KI, PID_KD    ! PID coefficients
      INTEGER :: SLIDING_WINDOW_SIZE = 10        ! Number of timesteps for averaging
      
      ! State machine
      LOGICAL :: CC_MODE_ACTIVE = .FALSE.        ! FALSE = CV mode, TRUE = CC mode
      
      ! Sliding window for current filtering
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: CURRENT_WINDOW_ION
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: CURRENT_WINDOW_ELEC
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: CURRENT_WINDOW_SEE
      INTEGER :: WINDOW_INDEX = 0
      
      ! Current statistics (per timestep, per boundary)
      REAL(KIND=8) :: TIMESTEP_CHARGE_ION = 0.d0   ! Ion charge accumulated this timestep [C]
      REAL(KIND=8) :: TIMESTEP_CHARGE_ELEC = 0.d0  ! Electron charge accumulated this timestep [C]
      REAL(KIND=8) :: TIMESTEP_CHARGE_SEE = 0.d0   ! SEE charge accumulated this timestep [C]
      
      ! PID controller state
      REAL(KIND=8) :: ERROR_INTEGRAL = 0.d0
      REAL(KIND=8) :: ERROR_PREV = 0.d0
      
      ! Voltage safety limits
      REAL(KIND=8) :: VOLTAGE_MIN = -1.0d10     ! Minimum allowed voltage [V]
      REAL(KIND=8) :: VOLTAGE_MAX = +1.0d10     ! Maximum allowed voltage [V]
      LOGICAL :: APPLY_VOLTAGE_LIMITS = .FALSE. ! Enable voltage clamping

      ! Boundary-specific wall reactions
      INTEGER :: N_WALL_REACTIONS = 0
      TYPE(WALL_REACTIONS_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: WALL_REACTIONS

   END TYPE BOUNDARY_CONDITION_DATA_STRUCTURE

   TYPE(BOUNDARY_CONDITION_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: GRID_BC



   TYPE EMIT_TASK_DATA_STRUCTURE
      REAL(KIND=8) :: NRHO
      REAL(KIND=8) :: UX, UY, UZ
      REAL(KIND=8) :: TTRA, TROT, TVIB
      REAL(KIND=8) :: U_NORM
      INTEGER      :: MIX_ID
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: NFS
      INTEGER      :: IC
      INTEGER      :: IV1, IV2
      INTEGER      :: IFACE
   END TYPE EMIT_TASK_DATA_STRUCTURE

   TYPE(EMIT_TASK_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: EMIT_TASKS
   INTEGER :: N_EMIT_TASKS
   



   TYPE INITIAL_PARTICLES_DATA_STRUCTURE
      REAL(KIND=8) :: NRHO
      REAL(KIND=8) :: UX, UY, UZ
      REAL(KIND=8) :: TTRAX, TTRAY, TTRAZ, TROT, TVIB
      INTEGER      :: MIX_ID
   END TYPE INITIAL_PARTICLES_DATA_STRUCTURE

   TYPE(INITIAL_PARTICLES_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: INITIAL_PARTICLES_TASKS
   INTEGER :: N_INITIAL_PARTICLES_TASKS
   

   TYPE VOLUME_INJECT_DATA_STRUCTURE
      REAL(KIND=8) :: NRHODOT
      REAL(KIND=8) :: UX, UY, UZ
      REAL(KIND=8) :: TTRAX, TTRAY, TTRAZ, TROT, TVIB
      INTEGER      :: MIX_ID
   END TYPE VOLUME_INJECT_DATA_STRUCTURE

   TYPE(VOLUME_INJECT_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: VOLUME_INJECT_TASKS
   INTEGER :: N_VOLUME_INJECT_TASKS
   

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Electromagnetic fields !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   ENUM, BIND(C)
      ENUMERATOR NONE, EXPLICIT, SEMIIMPLICIT, FULLYIMPLICIT, EXPLICITLIMITED, HYBRID
   END ENUM
   INTEGER(KIND(EXPLICIT)) :: PIC_TYPE = NONE

   INTEGER :: NPX, NPY
   REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: E_FIELD, B_FIELD
   REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: EBAR_FIELD
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: PHI_FIELD
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: PHIBAR_FIELD
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: PHI_FIELD_OLD
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: PHI_FIELD_NEW
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: RHS
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: DXLDRATIO
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: RHS_NEW
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: J_FIELD
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: MASS_MATRIX
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: DIRICHLET
   LOGICAL, DIMENSION(:), ALLOCATABLE :: IS_DIRICHLET
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: NEUMANN
   LOGICAL, DIMENSION(:), ALLOCATABLE :: IS_NEUMANN
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: SURFACE_CHARGE 


   REAL(KIND=8), DIMENSION(3) :: EXTERNAL_B_FIELD = 0
   ! This is used for defining static magnetic fields from solenoids.
   INTEGER         :: N_SOLENOIDS = 0

   TYPE SOLENOID
      REAL(KIND=8) :: X1, Y1, X2, Y2
      REAL(KIND=8) :: WIRE_CURRENT
      INTEGER :: N_WIRES_X, N_WIRES_Y
   END TYPE SOLENOID

   TYPE(SOLENOID), DIMENSION(:), ALLOCATABLE :: SOLENOIDS

   INTEGER         :: N_MAGNETS = 0

   TYPE MAGNET
      REAL(KIND=8) :: X1, Y1, X2, Y2
      REAL(KIND=8) :: STRENGTH
   END TYPE MAGNET

   TYPE(MAGNET), DIMENSION(:), ALLOCATABLE :: MAGNETS

   ! Magnetic field from magnetic dipole
   LOGICAL :: BOOL_MAGNETIC_DIPOLE = .FALSE.
   REAL(KIND=8), DIMENSION(3) :: DIPOLE_POSITION, DIPOLE_ORIENTATION
   REAL(KIND=8) :: MAGNETIC_MOMENT

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Numerical settings !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   INTEGER      :: NT, tID
   INTEGER      :: RESTART_TIMESTEP = 0
   LOGICAL      :: SAVE_INITIAL_TIMESTEP = .FALSE.
   REAL(KIND=8) :: FNUM, DT, START_CPU_TIME
   INTEGER(KIND=8) :: RNG_SEED_GLOBAL, RNG_SEED_LOCAL
   INTEGER      :: DUMP_PART_EVERY = 1
   INTEGER      :: DUMP_PART_START = 0
   INTEGER      :: DUMP_PART_BOUND_EVERY = 1
   INTEGER      :: DUMP_PART_BOUND_START = 0
   INTEGER      :: DUMP_GRID_AVG_EVERY = 1
   INTEGER      :: DUMP_GRID_START = 0
   INTEGER      :: DUMP_GRID_N_AVG = 1
   INTEGER      :: DUMP_BOUND_AVG_EVERY = 1
   INTEGER      :: DUMP_BOUND_START = 0
   INTEGER      :: DUMP_BOUND_N_AVG = 1
   REAL(KIND=8) :: PARTDUMP_FRACSAMPLE = 1
   LOGICAL      :: PERFORM_CHECKS = .FALSE.
   INTEGER      :: CHECKS_EVERY = 1
   INTEGER      :: STATS_EVERY = 1
   INTEGER      :: TIMING_STATS_EVERY = 100
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: CELL_FNUM
   LOGICAL      :: BOOL_RADIAL_WEIGHTING = .FALSE.
   INTEGER      :: JACOBIAN_TYPE = 1
   LOGICAL      :: COLOCATED_ELECTRONS = .FALSE.
   REAL(KIND=8) :: COLOCATED_ELECTRONS_TTRA = 11600.d0
   LOGICAL      :: RESIDUAL_AND_JACOBIAN_COMBINED = .FALSE.
   REAL(KIND=8) :: SNES_RTOL = 1.d-9

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Collisions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !LOGICAL           :: BOOL_MCC = .FALSE., BOOL_DSMC = .FALSE., BOOL_BGK = .FALSE.
   ENUM, BIND(C)
   ENUMERATOR NO_COLL, DSMC, BGK, MCC, MCC_VAHEDI, DSMC_VAHEDI
   END ENUM
   INTEGER(KIND(DSMC)) :: COLLISION_TYPE = NONE

   REAL(KIND=8)      :: MCC_BG_DENS, MCC_BG_TTRA
   INTEGER           :: MCC_BG_MIX
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: MCC_BG_CELL_NRHO
   LOGICAL           :: BOOL_BG_DENSITY_FILE = .FALSE.
   INTEGER           :: DSMC_COLL_MIX
   INTEGER           :: TIMESTEP_COLL
   INTEGER           :: TIMESTEP_REAC
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: VSS_GREFS ! Matrix of reference relative velocities for VSS
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: VSS_SIGMAS ! Matrix of reference cross sections for VSS
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: VSS_ALPHAS ! Matrix of reference scattering coeff. for VSS
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: VSS_OMEGAS ! Matrix of reference temperature exponent for VSS
   REAL(KIND=8) :: SIGMAMAX = 0

   LOGICAL           :: BOOL_THERMAL_BATH = .FALSE.
   REAL(KIND=8)      :: TBATH

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Spectral Diagnostics !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   LOGICAL :: BOOL_SPECTRAL_DIAGNOSTICS = .FALSE.
   REAL(KIND=8) :: SPECTRAL_SAMPLING_RATE = 1.0D0
   INTEGER :: N_SPECTRAL_BINS = 200
   REAL(KIND=8) :: SPECTRAL_VLOS_MIN = -2.3D6  ! m/s (corresponds to ~-5 nm from Hα center)
   REAL(KIND=8) :: SPECTRAL_VLOS_MAX = +2.3D6  ! m/s (corresponds to ~+5 nm from Hα center)
   REAL(KIND=8), DIMENSION(3) :: SPECTRAL_LOS_DIRECTION = (/0.d0, 0.d0, 1.d0/)  ! Default: Z direction
   REAL(KIND=8) :: SPECTRAL_BACKGROUND_FRACTION = 0.0D0  ! (unused, kept for compatibility)
   
   ! Per-reaction spectral histogram: (reaction_id, wavelength_bin)
   INTEGER(KIND=8), DIMENSION(:,:), ALLOCATABLE :: SPECTRAL_HISTOGRAM  ! (N_REACTIONS, N_SPECTRAL_BINS)
   INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: SPECTRAL_TOTAL_EVENTS ! Per-reaction event counts
   INTEGER :: N_HALPHA_REACTIONS = 0  ! Number of reactions that produce Hα
   
   ! ========== Hα Passive Diagnostic Event Buffer (Phase 2) ==========
   INTEGER, PARAMETER :: HA_BUFFER_SIZE = 10000
   INTEGER :: HA_EVENT_COUNT = 0
   LOGICAL :: BOOL_HA_PASSIVE_DIAGNOSTIC = .FALSE.
   
   TYPE HA_EVENT_RECORD
      INTEGER :: TIMESTEP
      REAL(KIND=8) :: TIME
      INTEGER :: REACTION_ID
      INTEGER :: PROJECTILE_SPECIES
      REAL(KIND=8) :: POS_X, POS_Y, POS_Z
      REAL(KIND=8) :: VEL_X, VEL_Y, VEL_Z
      REAL(KIND=8) :: MACRO_WEIGHT
      REAL(KIND=8) :: COLLISION_ENERGY  ! eV
   END TYPE HA_EVENT_RECORD
   
   TYPE(HA_EVENT_RECORD), DIMENSION(HA_BUFFER_SIZE) :: HA_EVENT_BUFFER
   ! ==================================================================
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Reaction Statistics for VTK Output !!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   LOGICAL :: BOOL_REACTION_STATISTICS = .FALSE.
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: REACTION_CELL_COUNTS            ! Window counts since last dump
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: REACTION_CELL_COUNTS_GLOBAL     ! Reduced window counts for output
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: REACTION_CELL_COUNTS_CUM        ! Cumulative counts since start
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: REACTION_CELL_COUNTS_CUM_GLOBAL ! Reduced cumulative counts for output

   INTEGER           :: BGK_MODEL_TYPE_INT = 0
   REAL(KIND=8)      :: BGK_BG_DENS, BGK_SIGMA, BGK_BG_MASS


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Particles injection from boundaries !!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   LOGICAL      :: BOOL_BOUNDINJECT = .FALSE. ! Assign default value!
   LOGICAL      :: BOOL_INJ_XMIN=.FALSE.,BOOL_INJ_XMAX=.FALSE. ! Assign default value!
   LOGICAL      :: BOOL_INJ_YMIN=.FALSE.,BOOL_INJ_YMAX=.FALSE. ! Assign default value!
   REAL(KIND=8) :: NRHO_BOUNDINJECT
   REAL(KIND=8) :: UX_BOUND, UY_BOUND, UZ_BOUND
   REAL(KIND=8) :: TTRA_BOUND, TROT_BOUND, TVIB_BOUND ! No TTRAX_ TTRAY_ TTRAZ_ for now! IMPLEMENT IT!
   INTEGER      :: MIX_BOUNDINJECT

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: nfs_XMIN, nfs_XMAX, nfs_YMIN, nfs_YMAX
   REAL(KIND=8) :: KAPPA_XMIN, KAPPA_XMAX, KAPPA_YMIN, KAPPA_YMAX
   REAL(KIND=8) :: S_NORM_XMIN, S_NORM_XMAX, S_NORM_YMIN, S_NORM_YMAX
   INTEGER           :: REMOVE_MIX = -1

   ! Energy-based particle removal
   LOGICAL           :: BOOL_REMOVE_LOW_ENERGY = .FALSE.
   INTEGER           :: N_ENERGY_REMOVAL_RULES = 0

   TYPE ENERGY_REMOVAL_RULE
      INTEGER :: N_SPECIES                              ! Number of species in this rule
      INTEGER, DIMENSION(:), ALLOCATABLE :: SPECIES_IDS ! Species IDs to check
      REAL(KIND=8) :: ENERGY_THRESHOLD_EV               ! Threshold in eV
   END TYPE ENERGY_REMOVAL_RULE

   TYPE(ENERGY_REMOVAL_RULE), DIMENSION(:), ALLOCATABLE :: ENERGY_REMOVAL_RULES


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Particles injection from line source !!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Not used anymore
   LOGICAL      :: BOOL_LINESOURCE = .FALSE. ! Assign default value!
  
   REAL(KIND=8) :: X_LINESOURCE, Y_LINESOURCE, L_LINESOURCE
   REAL(KIND=8) :: NRHO_LINESOURCE
   REAL(KIND=8) :: UX_LINESOURCE, UY_LINESOURCE, UZ_LINESOURCE
   REAL(KIND=8) :: TTRA_LINESOURCE, TROT_LINESOURCE, TVIB_LINESOURCE

   INTEGER      :: nfs_LINESOURCE
   REAL(KIND=8) :: KAPPA_LINESOURCE
   REAL(KIND=8) :: S_NORM_LINESOURCE

   ! This is used instead.
   INTEGER         :: N_LINESOURCES = 0

   TYPE LINESOURCE
      REAL(KIND=8) :: X1, Y1, X2, Y2
      REAL(KIND=8) :: NRHO
      REAL(KIND=8) :: UX, UY, UZ
      REAL(KIND=8) :: TTRA, TROT, TVIB
      INTEGER      :: MIX_ID

      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: nfs
      REAL(KIND=8) :: KAPPA
      REAL(KIND=8) :: S_NORM
      REAL(KIND=8) :: NORMX, NORMY
   END TYPE LINESOURCE

   TYPE(LINESOURCE), DIMENSION(:), ALLOCATABLE :: LINESOURCES

   CHARACTER*256 :: FLUXDUMP_SAVE_PATH
   LOGICAL :: BOOL_DUMP_FLUXES = .FALSE.

   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Walls !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   INTEGER         :: N_WALLS = 0

   TYPE WALL
      REAL(KIND=8) :: X1, Y1, X2, Y2
      REAL(KIND=8) :: TEMP
      LOGICAL      :: SPECULAR, DIFFUSE, POROUS, REACT
      REAL(KIND=8) :: TRANSMISSIVITY
      REAL(KIND=8) :: NORMX, NORMY
   END TYPE WALL

   TYPE(WALL), DIMENSION(:), ALLOCATABLE :: WALLS


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Multispecies !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   CHARACTER(LEN=256) :: SPECIES_FILENAME
   INTEGER            :: N_SPECIES = 0
   
   TYPE SPECIES_DATA_STRUCTURE
      CHARACTER*10 :: NAME
      REAL(KIND=8) :: MOLWT
      REAL(KIND=8) :: MOLECULAR_MASS
      INTEGER      :: ROTDOF
      REAL(KIND=8) :: ROTREL
      INTEGER      :: VIBDOF
      REAL(KIND=8) :: VIBREL
      REAL(KIND=8) :: VIBTEMP
      REAL(KIND=8) :: SPWT
      REAL(KIND=8) :: INVSPWT
      REAL(KIND=8) :: CHARGE
      REAL(KIND=8) :: DIAM
      REAL(KIND=8) :: OMEGA
      REAL(KIND=8) :: TREF
      REAL(KIND=8) :: ALPHA
      REAL(KIND=8) :: SIGMA
      !REAL(KIND=8) :: NU
      REAL(KIND=8) :: CREF
      
   END TYPE SPECIES_DATA_STRUCTURE

   TYPE(SPECIES_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: SPECIES, TEMP_SPECIES


   INTEGER            :: N_MIXTURES = 0

   TYPE MIXTURE_COMPONENT
      INTEGER :: ID
      CHARACTER*64 :: NAME
      REAL(KIND=8) :: MOLFRAC
   END TYPE MIXTURE_COMPONENT

   TYPE MIXTURE
      CHARACTER*64 :: NAME
      INTEGER      :: N_COMPONENTS
      TYPE(MIXTURE_COMPONENT), DIMENSION(:), ALLOCATABLE :: COMPONENTS
   END TYPE MIXTURE

   TYPE(MIXTURE), DIMENSION(:), ALLOCATABLE :: MIXTURES
 


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Chemical reactions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   CHARACTER(LEN=256) :: REACTIONS_FILENAME
   INTEGER :: N_REACTIONS = 0

   ENUM, BIND(C)
      ENUMERATOR FIXED_RATE, TCE, LXCAT
   END ENUM
   
   TYPE REACTIONS_DATA_STRUCTURE
      INTEGER(KIND(FIXED_RATE)) :: TYPE = FIXED_RATE
      INTEGER :: R1_SP_ID
      INTEGER :: R2_SP_ID
      INTEGER :: P1_SP_ID
      INTEGER :: P2_SP_ID
      INTEGER :: P3_SP_ID
      INTEGER :: P4_SP_ID
      REAL(KIND=8) :: A, N, EA
      REAL(KIND=8) :: C1, C2, C3
      INTEGER :: N_PROD
      LOGICAL :: IS_CEX
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TABLE_ENERGY
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TABLE_CS
      REAL(KIND=8) :: MAX_SIGMA
      INTEGER :: COUNTS
      ! Spectral diagnostics fields
      LOGICAL :: PRODUCES_HALPHA = .FALSE.
      INTEGER :: EMITTING_PRODUCT_ID = 0
      ! Passive Hα diagnostic (Phase 1)
      LOGICAL :: IS_HA_HEAVY = .FALSE.  ! Heavy-particle Ha(total) channel flag
   END TYPE REACTIONS_DATA_STRUCTURE

   TYPE(REACTIONS_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: REACTIONS, TEMP_REACTIONS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Wall reactions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   CHARACTER(LEN=256) :: WALL_REACTIONS_FILENAME
   INTEGER :: N_WALL_REACTIONS = 0
   
   TYPE WALL_REACTIONS_DATA_STRUCTURE
      INTEGER :: R_SP_ID
      INTEGER :: P1_SP_ID
      INTEGER :: P2_SP_ID
      REAL(KIND=8) :: PROB
      INTEGER :: N_PROD
   END TYPE WALL_REACTIONS_DATA_STRUCTURE

   TYPE(WALL_REACTIONS_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: WALL_REACTIONS, TEMP_WALL_REACTIONS


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Secondary Electron Emission (SEE) !!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   LOGICAL :: BOOL_SEE_ENABLED = .FALSE.    ! Enable secondary electron emission
   INTEGER :: SEE_ELECTRON_SPECIES_ID = -1  ! Species ID for electrons
   
   ! Multi-species support for ion and neutral SEE
   INTEGER, ALLOCATABLE, DIMENSION(:) :: SEE_ION_SPECIES_IDS      ! All positive ion species IDs
   INTEGER, ALLOCATABLE, DIMENSION(:) :: SEE_NEUTRAL_SPECIES_IDS  ! All neutral atom species IDs
   INTEGER :: N_SEE_ION_SPECIES = 0                               ! Number of ion species
   INTEGER :: N_SEE_NEUTRAL_SPECIES = 0                           ! Number of neutral species
   
   INTEGER :: N_SEE_MATERIALS = 0           ! Number of SEE materials defined
   
   ! SEE material properties structure
   TYPE SEE_MATERIAL_PROPERTIES
      CHARACTER*64 :: NAME              ! Material name
      ! Electron SEE parameters (Vaughan model)
      REAL(KIND=8) :: DELTA_MAX         ! Maximum electron SEE yield
      REAL(KIND=8) :: E_MAX             ! Energy at maximum electron yield [eV]
      REAL(KIND=8) :: E_TH              ! Electron threshold energy [eV]
      REAL(KIND=8) :: S_PARAMETER       ! Shape parameter (typically ~1.35)
      ! Ion SEE parameters (Hagstrum model)
      REAL(KIND=8) :: GAMMA_MAX         ! Maximum ion SEE yield
      REAL(KIND=8) :: ION_E_THRESHOLD   ! Ion threshold energy [eV]
      REAL(KIND=8) :: ION_ALPHA         ! Ion energy exponent
      REAL(KIND=8) :: ION_BETA          ! Ion decay parameter
      REAL(KIND=8) :: ION_CHARGE_FACTOR ! Ion charge state correction
      ! Neutral atom SEE parameters (pure kinetic emission)
      REAL(KIND=8) :: NEUTRAL_GAMMA_MAX       ! Maximum neutral SEE yield
      REAL(KIND=8) :: NEUTRAL_E_THRESHOLD     ! Neutral threshold energy [eV]
      REAL(KIND=8) :: NEUTRAL_ALPHA           ! Neutral energy exponent
      REAL(KIND=8) :: NEUTRAL_BETA            ! Neutral decay parameter
      ! Common parameters
      REAL(KIND=8) :: W_WORK_FUNCTION   ! Work function [eV]
      REAL(KIND=8) :: P1                ! Backscatter parameter 1
      REAL(KIND=8) :: P2                ! Backscatter parameter 2
      REAL(KIND=8) :: E1                ! Backscatter energy parameter 1 [eV]
      REAL(KIND=8) :: E2                ! Backscatter energy parameter 2 [eV]
      REAL(KIND=8) :: SIGMA             ! Surface roughness parameter
   END TYPE SEE_MATERIAL_PROPERTIES

   TYPE(SEE_MATERIAL_PROPERTIES), DIMENSION(:), ALLOCATABLE :: SEE_MATERIALS

   ! SEE boundary mapping
   TYPE SEE_BOUNDARY_MAPPING
      CHARACTER*64 :: BOUNDARY_NAME  ! Boundary group name (e.g., "anode", "cathode")
      INTEGER :: BOUNDARY_ID          ! Boundary group ID (for backward compatibility)
      INTEGER :: WALL_ID              ! Wall ID (-1 for grid boundaries)
      INTEGER :: MATERIAL_ID          ! SEE material ID
      LOGICAL :: ENABLED              ! Is SEE enabled for this boundary
      LOGICAL :: USE_NAME             ! Whether to use name or ID for boundary matching
   END TYPE SEE_BOUNDARY_MAPPING

   TYPE(SEE_BOUNDARY_MAPPING), DIMENSION(:), ALLOCATABLE :: SEE_BOUNDARY_MAP
   INTEGER :: N_SEE_BOUNDARIES = 0

   ! SEE statistics
   INTEGER(KIND=8) :: SEE_TOTAL_IMPACTS = 0        ! Total electron impacts
   INTEGER(KIND=8) :: SEE_TOTAL_EMISSIONS = 0      ! Total secondary emissions
   INTEGER(KIND=8) :: SEE_TOTAL_ION_IMPACTS = 0    ! Total ion impacts
   INTEGER(KIND=8) :: SEE_TOTAL_ION_EMISSIONS = 0  ! Total ion-induced emissions
   INTEGER(KIND=8) :: SEE_TOTAL_NEUTRAL_IMPACTS = 0    ! Total neutral atom impacts
   INTEGER(KIND=8) :: SEE_TOTAL_NEUTRAL_EMISSIONS = 0  ! Total neutral-induced emissions
   REAL(KIND=8) :: SEE_TOTAL_YIELD = 0.d0          ! Average electron yield
   REAL(KIND=8) :: SEE_TOTAL_ION_YIELD = 0.d0      ! Average ion yield (γ coefficient)
   REAL(KIND=8) :: SEE_TOTAL_NEUTRAL_YIELD = 0.d0  ! Average neutral yield
   INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: SEE_MATERIAL_IMPACTS    ! Per-material impacts
   INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: SEE_MATERIAL_EMISSIONS  ! Per-material emissions

   CHARACTER*256 :: SEE_MATERIALS_FILENAME = ''    ! SEE materials definition file
   CHARACTER*256 :: SEE_STATS_SAVE_PATH = ''       ! SEE statistics output path


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Average flowfield !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   CHARACTER*256                           :: FLOWFIELD_SAVE_PATH

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_N

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_NP

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_VX
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_VY
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_VZ

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_TTRX
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_TTRY
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_TTRZ
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_TTR

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_TROT
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_TVIB

   INTEGER                                 :: AVG_CUMULATED
   INTEGER, DIMENSION(:), ALLOCATABLE      :: AVG_CUMULATED_INTENSIVE_ONE
   INTEGER, DIMENSION(:), ALLOCATABLE      :: AVG_CUMULATED_INTENSIVE_TWO

   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: AVG_MOMENTS
   LOGICAL                                   :: BOOL_DUMP_MOMENTS = .FALSE.

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_PHI

   INTEGER                                 :: BOUNDARY_AVG_CUMULATED

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Average boundary !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   CHARACTER*256                           :: BOUNDARY_SAVE_PATH

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_NIN
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_NOUT

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_PXIN
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_PYIN
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_PZIN
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_PXOUT
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_PYOUT
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_PZOUT

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_PXEM
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_PYEM
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_PZEM

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_EIN
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_EOUT

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_PHI_BOUND
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_QRHO_BOUND

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: TIMESTEP_PFLUID_BOUND

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_NIN
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_NOUT

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_PXIN
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_PYIN
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_PZIN
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_PXOUT
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_PYOUT
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_PZOUT

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_PXEM
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_PYEM
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_PZEM

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_EIN
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_EOUT

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_PHI_BOUND
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_QRHO_BOUND

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: AVG_PFLUID_BOUND

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!! Timers !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   REAL(KIND=8), DIMENSION(6) :: TIMERS_START_TIME = 0.d0
   REAL(KIND=8), DIMENSION(6) :: TIMERS_ELAPSED = 0.d0
   

   REAL(KIND=8) :: FIELD_POWER
   REAL(KIND=8) :: COIL_CURRENT = 1.5d0
   REAL(KIND=8) :: FIELD_POWER_TARGET = 20.d0
   INTEGER :: FIELD_POWER_NUMAVG = 737
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: FIELD_POWER_AVG


CONTAINS  ! @@@@@@@@@@@@@@@@@@@@@ SUBROUTINES @@@@@@@@@@@@@@@@@@@@@@@@

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE NEWTYPE -> defines a new type needed by MPI             !
   ! to package messages in the "particle" format                       !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   SUBROUTINE NEWTYPE
   
      INTEGER :: ii, extent_dpr, extent_int, extent_int8, extent_logical
      INTEGER, DIMENSION(15) :: blocklengths, oldtypes, offsets
     
      CALL MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, extent_dpr,  ierr)  
      CALL MPI_TYPE_EXTENT(MPI_INTEGER,          extent_int,  ierr)
      CALL MPI_TYPE_EXTENT(MPI_INTEGER8,         extent_int8, ierr)
      CALL MPI_TYPE_EXTENT(MPI_LOGICAL,         extent_logical, ierr)
        
      blocklengths = 1
     
      oldtypes(1:9) = MPI_DOUBLE_PRECISION  
      oldtypes(10:11) = MPI_INTEGER
      oldtypes(12) = MPI_INTEGER8
      oldtypes(13) = MPI_LOGICAL
      oldtypes(14) = MPI_LOGICAL
      oldtypes(15) = MPI_DOUBLE_PRECISION
          
      offsets(1) = 0  
      DO ii = 2, 10
         offsets(ii) = offsets(ii - 1) + extent_dpr * blocklengths(ii - 1)
      END DO
      offsets(11) = offsets(10) + extent_int * blocklengths(10)
      offsets(12) = offsets(11) + extent_int * blocklengths(11)
      offsets(13) = offsets(12) + extent_int8 * blocklengths(12)
      offsets(14) = offsets(13) + extent_logical * blocklengths(13)
      offsets(15) = offsets(14) + extent_logical * blocklengths(14)
      
      CALL MPI_TYPE_STRUCT(15, blocklengths, offsets, oldtypes, MPI_PARTICLE_DATA_STRUCTURE, ierr)  
      CALL MPI_TYPE_COMMIT(MPI_PARTICLE_DATA_STRUCTURE, ierr)   
   
   END SUBROUTINE NEWTYPE


END MODULE global
