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
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

! This module contains secondary electron emission (SEE) functionality

MODULE secondary_electron_emission

   USE global
   USE particle
   USE mpi_common
   USE tools
   
   IMPLICIT NONE

   TYPE SEE_ION_TABLE_DEFINITION
      CHARACTER(LEN=64) :: MATERIAL_REF = ''
      CHARACTER(LEN=64) :: SPECIES_NAME = ''
      CHARACTER(LEN=256) :: FILENAME = ''
   END TYPE SEE_ION_TABLE_DEFINITION

   TYPE SEE_ION_YIELD_TABLE
      CHARACTER(LEN=64) :: MATERIAL_REF = ''
      CHARACTER(LEN=64) :: SPECIES_NAME = ''
      CHARACTER(LEN=256) :: SOURCE_FILENAME = ''
      INTEGER :: MATERIAL_ID = -1
      INTEGER :: SPECIES_ID = -1
      INTEGER :: N_ENERGIES = 0
      INTEGER :: N_ANGLES = 0
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ENERGY_GRID
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ANGLE_GRID
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: YIELD_GRID
   END TYPE SEE_ION_YIELD_TABLE

   TYPE(SEE_ION_TABLE_DEFINITION), DIMENSION(:), ALLOCATABLE :: SEE_ION_TABLE_DEFS
   TYPE(SEE_ION_YIELD_TABLE), DIMENSION(:), ALLOCATABLE :: SEE_ION_YIELD_TABLES
   INTEGER :: N_SEE_ION_TABLE_DEFS = 0
   INTEGER :: N_SEE_ION_YIELD_TABLES = 0
   
   CONTAINS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE INIT_SEE -> Initialize secondary electron emission system !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE INIT_SEE

      IMPLICIT NONE
      
      INTEGER :: IM

      IF (.NOT. BOOL_SEE_ENABLED) RETURN

      ! Validate species array
      IF (.NOT. ALLOCATED(SPECIES) .OR. N_SPECIES <= 0) THEN
         CALL ERROR_ABORT('SEE enabled but species array not initialized!')
      END IF

      ! Find electron species ID
      SEE_ELECTRON_SPECIES_ID = -1
      DO IM = 1, N_SPECIES
         IF (ABS(SPECIES(IM)%CHARGE + 1.d0) < 1.d-6) THEN  ! Electron has charge -1
            SEE_ELECTRON_SPECIES_ID = IM
            EXIT
         END IF
      END DO

      IF (SEE_ELECTRON_SPECIES_ID == -1) THEN
         CALL ERROR_ABORT('SEE enabled but no electron species found with charge -1!')
      END IF
      
      ! Validate electron species ID range
      IF (SEE_ELECTRON_SPECIES_ID < 1 .OR. SEE_ELECTRON_SPECIES_ID > N_SPECIES) THEN
         CALL ERROR_ABORT('SEE electron species ID out of valid range!')
      END IF

      ! Find all positive ion species (for glow discharge)
      N_SEE_ION_SPECIES = 0
      DO IM = 1, N_SPECIES
         IF (SPECIES(IM)%CHARGE > 0.5d0) THEN  ! Positive ion
            N_SEE_ION_SPECIES = N_SEE_ION_SPECIES + 1
         END IF
      END DO

      IF (N_SEE_ION_SPECIES > 0) THEN
         ALLOCATE(SEE_ION_SPECIES_IDS(N_SEE_ION_SPECIES))
         N_SEE_ION_SPECIES = 0  ! Reset counter for actual assignment
         DO IM = 1, N_SPECIES
            IF (SPECIES(IM)%CHARGE > 0.5d0) THEN
               N_SEE_ION_SPECIES = N_SEE_ION_SPECIES + 1
               SEE_ION_SPECIES_IDS(N_SEE_ION_SPECIES) = IM
            END IF
         END DO
         
         IF (PROC_ID == 0) THEN
            WRITE(*,*) 'Ion SEE enabled for', N_SEE_ION_SPECIES, 'species:'
            DO IM = 1, N_SEE_ION_SPECIES
               WRITE(*,'(A,I3,A,A)') '  - Species ID ', SEE_ION_SPECIES_IDS(IM), ': ', &
                                     TRIM(SPECIES(SEE_ION_SPECIES_IDS(IM))%NAME)
            END DO
         END IF
      ELSE
         WRITE(*,*) 'Warning: No positive ion species found for ion SEE!'
         WRITE(*,*) 'Only electron SEE will be active.'
      END IF

      ! Find all neutral atom species
      N_SEE_NEUTRAL_SPECIES = 0
      DO IM = 1, N_SPECIES
         IF (ABS(SPECIES(IM)%CHARGE) < 1.d-6) THEN  ! Neutral atom (charge = 0)
            N_SEE_NEUTRAL_SPECIES = N_SEE_NEUTRAL_SPECIES + 1
         END IF
      END DO

      IF (N_SEE_NEUTRAL_SPECIES > 0) THEN
         ALLOCATE(SEE_NEUTRAL_SPECIES_IDS(N_SEE_NEUTRAL_SPECIES))
         N_SEE_NEUTRAL_SPECIES = 0  ! Reset counter for actual assignment
         DO IM = 1, N_SPECIES
            IF (ABS(SPECIES(IM)%CHARGE) < 1.d-6) THEN
               N_SEE_NEUTRAL_SPECIES = N_SEE_NEUTRAL_SPECIES + 1
               SEE_NEUTRAL_SPECIES_IDS(N_SEE_NEUTRAL_SPECIES) = IM
            END IF
         END DO
         
         IF (PROC_ID == 0) THEN
            WRITE(*,*) 'Neutral atom SEE enabled for', N_SEE_NEUTRAL_SPECIES, 'species:'
            DO IM = 1, N_SEE_NEUTRAL_SPECIES
               WRITE(*,'(A,I3,A,A)') '  - Species ID ', SEE_NEUTRAL_SPECIES_IDS(IM), ': ', &
                                     TRIM(SPECIES(SEE_NEUTRAL_SPECIES_IDS(IM))%NAME)
            END DO
         END IF
      ELSE
         WRITE(*,*) 'Warning: No neutral species found for neutral atom SEE!'
         WRITE(*,*) 'Only electron and ion SEE will be active.'
      END IF

      CALL LOAD_SEE_ION_YIELD_TABLES

      ! Allocate statistics arrays
      IF (ALLOCATED(SEE_MATERIAL_IMPACTS)) DEALLOCATE(SEE_MATERIAL_IMPACTS)
      IF (ALLOCATED(SEE_MATERIAL_EMISSIONS)) DEALLOCATE(SEE_MATERIAL_EMISSIONS)
      IF (N_SEE_MATERIALS > 0) THEN
         ALLOCATE(SEE_MATERIAL_IMPACTS(N_SEE_MATERIALS))
         ALLOCATE(SEE_MATERIAL_EMISSIONS(N_SEE_MATERIALS))
         SEE_MATERIAL_IMPACTS = 0
         SEE_MATERIAL_EMISSIONS = 0
      END IF

      IF (ALLOCATED(SEE_SPECIES_IMPACTS)) DEALLOCATE(SEE_SPECIES_IMPACTS)
      IF (ALLOCATED(SEE_SPECIES_EMISSIONS)) DEALLOCATE(SEE_SPECIES_EMISSIONS)
      IF (N_SPECIES > 0) THEN
         ALLOCATE(SEE_SPECIES_IMPACTS(N_SPECIES))
         ALLOCATE(SEE_SPECIES_EMISSIONS(N_SPECIES))
         SEE_SPECIES_IMPACTS = 0
         SEE_SPECIES_EMISSIONS = 0
      END IF

      ! Initialize global statistics
      SEE_TOTAL_IMPACTS = 0
      SEE_TOTAL_EMISSIONS = 0
      SEE_TOTAL_YIELD = 0.d0
      SEE_TOTAL_ION_IMPACTS = 0
      SEE_TOTAL_ION_EMISSIONS = 0
      SEE_TOTAL_ION_YIELD = 0.d0
      SEE_TOTAL_NEUTRAL_IMPACTS = 0
      SEE_TOTAL_NEUTRAL_EMISSIONS = 0
      SEE_TOTAL_NEUTRAL_YIELD = 0.d0

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '> SEE initialized with ', N_SEE_MATERIALS, ' materials'
         WRITE(*,*) '> Electron species ID: ', SEE_ELECTRON_SPECIES_ID
      END IF

   END SUBROUTINE INIT_SEE

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE READ_SEE_MATERIALS -> Read SEE material properties       !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE READ_SEE_MATERIALS(FILENAME)

      IMPLICIT NONE
      
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      INTEGER :: in_unit = 20
      INTEGER :: ios, IM, N_STR
      CHARACTER(LEN=512) :: line
      CHARACTER(LEN=64) :: name
      REAL(KIND=8) :: delta_max, e_max, e_th, s_param, w_func
      REAL(KIND=8) :: gamma_max, ion_e_th, ion_alpha, ion_beta, ion_charge_factor
      REAL(KIND=8) :: neutral_gamma_max, neutral_e_th, neutral_alpha, neutral_beta
      REAL(KIND=8) :: p1, p2, e1, e2, sigma
      REAL(KIND=8) :: ion_kin_angle_exp, ion_kin_angle_max, ion_pot_angle_exp, ion_pot_angle_max
      REAL(KIND=8) :: neutral_angle_exp, neutral_angle_max
      CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE :: STRARRAY
      TYPE(SEE_MATERIAL_PROPERTIES), DIMENSION(:), ALLOCATABLE :: TEMP_MATERIALS

      IF (TRIM(FILENAME) == '') RETURN

      OPEN(UNIT=in_unit, FILE=TRIM(FILENAME), STATUS='old', IOSTAT=ios)
      IF (ios /= 0) THEN
         CALL ERROR_ABORT('Cannot open SEE materials file: ' // TRIM(FILENAME))
      END IF

      N_SEE_MATERIALS = 0

      ! Read materials
      DO WHILE (.TRUE.)
         READ(in_unit, '(A)', iostat=ios) line
         IF (ios /= 0) EXIT
         
         line = ADJUSTL(line)
         IF (line(1:1) == '#' .OR. LEN_TRIM(line) == 0) CYCLE  ! Skip comments and empty lines

         CALL SPLIT_STR(TRIM(line), ' ', STRARRAY, N_STR)
         IF (N_STR < 20) THEN
            WRITE(*,*) 'Error parsing SEE material line: ', TRIM(line)
            WRITE(*,*) 'Expected at least 20 tokens, got: ', N_STR
            IF (ALLOCATED(STRARRAY)) DEALLOCATE(STRARRAY)
            CYCLE
         END IF

         ios = 0
         READ(STRARRAY(1), *, IOSTAT=ios) name
         IF (ios == 0) READ(STRARRAY(2), *, IOSTAT=ios) delta_max
         IF (ios == 0) READ(STRARRAY(3), *, IOSTAT=ios) e_max
         IF (ios == 0) READ(STRARRAY(4), *, IOSTAT=ios) e_th
         IF (ios == 0) READ(STRARRAY(5), *, IOSTAT=ios) s_param
         IF (ios == 0) READ(STRARRAY(6), *, IOSTAT=ios) gamma_max
         IF (ios == 0) READ(STRARRAY(7), *, IOSTAT=ios) ion_e_th
         IF (ios == 0) READ(STRARRAY(8), *, IOSTAT=ios) ion_alpha
         IF (ios == 0) READ(STRARRAY(9), *, IOSTAT=ios) ion_beta
         IF (ios == 0) READ(STRARRAY(10), *, IOSTAT=ios) ion_charge_factor
         IF (ios == 0) READ(STRARRAY(11), *, IOSTAT=ios) w_func
         IF (ios == 0) READ(STRARRAY(12), *, IOSTAT=ios) p1
         IF (ios == 0) READ(STRARRAY(13), *, IOSTAT=ios) p2
         IF (ios == 0) READ(STRARRAY(14), *, IOSTAT=ios) e1
         IF (ios == 0) READ(STRARRAY(15), *, IOSTAT=ios) e2
         IF (ios == 0) READ(STRARRAY(16), *, IOSTAT=ios) sigma
         IF (ios == 0) READ(STRARRAY(17), *, IOSTAT=ios) neutral_gamma_max
         IF (ios == 0) READ(STRARRAY(18), *, IOSTAT=ios) neutral_e_th
         IF (ios == 0) READ(STRARRAY(19), *, IOSTAT=ios) neutral_alpha
         IF (ios == 0) READ(STRARRAY(20), *, IOSTAT=ios) neutral_beta
         IF (ios /= 0) THEN
            WRITE(*,*) 'Error parsing SEE material line: ', TRIM(line)
            WRITE(*,*) 'Malformed numeric token, error code: ', ios
            IF (ALLOCATED(STRARRAY)) DEALLOCATE(STRARRAY)
            CYCLE
         END IF

         ! Optional material-specific angle parameters. Defaults reproduce the
         ! previously hard-coded behavior.
         ion_kin_angle_exp = 0.25d0
         ion_kin_angle_max = 3.0d0
         ion_pot_angle_exp = 0.30d0
         ion_pot_angle_max = 2.0d0
         neutral_angle_exp = 0.15d0
         neutral_angle_max = 3.0d0

         IF (N_STR >= 26) THEN
            READ(STRARRAY(21), *, IOSTAT=ios) ion_kin_angle_exp
            IF (ios == 0) READ(STRARRAY(22), *, IOSTAT=ios) ion_kin_angle_max
            IF (ios == 0) READ(STRARRAY(23), *, IOSTAT=ios) ion_pot_angle_exp
            IF (ios == 0) READ(STRARRAY(24), *, IOSTAT=ios) ion_pot_angle_max
            IF (ios == 0) READ(STRARRAY(25), *, IOSTAT=ios) neutral_angle_exp
            IF (ios == 0) READ(STRARRAY(26), *, IOSTAT=ios) neutral_angle_max
            IF (ios /= 0) THEN
               WRITE(*,*) 'Error parsing optional SEE angle parameters: ', TRIM(line)
               WRITE(*,*) 'Malformed optional token, error code: ', ios
               IF (ALLOCATED(STRARRAY)) DEALLOCATE(STRARRAY)
               CYCLE
            END IF
         END IF

         ! Add material to array
         IF (ALLOCATED(SEE_MATERIALS)) THEN
            ALLOCATE(TEMP_MATERIALS(N_SEE_MATERIALS+1))
            TEMP_MATERIALS(1:N_SEE_MATERIALS) = SEE_MATERIALS(1:N_SEE_MATERIALS)
            CALL MOVE_ALLOC(TEMP_MATERIALS, SEE_MATERIALS)
         ELSE
            ALLOCATE(SEE_MATERIALS(1))
         END IF

         N_SEE_MATERIALS = N_SEE_MATERIALS + 1
         
         SEE_MATERIALS(N_SEE_MATERIALS)%NAME = TRIM(name)
         ! Electron SEE parameters
         SEE_MATERIALS(N_SEE_MATERIALS)%DELTA_MAX = delta_max
         SEE_MATERIALS(N_SEE_MATERIALS)%E_MAX = e_max
         SEE_MATERIALS(N_SEE_MATERIALS)%E_TH = e_th
         SEE_MATERIALS(N_SEE_MATERIALS)%S_PARAMETER = s_param
         ! Ion SEE parameters
         SEE_MATERIALS(N_SEE_MATERIALS)%GAMMA_MAX = gamma_max
         SEE_MATERIALS(N_SEE_MATERIALS)%ION_E_THRESHOLD = ion_e_th
         SEE_MATERIALS(N_SEE_MATERIALS)%ION_ALPHA = ion_alpha
         SEE_MATERIALS(N_SEE_MATERIALS)%ION_BETA = ion_beta
         SEE_MATERIALS(N_SEE_MATERIALS)%ION_CHARGE_FACTOR = ion_charge_factor
         ! Neutral atom SEE parameters
         SEE_MATERIALS(N_SEE_MATERIALS)%NEUTRAL_GAMMA_MAX = neutral_gamma_max
         SEE_MATERIALS(N_SEE_MATERIALS)%NEUTRAL_E_THRESHOLD = neutral_e_th
         SEE_MATERIALS(N_SEE_MATERIALS)%NEUTRAL_ALPHA = neutral_alpha
         SEE_MATERIALS(N_SEE_MATERIALS)%NEUTRAL_BETA = neutral_beta
         SEE_MATERIALS(N_SEE_MATERIALS)%ION_KIN_ANGLE_EXP = ion_kin_angle_exp
         SEE_MATERIALS(N_SEE_MATERIALS)%ION_KIN_ANGLE_MAX = ion_kin_angle_max
         SEE_MATERIALS(N_SEE_MATERIALS)%ION_POT_ANGLE_EXP = ion_pot_angle_exp
         SEE_MATERIALS(N_SEE_MATERIALS)%ION_POT_ANGLE_MAX = ion_pot_angle_max
         SEE_MATERIALS(N_SEE_MATERIALS)%NEUTRAL_ANGLE_EXP = neutral_angle_exp
         SEE_MATERIALS(N_SEE_MATERIALS)%NEUTRAL_ANGLE_MAX = neutral_angle_max
         ! Common parameters
         SEE_MATERIALS(N_SEE_MATERIALS)%W_WORK_FUNCTION = w_func
         SEE_MATERIALS(N_SEE_MATERIALS)%P1 = p1
         SEE_MATERIALS(N_SEE_MATERIALS)%P2 = p2
         SEE_MATERIALS(N_SEE_MATERIALS)%E1 = e1
         SEE_MATERIALS(N_SEE_MATERIALS)%E2 = e2
         SEE_MATERIALS(N_SEE_MATERIALS)%SIGMA = sigma

         IF (ALLOCATED(STRARRAY)) DEALLOCATE(STRARRAY)

      END DO

      CLOSE(in_unit)

            IF (PROC_ID == 0 .AND. N_SEE_MATERIALS > 0) THEN
         WRITE(*,*) '> Read ', N_SEE_MATERIALS, ' SEE materials:'
         DO IM = 1, N_SEE_MATERIALS
            WRITE(*,'(A,I3,A,A,A,F6.3,A,F6.1,A,F6.3,A)') '  Material ', IM, ': ', &
                  TRIM(SEE_MATERIALS(IM)%NAME), ' (δ_max=', SEE_MATERIALS(IM)%DELTA_MAX, &
                  ', γ_max=', SEE_MATERIALS(IM)%GAMMA_MAX, ', E_max=', SEE_MATERIALS(IM)%E_MAX, ' eV)'
         END DO
         WRITE(*,*) '> Material ID range: 1 to ', N_SEE_MATERIALS
      END IF

   END SUBROUTINE READ_SEE_MATERIALS


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE DEF_SEE_ION_YIELD_TABLE -> Register one ion SEE table    !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE DEF_SEE_ION_YIELD_TABLE(DEFINITION)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: DEFINITION
      CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE :: STRARRAY
      TYPE(SEE_ION_TABLE_DEFINITION), DIMENSION(:), ALLOCATABLE :: TEMP_DEFS
      INTEGER :: N_STR

      CALL SPLIT_STR(TRIM(DEFINITION), ' ', STRARRAY, N_STR)
      IF (N_STR < 3) THEN
         CALL ERROR_ABORT('SEE ion yield table requires: material_ref ion_species table_file')
      END IF

      IF (ALLOCATED(SEE_ION_TABLE_DEFS)) THEN
         ALLOCATE(TEMP_DEFS(N_SEE_ION_TABLE_DEFS + 1))
         TEMP_DEFS(1:N_SEE_ION_TABLE_DEFS) = SEE_ION_TABLE_DEFS(1:N_SEE_ION_TABLE_DEFS)
         CALL MOVE_ALLOC(TEMP_DEFS, SEE_ION_TABLE_DEFS)
      ELSE
         ALLOCATE(SEE_ION_TABLE_DEFS(1))
      END IF

      N_SEE_ION_TABLE_DEFS = N_SEE_ION_TABLE_DEFS + 1
      SEE_ION_TABLE_DEFS(N_SEE_ION_TABLE_DEFS)%MATERIAL_REF = TRIM(STRARRAY(1))
      SEE_ION_TABLE_DEFS(N_SEE_ION_TABLE_DEFS)%SPECIES_NAME = TRIM(STRARRAY(2))
      SEE_ION_TABLE_DEFS(N_SEE_ION_TABLE_DEFS)%FILENAME = TRIM(STRIP_QUOTES(STRARRAY(3)))

      IF (PROC_ID == 0) THEN
         WRITE(*,'(A,A,A,A,A,A)') '> SEE ion table definition added: material=', &
              TRIM(SEE_ION_TABLE_DEFS(N_SEE_ION_TABLE_DEFS)%MATERIAL_REF), ', species=', &
              TRIM(SEE_ION_TABLE_DEFS(N_SEE_ION_TABLE_DEFS)%SPECIES_NAME), ', file=', &
              TRIM(SEE_ION_TABLE_DEFS(N_SEE_ION_TABLE_DEFS)%FILENAME)
      END IF

      IF (ALLOCATED(STRARRAY)) DEALLOCATE(STRARRAY)

   END SUBROUTINE DEF_SEE_ION_YIELD_TABLE

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE LOAD_SEE_ION_YIELD_TABLES -> Load configured ion tables   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE LOAD_SEE_ION_YIELD_TABLES

      IMPLICIT NONE

      INTEGER :: IDEF

      IF (ALLOCATED(SEE_ION_YIELD_TABLES)) DEALLOCATE(SEE_ION_YIELD_TABLES)
      N_SEE_ION_YIELD_TABLES = 0

      IF (N_SEE_ION_TABLE_DEFS <= 0) THEN
         IF (PROC_ID == 0) THEN
            WRITE(*,*) '> No ion SEE lookup tables configured; using analytic fallback model.'
         END IF
         RETURN
      END IF

      DO IDEF = 1, N_SEE_ION_TABLE_DEFS
         CALL LOAD_SINGLE_SEE_ION_YIELD_TABLE(SEE_ION_TABLE_DEFS(IDEF)%MATERIAL_REF, &
                                              SEE_ION_TABLE_DEFS(IDEF)%SPECIES_NAME, &
                                              SEE_ION_TABLE_DEFS(IDEF)%FILENAME)
      END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '> Loaded ', N_SEE_ION_YIELD_TABLES, ' ion SEE lookup tables'
      END IF

   END SUBROUTINE LOAD_SEE_ION_YIELD_TABLES

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE LOAD_SINGLE_SEE_ION_YIELD_TABLE -> Load one lookup table  !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE LOAD_SINGLE_SEE_ION_YIELD_TABLE(MATERIAL_REF, SPECIES_NAME, FILENAME)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: MATERIAL_REF, SPECIES_NAME, FILENAME
      TYPE(SEE_ION_YIELD_TABLE), DIMENSION(:), ALLOCATABLE :: TEMP_TABLES
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ENERGY_GRID, ANGLE_GRID
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: YIELD_GRID
      INTEGER :: MATERIAL_ID, SPECIES_ID, N_ENERGIES, N_ANGLES

      MATERIAL_ID = RESOLVE_SEE_MATERIAL_ID(MATERIAL_REF)
      IF (MATERIAL_ID < 1) THEN
         CALL ERROR_ABORT('Unknown SEE material reference in ion table definition: ' // TRIM(MATERIAL_REF))
      END IF

      SPECIES_ID = SPECIES_NAME_TO_ID(TRIM(SPECIES_NAME))
      IF (SPECIES_ID < 1) THEN
         CALL ERROR_ABORT('Unknown ion species in SEE ion table definition: ' // TRIM(SPECIES_NAME))
      END IF
      IF (SPECIES(SPECIES_ID)%CHARGE <= 0.d0) THEN
         CALL ERROR_ABORT('SEE ion table species is not positively charged: ' // TRIM(SPECIES_NAME))
      END IF

      CALL READ_SEE_ION_YIELD_TABLE_FILE(FILENAME, ENERGY_GRID, ANGLE_GRID, YIELD_GRID, &
                                         N_ENERGIES, N_ANGLES)

      IF (ALLOCATED(SEE_ION_YIELD_TABLES)) THEN
         ALLOCATE(TEMP_TABLES(N_SEE_ION_YIELD_TABLES + 1))
         TEMP_TABLES(1:N_SEE_ION_YIELD_TABLES) = SEE_ION_YIELD_TABLES(1:N_SEE_ION_YIELD_TABLES)
         CALL MOVE_ALLOC(TEMP_TABLES, SEE_ION_YIELD_TABLES)
      ELSE
         ALLOCATE(SEE_ION_YIELD_TABLES(1))
      END IF

      N_SEE_ION_YIELD_TABLES = N_SEE_ION_YIELD_TABLES + 1
      SEE_ION_YIELD_TABLES(N_SEE_ION_YIELD_TABLES)%MATERIAL_REF = TRIM(MATERIAL_REF)
      SEE_ION_YIELD_TABLES(N_SEE_ION_YIELD_TABLES)%SPECIES_NAME = TRIM(SPECIES_NAME)
      SEE_ION_YIELD_TABLES(N_SEE_ION_YIELD_TABLES)%SOURCE_FILENAME = TRIM(FILENAME)
      SEE_ION_YIELD_TABLES(N_SEE_ION_YIELD_TABLES)%MATERIAL_ID = MATERIAL_ID
      SEE_ION_YIELD_TABLES(N_SEE_ION_YIELD_TABLES)%SPECIES_ID = SPECIES_ID
      SEE_ION_YIELD_TABLES(N_SEE_ION_YIELD_TABLES)%N_ENERGIES = N_ENERGIES
      SEE_ION_YIELD_TABLES(N_SEE_ION_YIELD_TABLES)%N_ANGLES = N_ANGLES
      CALL MOVE_ALLOC(ENERGY_GRID, SEE_ION_YIELD_TABLES(N_SEE_ION_YIELD_TABLES)%ENERGY_GRID)
      CALL MOVE_ALLOC(ANGLE_GRID, SEE_ION_YIELD_TABLES(N_SEE_ION_YIELD_TABLES)%ANGLE_GRID)
      CALL MOVE_ALLOC(YIELD_GRID, SEE_ION_YIELD_TABLES(N_SEE_ION_YIELD_TABLES)%YIELD_GRID)

      IF (PROC_ID == 0) THEN
         WRITE(*,'(A,A,A,A,A,I5,A,I5,A)') '> Loaded ion SEE table: material=', &
              TRIM(SEE_MATERIALS(MATERIAL_ID)%NAME), ', species=', TRIM(SPECIES(SPECIES_ID)%NAME), &
              ', nE=', N_ENERGIES, ', nTheta=', N_ANGLES, ''
      END IF

   END SUBROUTINE LOAD_SINGLE_SEE_ION_YIELD_TABLE

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE READ_SEE_ION_YIELD_TABLE_FILE -> Read gamma(E,theta) data !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE READ_SEE_ION_YIELD_TABLE_FILE(FILENAME, ENERGY_GRID, ANGLE_GRID, YIELD_GRID, &
                                            N_ENERGIES, N_ANGLES)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: ENERGY_GRID, ANGLE_GRID
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: YIELD_GRID
      INTEGER, INTENT(OUT) :: N_ENERGIES, N_ANGLES

      INTEGER, PARAMETER :: IN_UNIT = 49
      CHARACTER(LEN=512) :: LINE
      CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE :: STRARRAY
      INTEGER :: IOS, N_STR, IE, IA
      LOGICAL :: FOUND_ANGLES

      N_ENERGIES = 0
      N_ANGLES = 0
      FOUND_ANGLES = .FALSE.

      OPEN(UNIT=IN_UNIT, FILE=TRIM(FILENAME), STATUS='old', IOSTAT=IOS)
      IF (IOS /= 0) THEN
         CALL ERROR_ABORT('Cannot open SEE ion lookup table: ' // TRIM(FILENAME))
      END IF

      DO
         READ(IN_UNIT, '(A)', IOSTAT=IOS) LINE
         IF (IOS /= 0) EXIT
         LINE = TRIM(STRIP_INLINE_COMMENT(LINE))
         IF (LEN_TRIM(LINE) == 0) CYCLE

         CALL SPLIT_STR(TRIM(LINE), ' ', STRARRAY, N_STR)
         IF (N_STR <= 0) CYCLE

         IF (.NOT. FOUND_ANGLES) THEN
            IF (N_STR < 2) THEN
               CALL ERROR_ABORT('SEE ion lookup table angle header is incomplete: ' // TRIM(FILENAME))
            END IF
            IF (IS_NUMERIC_TOKEN(STRARRAY(1))) THEN
               CALL ERROR_ABORT('First data line in SEE ion lookup table must start with an angle label: ' // TRIM(FILENAME))
            END IF
            N_ANGLES = N_STR - 1
            FOUND_ANGLES = .TRUE.
         ELSE
            IF (N_STR /= N_ANGLES + 1) THEN
               CALL ERROR_ABORT('SEE ion lookup table row width mismatch: ' // TRIM(FILENAME))
            END IF
            N_ENERGIES = N_ENERGIES + 1
         END IF

         IF (ALLOCATED(STRARRAY)) DEALLOCATE(STRARRAY)
      END DO
      CLOSE(IN_UNIT)

      IF (.NOT. FOUND_ANGLES .OR. N_ANGLES <= 0 .OR. N_ENERGIES <= 0) THEN
         CALL ERROR_ABORT('SEE ion lookup table is empty or incomplete: ' // TRIM(FILENAME))
      END IF

      ALLOCATE(ANGLE_GRID(N_ANGLES), ENERGY_GRID(N_ENERGIES), YIELD_GRID(N_ENERGIES, N_ANGLES))

      OPEN(UNIT=IN_UNIT, FILE=TRIM(FILENAME), STATUS='old', IOSTAT=IOS)
      IF (IOS /= 0) THEN
         CALL ERROR_ABORT('Cannot reopen SEE ion lookup table: ' // TRIM(FILENAME))
      END IF

      FOUND_ANGLES = .FALSE.
      IE = 0
      DO
         READ(IN_UNIT, '(A)', IOSTAT=IOS) LINE
         IF (IOS /= 0) EXIT
         LINE = TRIM(STRIP_INLINE_COMMENT(LINE))
         IF (LEN_TRIM(LINE) == 0) CYCLE

         CALL SPLIT_STR(TRIM(LINE), ' ', STRARRAY, N_STR)
         IF (N_STR <= 0) CYCLE

         IF (.NOT. FOUND_ANGLES) THEN
            DO IA = 1, N_ANGLES
               CALL READ_TABLE_REAL(STRARRAY(IA + 1), ANGLE_GRID(IA), 'angle header', FILENAME)
               IF (IA > 1 .AND. ANGLE_GRID(IA) <= ANGLE_GRID(IA - 1)) THEN
                  CALL ERROR_ABORT('SEE ion lookup table angles must be strictly increasing: ' // TRIM(FILENAME))
               END IF
            END DO
            FOUND_ANGLES = .TRUE.
         ELSE
            IE = IE + 1
            CALL READ_TABLE_REAL(STRARRAY(1), ENERGY_GRID(IE), 'energy row', FILENAME)
            IF (IE > 1 .AND. ENERGY_GRID(IE) <= ENERGY_GRID(IE - 1)) THEN
               CALL ERROR_ABORT('SEE ion lookup table energies must be strictly increasing: ' // TRIM(FILENAME))
            END IF
            DO IA = 1, N_ANGLES
               CALL READ_TABLE_REAL(STRARRAY(IA + 1), YIELD_GRID(IE, IA), 'yield value', FILENAME)
            END DO
         END IF

         IF (ALLOCATED(STRARRAY)) DEALLOCATE(STRARRAY)
      END DO
      CLOSE(IN_UNIT)

   END SUBROUTINE READ_SEE_ION_YIELD_TABLE_FILE

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE READ_TABLE_REAL -> Parse one lookup-table number          !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE READ_TABLE_REAL(TOKEN, VALUE, CONTEXT, FILENAME)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: TOKEN, CONTEXT, FILENAME
      REAL(KIND=8), INTENT(OUT) :: VALUE
      INTEGER :: IOS

      READ(TOKEN, *, IOSTAT=IOS) VALUE
      IF (IOS /= 0) THEN
         CALL ERROR_ABORT('Cannot parse ' // TRIM(CONTEXT) // ' in SEE ion lookup table: ' // TRIM(FILENAME))
      END IF

   END SUBROUTINE READ_TABLE_REAL

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION RESOLVE_SEE_MATERIAL_ID -> Resolve material name or ID      !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION RESOLVE_SEE_MATERIAL_ID(MATERIAL_REF) RESULT(MATERIAL_ID)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: MATERIAL_REF
      INTEGER :: MATERIAL_ID
      INTEGER :: IOS, IM

      READ(MATERIAL_REF, *, IOSTAT=IOS) MATERIAL_ID
      IF (IOS == 0) THEN
         IF (MATERIAL_ID < 1 .OR. MATERIAL_ID > N_SEE_MATERIALS) MATERIAL_ID = -1
         RETURN
      END IF

      MATERIAL_ID = -1
      IF (.NOT. ALLOCATED(SEE_MATERIALS)) RETURN
      DO IM = 1, N_SEE_MATERIALS
         IF (TRIM(SEE_MATERIALS(IM)%NAME) == TRIM(MATERIAL_REF)) THEN
            MATERIAL_ID = IM
            RETURN
         END IF
      END DO

   END FUNCTION RESOLVE_SEE_MATERIAL_ID

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION FIND_SEE_ION_YIELD_TABLE_INDEX -> Match material+species    !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION FIND_SEE_ION_YIELD_TABLE_INDEX(MATERIAL_ID, SPECIES_ID) RESULT(TABLE_INDEX)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: MATERIAL_ID, SPECIES_ID
      INTEGER :: TABLE_INDEX
      INTEGER :: IT

      TABLE_INDEX = -1
      IF (.NOT. ALLOCATED(SEE_ION_YIELD_TABLES) .OR. N_SEE_ION_YIELD_TABLES <= 0) RETURN

      DO IT = N_SEE_ION_YIELD_TABLES, 1, -1
         IF (SEE_ION_YIELD_TABLES(IT)%MATERIAL_ID == MATERIAL_ID .AND. &
             SEE_ION_YIELD_TABLES(IT)%SPECIES_ID == SPECIES_ID) THEN
            TABLE_INDEX = IT
            RETURN
         END IF
      END DO

   END FUNCTION FIND_SEE_ION_YIELD_TABLE_INDEX

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION LOOKUP_SEE_ION_YIELD -> Bilinear interpolation of gamma     !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION LOOKUP_SEE_ION_YIELD(TABLE_INDEX, ION_ENERGY_EV, ION_ANGLE_DEGREES) RESULT(GAMMA)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: TABLE_INDEX
      REAL(KIND=8), INTENT(IN) :: ION_ENERGY_EV, ION_ANGLE_DEGREES
      REAL(KIND=8) :: GAMMA
      REAL(KIND=8) :: CLAMPED_ANGLE, W_ENERGY, W_ANGLE
      REAL(KIND=8) :: G00, G01, G10, G11, G0, G1
      INTEGER :: IE0, IE1, IA0, IA1

      GAMMA = 0.d0
      IF (TABLE_INDEX < 1 .OR. TABLE_INDEX > N_SEE_ION_YIELD_TABLES) RETURN
      IF (ION_ENERGY_EV <= 0.d0) RETURN
      IF (SEE_ION_YIELD_TABLES(TABLE_INDEX)%N_ENERGIES <= 0) RETURN
      IF (SEE_ION_YIELD_TABLES(TABLE_INDEX)%N_ANGLES <= 0) RETURN

      IF (ION_ENERGY_EV < SEE_ION_YIELD_TABLES(TABLE_INDEX)%ENERGY_GRID(1)) RETURN

      CLAMPED_ANGLE = MAX(SEE_ION_YIELD_TABLES(TABLE_INDEX)%ANGLE_GRID(1), &
                          MIN(ION_ANGLE_DEGREES, &
                              SEE_ION_YIELD_TABLES(TABLE_INDEX)%ANGLE_GRID(SEE_ION_YIELD_TABLES(TABLE_INDEX)%N_ANGLES)))

      CALL FIND_BRACKETING_INTERVAL(SEE_ION_YIELD_TABLES(TABLE_INDEX)%ENERGY_GRID, &
                                    SEE_ION_YIELD_TABLES(TABLE_INDEX)%N_ENERGIES, ION_ENERGY_EV, &
                                    IE0, IE1, W_ENERGY)
      CALL FIND_BRACKETING_INTERVAL(SEE_ION_YIELD_TABLES(TABLE_INDEX)%ANGLE_GRID, &
                                    SEE_ION_YIELD_TABLES(TABLE_INDEX)%N_ANGLES, CLAMPED_ANGLE, &
                                    IA0, IA1, W_ANGLE)

      G00 = SEE_ION_YIELD_TABLES(TABLE_INDEX)%YIELD_GRID(IE0, IA0)
      G01 = SEE_ION_YIELD_TABLES(TABLE_INDEX)%YIELD_GRID(IE0, IA1)
      G10 = SEE_ION_YIELD_TABLES(TABLE_INDEX)%YIELD_GRID(IE1, IA0)
      G11 = SEE_ION_YIELD_TABLES(TABLE_INDEX)%YIELD_GRID(IE1, IA1)

      IF (IA0 == IA1) THEN
         G0 = G00
         G1 = G10
      ELSE
         G0 = (1.d0 - W_ANGLE) * G00 + W_ANGLE * G01
         G1 = (1.d0 - W_ANGLE) * G10 + W_ANGLE * G11
      END IF

      IF (IE0 == IE1) THEN
         GAMMA = G0
      ELSE
         GAMMA = (1.d0 - W_ENERGY) * G0 + W_ENERGY * G1
      END IF

      GAMMA = MAX(0.d0, GAMMA)

   END FUNCTION LOOKUP_SEE_ION_YIELD

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION LOOKUP_SEE_ION_YIELD_FIRST_ROW -> First-row angle lookup    !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION LOOKUP_SEE_ION_YIELD_FIRST_ROW(TABLE_INDEX, ION_ANGLE_DEGREES) RESULT(GAMMA)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: TABLE_INDEX
      REAL(KIND=8), INTENT(IN) :: ION_ANGLE_DEGREES
      REAL(KIND=8) :: GAMMA
      REAL(KIND=8) :: CLAMPED_ANGLE, W_ANGLE
      REAL(KIND=8) :: G0, G1
      INTEGER :: IA0, IA1

      GAMMA = 0.d0
      IF (TABLE_INDEX < 1 .OR. TABLE_INDEX > N_SEE_ION_YIELD_TABLES) RETURN
      IF (SEE_ION_YIELD_TABLES(TABLE_INDEX)%N_ANGLES <= 0) RETURN

      CLAMPED_ANGLE = MAX(SEE_ION_YIELD_TABLES(TABLE_INDEX)%ANGLE_GRID(1), &
                          MIN(ION_ANGLE_DEGREES, &
                              SEE_ION_YIELD_TABLES(TABLE_INDEX)%ANGLE_GRID(SEE_ION_YIELD_TABLES(TABLE_INDEX)%N_ANGLES)))

      CALL FIND_BRACKETING_INTERVAL(SEE_ION_YIELD_TABLES(TABLE_INDEX)%ANGLE_GRID, &
                                    SEE_ION_YIELD_TABLES(TABLE_INDEX)%N_ANGLES, CLAMPED_ANGLE, &
                                    IA0, IA1, W_ANGLE)

      G0 = SEE_ION_YIELD_TABLES(TABLE_INDEX)%YIELD_GRID(1, IA0)
      G1 = SEE_ION_YIELD_TABLES(TABLE_INDEX)%YIELD_GRID(1, IA1)

      IF (IA0 == IA1) THEN
         GAMMA = G0
      ELSE
         GAMMA = (1.d0 - W_ANGLE) * G0 + W_ANGLE * G1
      END IF

      GAMMA = MAX(0.d0, GAMMA)

   END FUNCTION LOOKUP_SEE_ION_YIELD_FIRST_ROW

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION CALCULATE_SUBTHRESHOLD_ION_YIELD -> low-energy continuation !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION CALCULATE_SUBTHRESHOLD_ION_YIELD(ION_SPECIES_ID, ION_ENERGY_EV, &
                                             ION_ANGLE_DEGREES, ANCHOR_ENERGY_EV, &
                                             ANCHOR_GAMMA) RESULT(GAMMA)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ION_SPECIES_ID
      REAL(KIND=8), INTENT(IN) :: ION_ENERGY_EV, ION_ANGLE_DEGREES
      REAL(KIND=8), INTENT(IN) :: ANCHOR_ENERGY_EV, ANCHOR_GAMMA
      REAL(KIND=8) :: GAMMA
      REAL(KIND=8) :: VELOCITY, ANCHOR_VELOCITY, COS_ANGLE
      REAL(KIND=8) :: V_PERP, V_ANCHOR_PERP, A_PREFACTOR, SLOPE_B

      GAMMA = 0.d0
      IF (ION_SPECIES_ID < 1 .OR. ION_SPECIES_ID > N_SPECIES) RETURN
      IF (ION_ENERGY_EV <= 0.d0 .OR. ANCHOR_ENERGY_EV <= 0.d0) RETURN
      IF (ANCHOR_GAMMA <= 0.d0) RETURN
      IF (SPECIES(ION_SPECIES_ID)%MOLECULAR_MASS <= 0.d0) RETURN

      COS_ANGLE = MAX(COS(ION_ANGLE_DEGREES * PI / 180.d0), 0.d0)
      IF (COS_ANGLE <= 0.d0) RETURN

      VELOCITY = SQRT(2.d0 * ION_ENERGY_EV * QE / SPECIES(ION_SPECIES_ID)%MOLECULAR_MASS)
      ANCHOR_VELOCITY = SQRT(2.d0 * ANCHOR_ENERGY_EV * QE / SPECIES(ION_SPECIES_ID)%MOLECULAR_MASS)

      V_PERP = VELOCITY * COS_ANGLE
      V_ANCHOR_PERP = ANCHOR_VELOCITY * COS_ANGLE
      IF (V_PERP <= 0.d0 .OR. V_ANCHOR_PERP <= 0.d0) RETURN

      ! Literature-inspired subthreshold continuation: Gamma = A exp(-B / v_perp),
      ! consistent with the 1/v dependence reported by Sroubek et al.
      ! (NIMB 268, 3377-3380, 2010; Surface Science 625, 7-9, 2014).
      ! A and B are inferred from the anchor point rather than assumed universal.
      A_PREFACTOR = 1.d0
      IF (ANCHOR_GAMMA >= 0.999d0) THEN
         A_PREFACTOR = MAX(1.05d0 * ANCHOR_GAMMA, ANCHOR_GAMMA + 1.d-6)
      END IF

      SLOPE_B = V_ANCHOR_PERP * LOG(A_PREFACTOR / ANCHOR_GAMMA)
      IF (SLOPE_B <= 0.d0) THEN
         GAMMA = MAX(0.d0, ANCHOR_GAMMA * V_PERP / V_ANCHOR_PERP)
      ELSE
         GAMMA = A_PREFACTOR * EXP(-SLOPE_B / V_PERP)
      END IF

      GAMMA = MIN(MAX(0.d0, GAMMA), MAX(A_PREFACTOR, ANCHOR_GAMMA))

   END FUNCTION CALCULATE_SUBTHRESHOLD_ION_YIELD

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION GET_PROJECTILE_EFFECTIVE_IONIZATION_ENERGY -> potential energy estimate
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION GET_PROJECTILE_EFFECTIVE_IONIZATION_ENERGY(ION_SPECIES_ID, ION_CHARGE) RESULT(IONIZATION_ENERGY_EV)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ION_SPECIES_ID
      REAL(KIND=8), INTENT(IN) :: ION_CHARGE
      REAL(KIND=8) :: IONIZATION_ENERGY_EV
      CHARACTER(LEN=64) :: SPECIES_NAME
      REAL(KIND=8) :: CHARGE_MAG

      IONIZATION_ENERGY_EV = 0.d0
      CHARGE_MAG = MAX(ABS(ION_CHARGE), 1.d0)

      IF (ION_SPECIES_ID < 1 .OR. ION_SPECIES_ID > N_SPECIES) THEN
         IONIZATION_ENERGY_EV = 13.6d0 * CHARGE_MAG
         RETURN
      END IF

      SPECIES_NAME = ADJUSTL(TRIM(SPECIES(ION_SPECIES_ID)%NAME))

      SELECT CASE (TRIM(SPECIES_NAME))
      CASE ('H+', 'D+')
         IONIZATION_ENERGY_EV = 13.598d0
      CASE ('H2+', 'D2+')
         IONIZATION_ENERGY_EV = 15.426d0
      CASE ('H3+', 'D3+')
         IONIZATION_ENERGY_EV = 13.6d0
      CASE ('He+')
         IONIZATION_ENERGY_EV = 24.587d0
      CASE ('Ne+')
         IONIZATION_ENERGY_EV = 21.565d0
      CASE ('Ar+')
         IONIZATION_ENERGY_EV = 15.760d0
      CASE ('Kr+')
         IONIZATION_ENERGY_EV = 13.999d0
      CASE ('Xe+')
         IONIZATION_ENERGY_EV = 12.130d0
      CASE ('N+', 'N2+')
         IONIZATION_ENERGY_EV = 15.580d0
      CASE ('O+')
         IONIZATION_ENERGY_EV = 13.618d0
      CASE ('O2+')
         IONIZATION_ENERGY_EV = 12.070d0
      CASE DEFAULT
         ! Fallback for unsupported projectiles: use a first-ionization-scale estimate.
         IONIZATION_ENERGY_EV = 13.6d0 * CHARGE_MAG
      END SELECT

      IONIZATION_ENERGY_EV = IONIZATION_ENERGY_EV * CHARGE_MAG

   END FUNCTION GET_PROJECTILE_EFFECTIVE_IONIZATION_ENERGY

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION CALCULATE_LOGNORMAL_BRANCH -> peaked shape in impact energy !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION CALCULATE_LOGNORMAL_BRANCH(IMPACT_ENERGY_EV, PEAK_ENERGY_EV, SIGMA_LOG) RESULT(BRANCH_VALUE)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: IMPACT_ENERGY_EV, PEAK_ENERGY_EV, SIGMA_LOG
      REAL(KIND=8) :: BRANCH_VALUE
      REAL(KIND=8) :: CLAMPED_SIGMA

      BRANCH_VALUE = 0.d0
      IF (IMPACT_ENERGY_EV <= 0.d0 .OR. PEAK_ENERGY_EV <= 0.d0) RETURN

      CLAMPED_SIGMA = MAX(SIGMA_LOG, 0.2d0)
      BRANCH_VALUE = EXP(-0.5d0 * (LOG(IMPACT_ENERGY_EV / PEAK_ENERGY_EV) / CLAMPED_SIGMA)**2)

   END FUNCTION CALCULATE_LOGNORMAL_BRANCH

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION CALCULATE_SOFT_ONSET_FACTOR -> smooth low-energy suppression !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION CALCULATE_SOFT_ONSET_FACTOR(IMPACT_ENERGY_EV, ONSET_ENERGY_EV) RESULT(ONSET_FACTOR)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: IMPACT_ENERGY_EV, ONSET_ENERGY_EV
      REAL(KIND=8) :: ONSET_FACTOR

      IF (IMPACT_ENERGY_EV <= 0.d0) THEN
         ONSET_FACTOR = 0.d0
      ELSE IF (ONSET_ENERGY_EV <= 0.d0) THEN
         ONSET_FACTOR = 1.d0
      ELSE
         ONSET_FACTOR = 1.d0 - EXP(-(IMPACT_ENERGY_EV / ONSET_ENERGY_EV)**2)
      END IF

   END FUNCTION CALCULATE_SOFT_ONSET_FACTOR

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION CALCULATE_SEE_ANGLE_FACTOR -> capped angular enhancement    !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION CALCULATE_SEE_ANGLE_FACTOR(IMPACT_ANGLE_DEGREES, ANGLE_EXPONENT, MAX_FACTOR) RESULT(ANGLE_FACTOR)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: IMPACT_ANGLE_DEGREES, ANGLE_EXPONENT, MAX_FACTOR
      REAL(KIND=8) :: ANGLE_FACTOR
      REAL(KIND=8) :: COS_ANGLE, CLAMPED_COS

      ANGLE_FACTOR = 0.d0

      COS_ANGLE = COS(IMPACT_ANGLE_DEGREES * PI / 180.d0)
      CLAMPED_COS = MAX(COS_ANGLE, 5.d-2)
      IF (COS_ANGLE <= 0.d0) RETURN

      ! Stronger-than-linear angular enhancement with a hard cap to avoid
      ! unphysical divergence at grazing incidence.
      ANGLE_FACTOR = CLAMPED_COS**(-MAX(ANGLE_EXPONENT, 0.d0))
      ANGLE_FACTOR = MIN(MAX(ANGLE_FACTOR, 1.d0), MAX(MAX_FACTOR, 1.d0))

   END FUNCTION CALCULATE_SEE_ANGLE_FACTOR

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION CALCULATE_ENERGY_MODULATED_ANGLE_FACTOR -> weak at low E    !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION CALCULATE_ENERGY_MODULATED_ANGLE_FACTOR(IMPACT_ANGLE_DEGREES, ANGLE_EXPONENT, &
                                                    MAX_FACTOR, IMPACT_ENERGY_EV, &
                                                    REFERENCE_ENERGY_EV, RISE_POWER) RESULT(ANGLE_FACTOR)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: IMPACT_ANGLE_DEGREES, ANGLE_EXPONENT, MAX_FACTOR
      REAL(KIND=8), INTENT(IN) :: IMPACT_ENERGY_EV, REFERENCE_ENERGY_EV, RISE_POWER
      REAL(KIND=8) :: ANGLE_FACTOR
      REAL(KIND=8) :: ENERGY_WEIGHT, ENERGY_RATIO, EFFECTIVE_EXPONENT, EFFECTIVE_MAX

      ANGLE_FACTOR = 1.d0
      IF (IMPACT_ENERGY_EV <= 0.d0) RETURN

      IF (REFERENCE_ENERGY_EV <= 0.d0) THEN
         ANGLE_FACTOR = CALCULATE_SEE_ANGLE_FACTOR(IMPACT_ANGLE_DEGREES, ANGLE_EXPONENT, MAX_FACTOR)
         RETURN
      END IF

      ENERGY_RATIO = MAX(IMPACT_ENERGY_EV / REFERENCE_ENERGY_EV, 0.d0)

      ! Conservative low-energy behavior: angular enhancement grows only
      ! gradually with impact energy and approaches the material-defined
      ! high-energy limit asymptotically.
      ENERGY_WEIGHT = 1.d0 - EXP(-(ENERGY_RATIO**MAX(RISE_POWER, 1.d0)))
      ENERGY_WEIGHT = MIN(MAX(ENERGY_WEIGHT, 0.d0), 1.d0)

      EFFECTIVE_EXPONENT = ANGLE_EXPONENT * ENERGY_WEIGHT
      EFFECTIVE_MAX = 1.d0 + (MAX_FACTOR - 1.d0) * ENERGY_WEIGHT

      ANGLE_FACTOR = CALCULATE_SEE_ANGLE_FACTOR(IMPACT_ANGLE_DEGREES, EFFECTIVE_EXPONENT, EFFECTIVE_MAX)

   END FUNCTION CALCULATE_ENERGY_MODULATED_ANGLE_FACTOR

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION CALCULATE_KINETIC_SEE_BRANCH -> high-energy kinetic peak     !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION CALCULATE_KINETIC_SEE_BRANCH(PEAK_YIELD, PEAK_ENERGY_EV, SIGMA_LOG, &
                                         ONSET_ENERGY_EV, IMPACT_ENERGY_EV, &
                                         IMPACT_ANGLE_DEGREES, ANGLE_EXPONENT, &
                                         ANGLE_MAX_FACTOR) RESULT(GAMMA)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: PEAK_YIELD, PEAK_ENERGY_EV, SIGMA_LOG, ONSET_ENERGY_EV
      REAL(KIND=8), INTENT(IN) :: IMPACT_ENERGY_EV, IMPACT_ANGLE_DEGREES
      REAL(KIND=8), INTENT(IN) :: ANGLE_EXPONENT, ANGLE_MAX_FACTOR
      REAL(KIND=8) :: GAMMA
      REAL(KIND=8) :: COS_ANGLE, ANGLE_FACTOR, SHAPE_FACTOR, ONSET_FACTOR
      REAL(KIND=8) :: ANGLE_REFERENCE_ENERGY

      GAMMA = 0.d0
      IF (PEAK_YIELD <= 0.d0 .OR. IMPACT_ENERGY_EV <= 0.d0 .OR. PEAK_ENERGY_EV <= 0.d0) RETURN

      COS_ANGLE = MAX(COS(IMPACT_ANGLE_DEGREES * PI / 180.d0), 0.d0)
      IF (COS_ANGLE <= 0.d0) RETURN

      ! Let angular enhancement ramp up when the kinetic branch starts to
      ! become relevant, not only when the full high-energy peak is reached.
      ! This keeps low-energy ions conservative, while making keV-scale
      ! oblique impacts measurably larger than normal incidence.
      ANGLE_REFERENCE_ENERGY = MAX(ONSET_ENERGY_EV, 1.d0)
      ANGLE_REFERENCE_ENERGY = SQRT(MAX(PEAK_ENERGY_EV, ANGLE_REFERENCE_ENERGY) * ANGLE_REFERENCE_ENERGY)

      ANGLE_FACTOR = CALCULATE_ENERGY_MODULATED_ANGLE_FACTOR(IMPACT_ANGLE_DEGREES, ANGLE_EXPONENT, &
                                                             ANGLE_MAX_FACTOR, IMPACT_ENERGY_EV, &
                                                             ANGLE_REFERENCE_ENERGY, 1.5d0)
      SHAPE_FACTOR = CALCULATE_LOGNORMAL_BRANCH(IMPACT_ENERGY_EV, PEAK_ENERGY_EV, SIGMA_LOG)
      ONSET_FACTOR = CALCULATE_SOFT_ONSET_FACTOR(IMPACT_ENERGY_EV, ONSET_ENERGY_EV)

      GAMMA = PEAK_YIELD * SHAPE_FACTOR * ONSET_FACTOR * ANGLE_FACTOR

   END FUNCTION CALCULATE_KINETIC_SEE_BRANCH

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION CALCULATE_POTENTIAL_ION_BRANCH -> low-energy Hagstrum-like branch
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION CALCULATE_POTENTIAL_ION_BRANCH(MATERIAL_ID, ION_SPECIES_ID, ION_ENERGY_EV, &
                                           ION_CHARGE, ION_ANGLE_DEGREES) RESULT(GAMMA)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: MATERIAL_ID, ION_SPECIES_ID
      REAL(KIND=8), INTENT(IN) :: ION_ENERGY_EV, ION_CHARGE, ION_ANGLE_DEGREES
      REAL(KIND=8) :: GAMMA
      REAL(KIND=8) :: WORK_FUNCTION, IONIZATION_ENERGY_EV, AVAILABLE_POTENTIAL_ENERGY
      REAL(KIND=8) :: PEAK_YIELD, PEAK_ENERGY_EV, SIGMA_LOG
      REAL(KIND=8) :: COS_ANGLE, ANGLE_FACTOR, CHARGE_SCALING

      GAMMA = 0.d0
      IF (MATERIAL_ID < 1 .OR. MATERIAL_ID > N_SEE_MATERIALS) RETURN
      IF (ION_ENERGY_EV <= 0.d0 .OR. ION_CHARGE <= 0.d0) RETURN

      WORK_FUNCTION = MAX(SEE_MATERIALS(MATERIAL_ID)%W_WORK_FUNCTION, 0.d0)
      IONIZATION_ENERGY_EV = GET_PROJECTILE_EFFECTIVE_IONIZATION_ENERGY(ION_SPECIES_ID, ION_CHARGE)
      AVAILABLE_POTENTIAL_ENERGY = MAX(0.d0, 0.78d0 * IONIZATION_ENERGY_EV - 2.d0 * WORK_FUNCTION)
      IF (AVAILABLE_POTENTIAL_ENERGY <= 0.d0) RETURN

      ! Empirical peak-yield estimate for singly charged ions from Baragiola et al.
      ! as quoted in FIB/SEE reviews: gamma_p ≈ 0.032 * (0.78 * I - 2 * Phi),
      ! where I is projectile ionization energy and Phi is the work function.
      CHARGE_SCALING = MAX(SEE_MATERIALS(MATERIAL_ID)%ION_CHARGE_FACTOR, 0.d0)
      CHARGE_SCALING = MAX(1.d0, CHARGE_SCALING * MAX(ION_CHARGE, 1.d0))
      PEAK_YIELD = 0.032d0 * AVAILABLE_POTENTIAL_ENERGY * CHARGE_SCALING
      PEAK_YIELD = MIN(MAX(PEAK_YIELD, 0.d0), 2.d0)

      PEAK_ENERGY_EV = MAX(5.d0, MIN(MAX(SEE_MATERIALS(MATERIAL_ID)%ION_BETA, 5.d0), &
                                     0.7d0 * MAX(SEE_MATERIALS(MATERIAL_ID)%ION_E_THRESHOLD, 5.d0)))
      SIGMA_LOG = MAX(SEE_MATERIALS(MATERIAL_ID)%ION_ALPHA, 0.2d0)

      COS_ANGLE = MAX(COS(ION_ANGLE_DEGREES * PI / 180.d0), 0.d0)
      IF (COS_ANGLE <= 0.d0) RETURN
      ANGLE_FACTOR = CALCULATE_ENERGY_MODULATED_ANGLE_FACTOR(ION_ANGLE_DEGREES, &
                     SEE_MATERIALS(MATERIAL_ID)%ION_POT_ANGLE_EXP, &
                     SEE_MATERIALS(MATERIAL_ID)%ION_POT_ANGLE_MAX, ION_ENERGY_EV, &
                     PEAK_ENERGY_EV, 2.0d0)

      GAMMA = PEAK_YIELD * CALCULATE_LOGNORMAL_BRANCH(ION_ENERGY_EV, PEAK_ENERGY_EV, SIGMA_LOG) * ANGLE_FACTOR

   END FUNCTION CALCULATE_POTENTIAL_ION_BRANCH

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION CALCULATE_EMPIRICAL_ION_SEE_YIELD -> piecewise potential + kinetic
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION CALCULATE_EMPIRICAL_ION_SEE_YIELD(MATERIAL_ID, ION_SPECIES_ID, ION_ENERGY_EV, &
                                              ION_CHARGE, ION_ANGLE_DEGREES) RESULT(GAMMA)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: MATERIAL_ID, ION_SPECIES_ID
      REAL(KIND=8), INTENT(IN) :: ION_ENERGY_EV, ION_CHARGE, ION_ANGLE_DEGREES
      REAL(KIND=8) :: GAMMA
      REAL(KIND=8) :: POTENTIAL_GAMMA, KINETIC_GAMMA
      REAL(KIND=8) :: POTENTIAL_PEAK_ENERGY, KINETIC_PEAK_ENERGY, TRANSITION_ENERGY
      REAL(KIND=8) :: TRANSITION_LOW, TRANSITION_HIGH, BLEND_WEIGHT
      REAL(KIND=8) :: KINETIC_ONSET, KINETIC_SIGMA, KINETIC_PEAK_YIELD

      GAMMA = 0.d0
      IF (MATERIAL_ID < 1 .OR. MATERIAL_ID > N_SEE_MATERIALS) RETURN
      IF (ION_ENERGY_EV <= 0.d0 .OR. ION_CHARGE <= 0.d0) RETURN

      POTENTIAL_GAMMA = CALCULATE_POTENTIAL_ION_BRANCH(MATERIAL_ID, ION_SPECIES_ID, ION_ENERGY_EV, &
                                                       ION_CHARGE, ION_ANGLE_DEGREES)

      POTENTIAL_PEAK_ENERGY = MAX(5.d0, MIN(MAX(SEE_MATERIALS(MATERIAL_ID)%ION_BETA, 5.d0), &
                                            0.7d0 * MAX(SEE_MATERIALS(MATERIAL_ID)%ION_E_THRESHOLD, 5.d0)))
      KINETIC_PEAK_ENERGY = MAX(SEE_MATERIALS(MATERIAL_ID)%ION_E_THRESHOLD, 1.5d0 * POTENTIAL_PEAK_ENERGY)
      KINETIC_SIGMA = MAX(SEE_MATERIALS(MATERIAL_ID)%ION_ALPHA, 0.2d0)
      KINETIC_ONSET = MAX(5.d0, MIN(0.6d0 * KINETIC_PEAK_ENERGY, MAX(SEE_MATERIALS(MATERIAL_ID)%ION_BETA, 5.d0)))
      KINETIC_PEAK_YIELD = MAX(SEE_MATERIALS(MATERIAL_ID)%GAMMA_MAX, 0.d0)
      KINETIC_GAMMA = CALCULATE_KINETIC_SEE_BRANCH(KINETIC_PEAK_YIELD, KINETIC_PEAK_ENERGY, &
                                                   KINETIC_SIGMA, KINETIC_ONSET, ION_ENERGY_EV, &
                                                   ION_ANGLE_DEGREES, &
                                                   SEE_MATERIALS(MATERIAL_ID)%ION_KIN_ANGLE_EXP, &
                                                   SEE_MATERIALS(MATERIAL_ID)%ION_KIN_ANGLE_MAX)

      TRANSITION_ENERGY = SQRT(POTENTIAL_PEAK_ENERGY * KINETIC_PEAK_ENERGY)
      TRANSITION_LOW = MAX(POTENTIAL_PEAK_ENERGY, TRANSITION_ENERGY / 1.5d0)
      TRANSITION_HIGH = MAX(TRANSITION_LOW + 1.d-9, TRANSITION_ENERGY * 1.5d0)

      IF (ION_ENERGY_EV <= TRANSITION_LOW) THEN
         GAMMA = POTENTIAL_GAMMA
      ELSE IF (ION_ENERGY_EV >= TRANSITION_HIGH) THEN
         GAMMA = KINETIC_GAMMA
      ELSE
         BLEND_WEIGHT = LOG(ION_ENERGY_EV / TRANSITION_LOW) / LOG(TRANSITION_HIGH / TRANSITION_LOW)
         BLEND_WEIGHT = MIN(MAX(BLEND_WEIGHT, 0.d0), 1.d0)
         GAMMA = (1.d0 - BLEND_WEIGHT) * POTENTIAL_GAMMA + BLEND_WEIGHT * KINETIC_GAMMA
      END IF

      GAMMA = MAX(0.d0, MIN(GAMMA, 10.d0))

   END FUNCTION CALCULATE_EMPIRICAL_ION_SEE_YIELD
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE FIND_BRACKETING_INTERVAL -> Locate interpolation segment  !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE FIND_BRACKETING_INTERVAL(GRID, N_GRID, VALUE, LOWER_INDEX, UPPER_INDEX, WEIGHT)

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: GRID
      INTEGER, INTENT(IN) :: N_GRID
      REAL(KIND=8), INTENT(IN) :: VALUE
      INTEGER, INTENT(OUT) :: LOWER_INDEX, UPPER_INDEX
      REAL(KIND=8), INTENT(OUT) :: WEIGHT
      INTEGER :: I
      REAL(KIND=8) :: SPAN

      LOWER_INDEX = 1
      UPPER_INDEX = 1
      WEIGHT = 0.d0

      IF (N_GRID <= 1 .OR. VALUE <= GRID(1)) RETURN
      IF (VALUE >= GRID(N_GRID)) THEN
         LOWER_INDEX = N_GRID
         UPPER_INDEX = N_GRID
         RETURN
      END IF

      DO I = 1, N_GRID - 1
         IF (VALUE <= GRID(I + 1)) THEN
            LOWER_INDEX = I
            UPPER_INDEX = I + 1
            SPAN = GRID(I + 1) - GRID(I)
            IF (SPAN > 0.d0) WEIGHT = (VALUE - GRID(I)) / SPAN
            RETURN
         END IF
      END DO

      LOWER_INDEX = N_GRID
      UPPER_INDEX = N_GRID

   END SUBROUTINE FIND_BRACKETING_INTERVAL

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION STRIP_INLINE_COMMENT -> Remove trailing ! or # comments      !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION STRIP_INLINE_COMMENT(TEXT) RESULT(CLEAN)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: TEXT
      CHARACTER(LEN=LEN(TEXT)) :: CLEAN
      INTEGER :: POS_EXCL, POS_HASH, POS_CUT

      CLEAN = TEXT
      POS_EXCL = INDEX(CLEAN, '!')
      POS_HASH = INDEX(CLEAN, '#')
      POS_CUT = 0

      IF (POS_EXCL > 0) POS_CUT = POS_EXCL
      IF (POS_HASH > 0) THEN
         IF (POS_CUT == 0 .OR. POS_HASH < POS_CUT) POS_CUT = POS_HASH
      END IF
      IF (POS_CUT > 0) CLEAN = CLEAN(:POS_CUT - 1)

   END FUNCTION STRIP_INLINE_COMMENT

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION STRIP_QUOTES -> Remove matching leading/trailing quotes      !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION STRIP_QUOTES(TEXT) RESULT(UNQUOTED)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: TEXT
      CHARACTER(LEN=LEN(TEXT)) :: UNQUOTED
      INTEGER :: LAST_CHAR

      UNQUOTED = ADJUSTL(TEXT)
      LAST_CHAR = LEN_TRIM(UNQUOTED)
      IF (LAST_CHAR >= 2) THEN
         IF ((UNQUOTED(1:1) == ''' .AND. UNQUOTED(LAST_CHAR:LAST_CHAR) == ''') .OR. &
             (UNQUOTED(1:1) == '"' .AND. UNQUOTED(LAST_CHAR:LAST_CHAR) == '"')) THEN
            UNQUOTED = UNQUOTED(2:LAST_CHAR - 1)
         END IF
      END IF

   END FUNCTION STRIP_QUOTES

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION IS_NUMERIC_TOKEN -> Check whether a token is numeric         !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION IS_NUMERIC_TOKEN(TOKEN) RESULT(IS_NUMERIC)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: TOKEN
      LOGICAL :: IS_NUMERIC
      REAL(KIND=8) :: VALUE
      INTEGER :: IOS

      READ(TOKEN, *, IOSTAT=IOS) VALUE
      IS_NUMERIC = IOS == 0

   END FUNCTION IS_NUMERIC_TOKEN

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION CALCULATE_ELECTRON_VAUGHAN_CORE -> Base electron SEE branch   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION CALCULATE_ELECTRON_VAUGHAN_CORE(MATERIAL_ID, IMPACT_ENERGY_EV) RESULT(DELTA_CORE)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: MATERIAL_ID
      REAL(KIND=8), INTENT(IN) :: IMPACT_ENERGY_EV
      REAL(KIND=8) :: DELTA_CORE
      
      REAL(KIND=8) :: delta_max, e_max, e_th, s_param
      REAL(KIND=8) :: energy_ratio
      
      ! Validate inputs
      IF (MATERIAL_ID < 1 .OR. MATERIAL_ID > N_SEE_MATERIALS) THEN
         DELTA_CORE = 0.d0
         RETURN
      END IF
      
      ! Check if SEE_MATERIALS array is allocated
      IF (.NOT. ALLOCATED(SEE_MATERIALS)) THEN
         DELTA_CORE = 0.d0
         RETURN
      END IF
      
      IF (IMPACT_ENERGY_EV <= 0.d0) THEN
         DELTA_CORE = 0.d0
         RETURN
      END IF

      ! Get material properties
      delta_max = SEE_MATERIALS(MATERIAL_ID)%DELTA_MAX
      e_max = SEE_MATERIALS(MATERIAL_ID)%E_MAX
      e_th = SEE_MATERIALS(MATERIAL_ID)%E_TH
      s_param = SEE_MATERIALS(MATERIAL_ID)%S_PARAMETER

      ! Check threshold
      IF (IMPACT_ENERGY_EV < e_th) THEN
         DELTA_CORE = 0.d0
         RETURN
      END IF

      ! Energy dependence (bounded Scholtz/Vaughan-type model).
      ! Keep DELTA_MAX as the true normal-incidence peak yield at E_MAX;
      ! angular effects are handled separately below.
      energy_ratio = IMPACT_ENERGY_EV / e_max
      
      IF (energy_ratio > 0.d0) THEN
         DELTA_CORE = delta_max * (s_param * energy_ratio) / &
                      (s_param - 1.d0 + energy_ratio**s_param)
      ELSE
         DELTA_CORE = 0.d0
      END IF

      ! Ensure reasonable bounds
      DELTA_CORE = MAX(0.d0, MIN(DELTA_CORE, 20.d0))  ! Cap at 20 for numerical stability

   END FUNCTION CALCULATE_ELECTRON_VAUGHAN_CORE

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION CALCULATE_ELECTRON_LOW_ENERGY_BRANCH -> conservative onset   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION CALCULATE_ELECTRON_LOW_ENERGY_BRANCH(MATERIAL_ID, IMPACT_ENERGY_EV) RESULT(DELTA_LOW)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: MATERIAL_ID
      REAL(KIND=8), INTENT(IN) :: IMPACT_ENERGY_EV
      REAL(KIND=8) :: DELTA_LOW
      REAL(KIND=8) :: E_TH, E_MAX, TRANSITION_ENERGY, ANCHOR_DELTA
      REAL(KIND=8) :: SHAPE_EXPONENT, ENERGY_PROGRESS

      DELTA_LOW = 0.d0
      IF (MATERIAL_ID < 1 .OR. MATERIAL_ID > N_SEE_MATERIALS) RETURN
      IF (IMPACT_ENERGY_EV <= 0.d0) RETURN

      E_TH = SEE_MATERIALS(MATERIAL_ID)%E_TH
      E_MAX = SEE_MATERIALS(MATERIAL_ID)%E_MAX
      IF (IMPACT_ENERGY_EV <= E_TH) RETURN

      TRANSITION_ENERGY = MAX(SEE_MATERIALS(MATERIAL_ID)%E1, E_TH + 5.d0)
      TRANSITION_ENERGY = MIN(TRANSITION_ENERGY, MAX(0.8d0 * E_MAX, E_TH + 5.d0))
      TRANSITION_ENERGY = MAX(TRANSITION_ENERGY, E_TH + 5.d0)

      ANCHOR_DELTA = CALCULATE_ELECTRON_VAUGHAN_CORE(MATERIAL_ID, TRANSITION_ENERGY)
      IF (ANCHOR_DELTA <= 0.d0) RETURN

      ENERGY_PROGRESS = (IMPACT_ENERGY_EV - E_TH) / MAX(TRANSITION_ENERGY - E_TH, 1.d-12)
      ENERGY_PROGRESS = MIN(MAX(ENERGY_PROGRESS, 0.d0), 1.d0)

      ! Smaller SIGMA gives a more conservative near-threshold rise.
      SHAPE_EXPONENT = 1.d0 / MAX(SEE_MATERIALS(MATERIAL_ID)%SIGMA, 0.15d0)
      DELTA_LOW = ANCHOR_DELTA * ENERGY_PROGRESS**SHAPE_EXPONENT

   END FUNCTION CALCULATE_ELECTRON_LOW_ENERGY_BRANCH

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION CALCULATE_SEE_YIELD -> segmented electron SEE with angle mod !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION CALCULATE_SEE_YIELD(MATERIAL_ID, IMPACT_ENERGY_EV, IMPACT_ANGLE_DEGREES) RESULT(DELTA)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: MATERIAL_ID
      REAL(KIND=8), INTENT(IN) :: IMPACT_ENERGY_EV, IMPACT_ANGLE_DEGREES
      REAL(KIND=8) :: DELTA
      REAL(KIND=8) :: E_TH, E_MAX, TRANSITION_ENERGY, BLEND_LOW, BLEND_HIGH, BLEND_WEIGHT
      REAL(KIND=8) :: DELTA_LOW, DELTA_CORE, ANGLE_FACTOR, ANGLE_REFERENCE_ENERGY
      
      DELTA = 0.d0

      IF (MATERIAL_ID < 1 .OR. MATERIAL_ID > N_SEE_MATERIALS) RETURN
      IF (.NOT. ALLOCATED(SEE_MATERIALS)) RETURN
      IF (IMPACT_ENERGY_EV <= 0.d0) RETURN

      E_TH = SEE_MATERIALS(MATERIAL_ID)%E_TH
      E_MAX = SEE_MATERIALS(MATERIAL_ID)%E_MAX
      IF (IMPACT_ENERGY_EV <= E_TH) RETURN

      TRANSITION_ENERGY = MAX(SEE_MATERIALS(MATERIAL_ID)%E1, E_TH + 5.d0)
      TRANSITION_ENERGY = MIN(TRANSITION_ENERGY, MAX(0.8d0 * E_MAX, E_TH + 5.d0))
      TRANSITION_ENERGY = MAX(TRANSITION_ENERGY, E_TH + 5.d0)

      DELTA_LOW = CALCULATE_ELECTRON_LOW_ENERGY_BRANCH(MATERIAL_ID, IMPACT_ENERGY_EV)
      DELTA_CORE = CALCULATE_ELECTRON_VAUGHAN_CORE(MATERIAL_ID, IMPACT_ENERGY_EV)

      BLEND_LOW = 0.8d0 * TRANSITION_ENERGY
      BLEND_HIGH = 1.2d0 * TRANSITION_ENERGY

      IF (IMPACT_ENERGY_EV <= BLEND_LOW) THEN
         DELTA = DELTA_LOW
      ELSE IF (IMPACT_ENERGY_EV >= BLEND_HIGH) THEN
         DELTA = DELTA_CORE
      ELSE
         BLEND_WEIGHT = (IMPACT_ENERGY_EV - BLEND_LOW) / MAX(BLEND_HIGH - BLEND_LOW, 1.d-12)
         BLEND_WEIGHT = BLEND_WEIGHT * BLEND_WEIGHT * (3.d0 - 2.d0 * BLEND_WEIGHT)
         DELTA = (1.d0 - BLEND_WEIGHT) * DELTA_LOW + BLEND_WEIGHT * DELTA_CORE
      END IF

      ANGLE_REFERENCE_ENERGY = MAX(SEE_MATERIALS(MATERIAL_ID)%E2, TRANSITION_ENERGY)
      ANGLE_FACTOR = CALCULATE_ENERGY_MODULATED_ANGLE_FACTOR(IMPACT_ANGLE_DEGREES, &
                     SEE_MATERIALS(MATERIAL_ID)%P1, SEE_MATERIALS(MATERIAL_ID)%P2, &
                     IMPACT_ENERGY_EV, ANGLE_REFERENCE_ENERGY, 1.5d0)

      DELTA = DELTA * ANGLE_FACTOR
      DELTA = MAX(0.d0, MIN(DELTA, 20.d0))

   END FUNCTION CALCULATE_SEE_YIELD

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION CALCULATE_ION_SEE_YIELD -> Calculate ion-induced SEE yield (γ coefficient)!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION CALCULATE_ION_SEE_YIELD(MATERIAL_ID, ION_SPECIES_ID, ION_ENERGY_EV, ION_CHARGE, ION_ANGLE_DEGREES) RESULT(GAMMA)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: MATERIAL_ID, ION_SPECIES_ID
      REAL(KIND=8), INTENT(IN) :: ION_ENERGY_EV, ION_CHARGE, ION_ANGLE_DEGREES
      REAL(KIND=8) :: GAMMA
      
      REAL(KIND=8) :: first_table_energy, anchor_gamma
      REAL(KIND=8) :: blend_start_energy, blend_weight, empirical_gamma
      INTEGER :: TABLE_INDEX
      
      ! Validate inputs
      IF (MATERIAL_ID < 1 .OR. MATERIAL_ID > N_SEE_MATERIALS) THEN
         GAMMA = 0.d0
         RETURN
      END IF
      
      ! Check if SEE_MATERIALS array is allocated
      IF (.NOT. ALLOCATED(SEE_MATERIALS)) THEN
         GAMMA = 0.d0
         RETURN
      END IF
      
      IF (ION_ENERGY_EV <= 0.d0 .OR. ION_CHARGE <= 0.d0) THEN
         GAMMA = 0.d0
         RETURN
      END IF

      TABLE_INDEX = FIND_SEE_ION_YIELD_TABLE_INDEX(MATERIAL_ID, ION_SPECIES_ID)
      IF (TABLE_INDEX > 0) THEN
         first_table_energy = SEE_ION_YIELD_TABLES(TABLE_INDEX)%ENERGY_GRID(1)
         IF (ION_ENERGY_EV >= first_table_energy) THEN
            GAMMA = LOOKUP_SEE_ION_YIELD(TABLE_INDEX, ION_ENERGY_EV, ION_ANGLE_DEGREES)
         ELSE
            empirical_gamma = CALCULATE_EMPIRICAL_ION_SEE_YIELD(MATERIAL_ID, ION_SPECIES_ID, &
                                                                ION_ENERGY_EV, ION_CHARGE, &
                                                                ION_ANGLE_DEGREES)
            anchor_gamma = LOOKUP_SEE_ION_YIELD_FIRST_ROW(TABLE_INDEX, ION_ANGLE_DEGREES)
            IF (first_table_energy > 0.d0 .AND. anchor_gamma > 0.d0) THEN
               blend_start_energy = 0.5d0 * first_table_energy
               IF (ION_ENERGY_EV <= blend_start_energy) THEN
                  GAMMA = empirical_gamma
               ELSE
                  GAMMA = CALCULATE_SUBTHRESHOLD_ION_YIELD(ION_SPECIES_ID, ION_ENERGY_EV, &
                                                           ION_ANGLE_DEGREES, first_table_energy, &
                                                           anchor_gamma)
                  blend_weight = LOG(ION_ENERGY_EV / blend_start_energy) / LOG(first_table_energy / blend_start_energy)
                  blend_weight = MIN(MAX(blend_weight, 0.d0), 1.d0)
                  GAMMA = (1.d0 - blend_weight) * empirical_gamma + blend_weight * GAMMA
               END IF
            ELSE
               GAMMA = empirical_gamma
            END IF
         END IF
         RETURN
      END IF
      GAMMA = CALCULATE_EMPIRICAL_ION_SEE_YIELD(MATERIAL_ID, ION_SPECIES_ID, ION_ENERGY_EV, &
                                                ION_CHARGE, ION_ANGLE_DEGREES)
      
   END FUNCTION CALCULATE_ION_SEE_YIELD

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION CALCULATE_NEUTRAL_SEE_YIELD -> Calculate neutral atom SEE yield!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION CALCULATE_NEUTRAL_SEE_YIELD(MATERIAL_ID, NEUTRAL_ENERGY_EV, NEUTRAL_ANGLE_DEGREES) RESULT(GAMMA)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: MATERIAL_ID
      REAL(KIND=8), INTENT(IN) :: NEUTRAL_ENERGY_EV, NEUTRAL_ANGLE_DEGREES
      REAL(KIND=8) :: GAMMA
      
      REAL(KIND=8) :: gamma_peak, peak_energy, sigma_log, onset_energy
      
      ! Validate inputs
      IF (MATERIAL_ID < 1 .OR. MATERIAL_ID > N_SEE_MATERIALS) THEN
         GAMMA = 0.d0
         RETURN
      END IF
      
      ! Check if SEE_MATERIALS array is allocated
      IF (.NOT. ALLOCATED(SEE_MATERIALS)) THEN
         GAMMA = 0.d0
         RETURN
      END IF
      
      IF (NEUTRAL_ENERGY_EV <= 0.d0) THEN
         GAMMA = 0.d0
         RETURN
      END IF
      
      ! Neutral bombardment is treated as a pure kinetic-emission process using
      ! the same peaked branch family as the high-energy part of the ion model.
      gamma_peak = MAX(SEE_MATERIALS(MATERIAL_ID)%NEUTRAL_GAMMA_MAX, 0.d0)
      peak_energy = MAX(SEE_MATERIALS(MATERIAL_ID)%NEUTRAL_E_THRESHOLD, 0.d0)
      sigma_log = MAX(SEE_MATERIALS(MATERIAL_ID)%NEUTRAL_ALPHA, 0.2d0)
      onset_energy = MAX(SEE_MATERIALS(MATERIAL_ID)%NEUTRAL_BETA, 0.d0)

      GAMMA = CALCULATE_KINETIC_SEE_BRANCH(gamma_peak, peak_energy, sigma_log, onset_energy, &
                                           NEUTRAL_ENERGY_EV, NEUTRAL_ANGLE_DEGREES, &
                                           SEE_MATERIALS(MATERIAL_ID)%NEUTRAL_ANGLE_EXP, &
                                           SEE_MATERIALS(MATERIAL_ID)%NEUTRAL_ANGLE_MAX)

      ! Physical limits (neutral yields should remain modest)
      GAMMA = MAX(0.d0, MIN(GAMMA, 1.0d0))
      
   END FUNCTION CALCULATE_NEUTRAL_SEE_YIELD

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE GENERATE_SINGLE_SECONDARY_ELECTRON -> Generate one secondary electron!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE GENERATE_SINGLE_SECONDARY_ELECTRON(X_IMPACT, Y_IMPACT, Z_IMPACT, &
                                                 SURFACE_NORMAL, MATERIAL_ID, &
                                                 CELL_ID, SECONDARY_PARTICLE)

      IMPLICIT NONE
      
      REAL(KIND=8), INTENT(IN) :: X_IMPACT, Y_IMPACT, Z_IMPACT
      REAL(KIND=8), DIMENSION(3), INTENT(IN) :: SURFACE_NORMAL
      INTEGER, INTENT(IN) :: MATERIAL_ID, CELL_ID
      TYPE(PARTICLE_DATA_STRUCTURE), INTENT(OUT) :: SECONDARY_PARTICLE
      
      REAL(KIND=8) :: energy_ev, theta, phi, speed
      REAL(KIND=8) :: cos_theta, sin_theta, cos_phi, sin_phi
      REAL(KIND=8), DIMENSION(3) :: tangent1, tangent2, velocity_vec
      REAL(KIND=8) :: w_func

      ! Validate material ID and get work function
      IF (MATERIAL_ID < 1 .OR. MATERIAL_ID > N_SEE_MATERIALS .OR. .NOT. ALLOCATED(SEE_MATERIALS)) THEN
         ! Invalid material, use default work function
         w_func = 4.5d0  ! Default work function in eV
      ELSE
         w_func = SEE_MATERIALS(MATERIAL_ID)%W_WORK_FUNCTION
      END IF

      ! Generate two orthogonal tangent vectors
      IF (ABS(SURFACE_NORMAL(1)) < 0.9d0) THEN
         tangent1 = [1.d0, 0.d0, 0.d0]
      ELSE
         tangent1 = [0.d0, 1.d0, 0.d0]
      END IF
      
      ! Make tangent1 orthogonal to normal
      tangent1 = tangent1 - DOT_PRODUCT(tangent1, SURFACE_NORMAL) * SURFACE_NORMAL
      tangent1 = tangent1 / SQRT(DOT_PRODUCT(tangent1, tangent1))
      
      ! Cross product for second tangent
      tangent2(1) = SURFACE_NORMAL(2)*tangent1(3) - SURFACE_NORMAL(3)*tangent1(2)
      tangent2(2) = SURFACE_NORMAL(3)*tangent1(1) - SURFACE_NORMAL(1)*tangent1(3)
      tangent2(3) = SURFACE_NORMAL(1)*tangent1(2) - SURFACE_NORMAL(2)*tangent1(1)

      ! Sample energy from exponential distribution (typical for SEE)
      energy_ev = w_func - 3.d0 * LOG(rf())  ! Mean energy ~3 eV above work function
      energy_ev = MAX(energy_ev, 0.1d0)  ! Minimum energy to avoid zero velocity
      
      ! Convert to speed (m/s)
      speed = SQRT(2.d0 * energy_ev * QE / ME)
      
      ! Sample emission angle (cosine distribution)
      cos_theta = SQRT(rf())  ! Cosine distribution for thermal emission
      sin_theta = SQRT(1.d0 - cos_theta**2)
      
      ! Random azimuthal angle
      phi = 2.d0 * PI * rf()
      cos_phi = COS(phi)
      sin_phi = SIN(phi)
      
      ! Construct velocity vector in surface coordinate system
      velocity_vec = speed * (cos_theta * SURFACE_NORMAL + &
                             sin_theta * cos_phi * tangent1 + &
                             sin_theta * sin_phi * tangent2)
      
      ! Create new particle
      CALL INIT_PARTICLE(X_IMPACT, Y_IMPACT, Z_IMPACT, &
                        velocity_vec(1), velocity_vec(2), velocity_vec(3), &
                        0.d0, 0.d0, SEE_ELECTRON_SPECIES_ID, CELL_ID, DT, SECONDARY_PARTICLE)

   END SUBROUTINE GENERATE_SINGLE_SECONDARY_ELECTRON

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE PROCESS_SEE_IMPACT -> Process an electron impact for SEE    !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE PROCESS_SEE_IMPACT(PARTICLE_INDEX, IMPACT_POSITION, SURFACE_NORMAL, &
                                 MATERIAL_ID, N_SECONDARY_OUT)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: PARTICLE_INDEX, MATERIAL_ID
      REAL(KIND=8), DIMENSION(3), INTENT(IN) :: IMPACT_POSITION, SURFACE_NORMAL
      INTEGER, INTENT(OUT) :: N_SECONDARY_OUT
      
      REAL(KIND=8) :: impact_energy_ev, impact_angle_degrees, delta_yield
      REAL(KIND=8) :: velocity_magnitude, cos_angle, ion_charge
      REAL(KIND=8), DIMENSION(3) :: velocity_vec
      
      N_SECONDARY_OUT = 0
      
      ! Basic validation - check if SEE is properly initialized
      IF (.NOT. BOOL_SEE_ENABLED .OR. SEE_ELECTRON_SPECIES_ID <= 0) THEN
         RETURN
      END IF
      
      ! Validate material ID early
      IF (MATERIAL_ID < 1 .OR. MATERIAL_ID > N_SEE_MATERIALS .OR. .NOT. ALLOCATED(SEE_MATERIALS)) THEN
         RETURN
      END IF
      
      ! Validate particle index
      IF (PARTICLE_INDEX < 1 .OR. PARTICLE_INDEX > NP_PROC) THEN
         RETURN
      END IF
      
      ! Check if this is an electron, ion, or neutral atom
      IF (particles(PARTICLE_INDEX)%S_ID /= SEE_ELECTRON_SPECIES_ID .AND. &
          .NOT. IS_SEE_ION(particles(PARTICLE_INDEX)%S_ID) .AND. &
          .NOT. IS_SEE_NEUTRAL(particles(PARTICLE_INDEX)%S_ID)) THEN
         RETURN
      END IF
      
      ! Calculate impact energy and angle
      velocity_vec = [particles(PARTICLE_INDEX)%VX, particles(PARTICLE_INDEX)%VY, particles(PARTICLE_INDEX)%VZ]
      velocity_magnitude = SQRT(DOT_PRODUCT(velocity_vec, velocity_vec))
      impact_energy_ev = RELATIVISTIC_KINETIC_ENERGY( &
                         SPECIES(particles(PARTICLE_INDEX)%S_ID)%MOLECULAR_MASS, &
                         particles(PARTICLE_INDEX)%VX, particles(PARTICLE_INDEX)%VY, &
                         particles(PARTICLE_INDEX)%VZ) / QE
      
      ! Calculate impact angle (angle between velocity and surface normal)
      cos_angle = -DOT_PRODUCT(velocity_vec, SURFACE_NORMAL) / velocity_magnitude
      cos_angle = MAX(-1.d0, MIN(1.d0, cos_angle))  ! Clamp to [-1,1]
      impact_angle_degrees = ACOS(ABS(cos_angle)) * 180.d0 / PI
      
      ! Count impacts for per-material statistics across all supported incident species.
      IF (MATERIAL_ID >= 1 .AND. MATERIAL_ID <= N_SEE_MATERIALS .AND. ALLOCATED(SEE_MATERIAL_IMPACTS)) THEN
         SEE_MATERIAL_IMPACTS(MATERIAL_ID) = SEE_MATERIAL_IMPACTS(MATERIAL_ID) + 1
      END IF
      IF (ALLOCATED(SEE_SPECIES_IMPACTS)) THEN
         SEE_SPECIES_IMPACTS(particles(PARTICLE_INDEX)%S_ID) = SEE_SPECIES_IMPACTS(particles(PARTICLE_INDEX)%S_ID) + 1
      END IF

      ! Calculate SEE yield based on particle type
      IF (particles(PARTICLE_INDEX)%S_ID == SEE_ELECTRON_SPECIES_ID) THEN
         ! Electron SEE (δ process)
         delta_yield = CALCULATE_SEE_YIELD(MATERIAL_ID, impact_energy_ev, impact_angle_degrees)
         
         ! Update electron statistics
         SEE_TOTAL_IMPACTS = SEE_TOTAL_IMPACTS + 1
         
      ELSE IF (IS_SEE_ION(particles(PARTICLE_INDEX)%S_ID)) THEN
         ! Ion SEE (γ process) - key for glow discharge!
         ! Works for all positive ions (H2+, H3+, etc.)
         ion_charge = ABS(SPECIES(particles(PARTICLE_INDEX)%S_ID)%CHARGE)
         delta_yield = CALCULATE_ION_SEE_YIELD(MATERIAL_ID, particles(PARTICLE_INDEX)%S_ID, &
                                               impact_energy_ev, ion_charge, impact_angle_degrees)
         
         ! Update ion statistics
         SEE_TOTAL_ION_IMPACTS = SEE_TOTAL_ION_IMPACTS + 1
         
      ELSE IF (IS_SEE_NEUTRAL(particles(PARTICLE_INDEX)%S_ID)) THEN
         ! Neutral atom SEE (pure kinetic emission)
         ! Works for all neutral species (H2, H, etc.)
         delta_yield = CALCULATE_NEUTRAL_SEE_YIELD(MATERIAL_ID, impact_energy_ev, impact_angle_degrees)
         
         ! Update neutral atom statistics
         SEE_TOTAL_NEUTRAL_IMPACTS = SEE_TOTAL_NEUTRAL_IMPACTS + 1
         
      ELSE
         ! Other particle types (should not reach here due to early return)
         delta_yield = 0.d0
      END IF
      
      ! Determine number of secondary electrons (Poisson process)
      N_SECONDARY_OUT = POISSON_SAMPLE(delta_yield)
      
      ! Update emission statistics
      IF (N_SECONDARY_OUT > 0) THEN
         SEE_TOTAL_EMISSIONS = SEE_TOTAL_EMISSIONS + N_SECONDARY_OUT
         IF (ALLOCATED(SEE_SPECIES_EMISSIONS)) THEN
            SEE_SPECIES_EMISSIONS(particles(PARTICLE_INDEX)%S_ID) = &
               SEE_SPECIES_EMISSIONS(particles(PARTICLE_INDEX)%S_ID) + N_SECONDARY_OUT
         END IF
         
         IF (IS_SEE_ION(particles(PARTICLE_INDEX)%S_ID)) THEN
            ! Ion-induced emissions
            SEE_TOTAL_ION_EMISSIONS = SEE_TOTAL_ION_EMISSIONS + N_SECONDARY_OUT
         ELSE IF (IS_SEE_NEUTRAL(particles(PARTICLE_INDEX)%S_ID)) THEN
            ! Neutral atom-induced emissions
            SEE_TOTAL_NEUTRAL_EMISSIONS = SEE_TOTAL_NEUTRAL_EMISSIONS + N_SECONDARY_OUT
         END IF
         
         IF (MATERIAL_ID >= 1 .AND. MATERIAL_ID <= N_SEE_MATERIALS .AND. ALLOCATED(SEE_MATERIAL_EMISSIONS)) THEN
            SEE_MATERIAL_EMISSIONS(MATERIAL_ID) = SEE_MATERIAL_EMISSIONS(MATERIAL_ID) + N_SECONDARY_OUT
         END IF
      END IF

   END SUBROUTINE PROCESS_SEE_IMPACT

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION POISSON_SAMPLE -> Sample from Poisson distribution           !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION POISSON_SAMPLE(LAMBDA) RESULT(K)

      IMPLICIT NONE
      
      REAL(KIND=8), INTENT(IN) :: LAMBDA
      INTEGER :: K
      
      REAL(KIND=8) :: L, P
      
      IF (LAMBDA <= 0.d0) THEN
         K = 0
         RETURN
      END IF
      
      ! For small lambda, use Knuth's algorithm
      IF (LAMBDA < 10.d0) THEN
         L = EXP(-LAMBDA)
         K = 0
         P = 1.d0
         
         DO WHILE (P > L)
            K = K + 1
            P = P * rf()
         END DO
         
         K = K - 1
      ELSE
         ! For large lambda, use normal approximation
         K = INT(LAMBDA + SQRT(LAMBDA) * GAUSSIAN_RANDOM() + 0.5d0)
         K = MAX(0, K)
      END IF

   END FUNCTION POISSON_SAMPLE

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION GAUSSIAN_RANDOM -> Generate Gaussian random number           !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION GAUSSIAN_RANDOM() RESULT(G)

      IMPLICIT NONE
      
      REAL(KIND=8) :: G
      REAL(KIND=8) :: U1, U2
      LOGICAL, SAVE :: HAVE_SPARE = .FALSE.
      REAL(KIND=8), SAVE :: SPARE
      
      IF (HAVE_SPARE) THEN
         G = SPARE
         HAVE_SPARE = .FALSE.
      ELSE
         U1 = rf()
         U2 = rf()
         G = SQRT(-2.d0 * LOG(U1)) * COS(2.d0 * PI * U2)
         SPARE = SQRT(-2.d0 * LOG(U1)) * SIN(2.d0 * PI * U2)
         HAVE_SPARE = .TRUE.
      END IF

   END FUNCTION GAUSSIAN_RANDOM

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE GATHER_SEE_STATISTICS -> Reduce SEE counters across MPI    !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE GATHER_SEE_STATISTICS(TOTAL_ELECTRON_IMPACTS, TOTAL_EMISSIONS, TOTAL_ION_IMPACTS, &
                                    TOTAL_ION_EMISSIONS, TOTAL_NEUTRAL_IMPACTS, TOTAL_NEUTRAL_EMISSIONS, &
                                    MATERIAL_IMPACTS_GLOBAL, MATERIAL_EMISSIONS_GLOBAL)

      IMPLICIT NONE

      INTEGER(KIND=8), INTENT(OUT) :: TOTAL_ELECTRON_IMPACTS, TOTAL_EMISSIONS
      INTEGER(KIND=8), INTENT(OUT) :: TOTAL_ION_IMPACTS, TOTAL_ION_EMISSIONS
      INTEGER(KIND=8), INTENT(OUT) :: TOTAL_NEUTRAL_IMPACTS, TOTAL_NEUTRAL_EMISSIONS
      INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: MATERIAL_IMPACTS_GLOBAL, MATERIAL_EMISSIONS_GLOBAL
      INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: MATERIAL_IMPACTS_LOCAL, MATERIAL_EMISSIONS_LOCAL
      INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: MATERIAL_IMPACTS_RECV, MATERIAL_EMISSIONS_RECV

      CALL MPI_REDUCE(SEE_TOTAL_IMPACTS, TOTAL_ELECTRON_IMPACTS, 1, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(SEE_TOTAL_EMISSIONS, TOTAL_EMISSIONS, 1, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(SEE_TOTAL_ION_IMPACTS, TOTAL_ION_IMPACTS, 1, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(SEE_TOTAL_ION_EMISSIONS, TOTAL_ION_EMISSIONS, 1, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(SEE_TOTAL_NEUTRAL_IMPACTS, TOTAL_NEUTRAL_IMPACTS, 1, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(SEE_TOTAL_NEUTRAL_EMISSIONS, TOTAL_NEUTRAL_EMISSIONS, 1, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

      IF (N_SEE_MATERIALS > 0) THEN
         ALLOCATE(MATERIAL_IMPACTS_LOCAL(N_SEE_MATERIALS), MATERIAL_EMISSIONS_LOCAL(N_SEE_MATERIALS))

         MATERIAL_IMPACTS_LOCAL = 0
         MATERIAL_EMISSIONS_LOCAL = 0

         IF (PROC_ID == 0) THEN
            ALLOCATE(MATERIAL_IMPACTS_GLOBAL(N_SEE_MATERIALS), MATERIAL_EMISSIONS_GLOBAL(N_SEE_MATERIALS))
            MATERIAL_IMPACTS_GLOBAL = 0
            MATERIAL_EMISSIONS_GLOBAL = 0
         ELSE
            ! The MPI receive buffer is significant only on the root rank.
            ! Keep non-root dummy buffers local so caller-owned allocatables are
            ! never allocated/deallocated on ranks that do not use the result.
            ALLOCATE(MATERIAL_IMPACTS_RECV(N_SEE_MATERIALS), MATERIAL_EMISSIONS_RECV(N_SEE_MATERIALS))
            MATERIAL_IMPACTS_RECV = 0
            MATERIAL_EMISSIONS_RECV = 0
         END IF

         IF (ALLOCATED(SEE_MATERIAL_IMPACTS)) MATERIAL_IMPACTS_LOCAL = SEE_MATERIAL_IMPACTS
         IF (ALLOCATED(SEE_MATERIAL_EMISSIONS)) MATERIAL_EMISSIONS_LOCAL = SEE_MATERIAL_EMISSIONS

         IF (PROC_ID == 0) THEN
            CALL MPI_REDUCE(MATERIAL_IMPACTS_LOCAL, MATERIAL_IMPACTS_GLOBAL, N_SEE_MATERIALS, &
                            MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_REDUCE(MATERIAL_EMISSIONS_LOCAL, MATERIAL_EMISSIONS_GLOBAL, N_SEE_MATERIALS, &
                            MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         ELSE
            CALL MPI_REDUCE(MATERIAL_IMPACTS_LOCAL, MATERIAL_IMPACTS_RECV, N_SEE_MATERIALS, &
                            MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_REDUCE(MATERIAL_EMISSIONS_LOCAL, MATERIAL_EMISSIONS_RECV, N_SEE_MATERIALS, &
                            MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         END IF

         DEALLOCATE(MATERIAL_IMPACTS_LOCAL, MATERIAL_EMISSIONS_LOCAL)
         IF (ALLOCATED(MATERIAL_IMPACTS_RECV)) DEALLOCATE(MATERIAL_IMPACTS_RECV)
         IF (ALLOCATED(MATERIAL_EMISSIONS_RECV)) DEALLOCATE(MATERIAL_EMISSIONS_RECV)
      END IF

   END SUBROUTINE GATHER_SEE_STATISTICS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE WRITE_SEE_STATISTICS -> Write SEE statistics to file       !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE WRITE_SEE_STATISTICS(TIMESTEP)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: TIMESTEP
      INTEGER :: IM, out_unit = 30
      CHARACTER(LEN=512) :: filename
      REAL(KIND=8) :: global_yield, electron_yield, ion_yield, neutral_yield
      INTEGER(KIND=8) :: TOTAL_ELECTRON_IMPACTS, TOTAL_EMISSIONS
      INTEGER(KIND=8) :: TOTAL_ION_IMPACTS, TOTAL_ION_EMISSIONS
      INTEGER(KIND=8) :: TOTAL_NEUTRAL_IMPACTS, TOTAL_NEUTRAL_EMISSIONS
      INTEGER(KIND=8) :: TOTAL_SUPPORTED_IMPACTS
      INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: MATERIAL_IMPACTS_GLOBAL, MATERIAL_EMISSIONS_GLOBAL
      
      IF (.NOT. BOOL_SEE_ENABLED) RETURN

      CALL GATHER_SEE_STATISTICS(TOTAL_ELECTRON_IMPACTS, TOTAL_EMISSIONS, TOTAL_ION_IMPACTS, TOTAL_ION_EMISSIONS, &
                                 TOTAL_NEUTRAL_IMPACTS, TOTAL_NEUTRAL_EMISSIONS, MATERIAL_IMPACTS_GLOBAL, MATERIAL_EMISSIONS_GLOBAL)

      IF (PROC_ID /= 0) THEN
         IF (ALLOCATED(MATERIAL_IMPACTS_GLOBAL)) DEALLOCATE(MATERIAL_IMPACTS_GLOBAL)
         IF (ALLOCATED(MATERIAL_EMISSIONS_GLOBAL)) DEALLOCATE(MATERIAL_EMISSIONS_GLOBAL)
         RETURN
      END IF
      
      ! Calculate global yield across all supported impacting species
      TOTAL_SUPPORTED_IMPACTS = TOTAL_ELECTRON_IMPACTS + TOTAL_ION_IMPACTS + TOTAL_NEUTRAL_IMPACTS
      IF (TOTAL_SUPPORTED_IMPACTS > 0) THEN
         global_yield = REAL(TOTAL_EMISSIONS, KIND=8) / REAL(TOTAL_SUPPORTED_IMPACTS, KIND=8)
      ELSE
         global_yield = 0.d0
      END IF
      SEE_TOTAL_YIELD = global_yield
      
      ! Write to file
      IF (TRIM(SEE_STATS_SAVE_PATH) /= '') THEN
         WRITE(filename, '(A,A,I8.8,A)') TRIM(SEE_STATS_SAVE_PATH), '/see_stats_', TIMESTEP, '.dat'
         OPEN(UNIT=out_unit, FILE=TRIM(filename), STATUS='replace', ACTION='write')
         
         WRITE(out_unit, '(A)') '# SEE Statistics'
         WRITE(out_unit, '(A,I12)') '# Timestep: ', TIMESTEP
         WRITE(out_unit, '(A,I12)') '# Total electron impacts: ', TOTAL_ELECTRON_IMPACTS
         WRITE(out_unit, '(A,I12)') '# Total ion impacts: ', TOTAL_ION_IMPACTS
         WRITE(out_unit, '(A,I12)') '# Total neutral atom impacts: ', TOTAL_NEUTRAL_IMPACTS
         WRITE(out_unit, '(A,I12)') '# Total emissions: ', TOTAL_EMISSIONS
         WRITE(out_unit, '(A,I12)') '# Ion-induced emissions: ', TOTAL_ION_EMISSIONS
         WRITE(out_unit, '(A,I12)') '# Neutral-induced emissions: ', TOTAL_NEUTRAL_EMISSIONS
         
         ! Calculate separate yields
         IF (TOTAL_ELECTRON_IMPACTS > 0) THEN
            electron_yield = REAL(TOTAL_EMISSIONS - TOTAL_ION_EMISSIONS - TOTAL_NEUTRAL_EMISSIONS, KIND=8) / &
                             REAL(TOTAL_ELECTRON_IMPACTS, KIND=8)
            WRITE(out_unit, '(A,F8.4)') '# Electron SEE yield (δ): ', electron_yield
         END IF
         
         IF (TOTAL_ION_IMPACTS > 0) THEN
            ion_yield = REAL(TOTAL_ION_EMISSIONS, KIND=8) / REAL(TOTAL_ION_IMPACTS, KIND=8)
            WRITE(out_unit, '(A,F8.4)') '# Ion SEE yield (γ): ', ion_yield
            SEE_TOTAL_ION_YIELD = ion_yield
         END IF
         
         IF (TOTAL_NEUTRAL_IMPACTS > 0) THEN
            neutral_yield = REAL(TOTAL_NEUTRAL_EMISSIONS, KIND=8) / REAL(TOTAL_NEUTRAL_IMPACTS, KIND=8)
            WRITE(out_unit, '(A,F8.4)') '# Neutral atom SEE yield: ', neutral_yield
            SEE_TOTAL_NEUTRAL_YIELD = neutral_yield
         END IF
         
         WRITE(out_unit, '(A,F8.4)') '# Global yield: ', global_yield
         WRITE(out_unit, '(A)') '#'
         WRITE(out_unit, '(A)') '# Material statistics:'
         WRITE(out_unit, '(A)') '# ID  Name                    Impacts      Emissions    Yield'
         
         DO IM = 1, N_SEE_MATERIALS
            IF (ALLOCATED(MATERIAL_IMPACTS_GLOBAL) .AND. MATERIAL_IMPACTS_GLOBAL(IM) > 0) THEN
               WRITE(out_unit, '(I4,2X,A20,2X,I12,2X,I12,2X,F8.4)') IM, &
                     SEE_MATERIALS(IM)%NAME, MATERIAL_IMPACTS_GLOBAL(IM), &
                     MATERIAL_EMISSIONS_GLOBAL(IM), &
                     REAL(MATERIAL_EMISSIONS_GLOBAL(IM), KIND=8) / REAL(MATERIAL_IMPACTS_GLOBAL(IM), KIND=8)
            ELSE
               WRITE(out_unit, '(I4,2X,A20,2X,I12,2X,I12,2X,F8.4)') IM, &
                     SEE_MATERIALS(IM)%NAME, 0, 0, 0.0
            END IF
         END DO
         
         CLOSE(out_unit)
      END IF

      IF (ALLOCATED(MATERIAL_IMPACTS_GLOBAL)) DEALLOCATE(MATERIAL_IMPACTS_GLOBAL)
      IF (ALLOCATED(MATERIAL_EMISSIONS_GLOBAL)) DEALLOCATE(MATERIAL_EMISSIONS_GLOBAL)

   END SUBROUTINE WRITE_SEE_STATISTICS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE CHECK_GLOW_DISCHARGE_CONDITIONS -> Check glow discharge physics!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE CHECK_GLOW_DISCHARGE_CONDITIONS(TIMESTEP)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: TIMESTEP
      REAL(KIND=8) :: gamma_coefficient, electron_yield, neutral_yield
      INTEGER(KIND=8) :: TOTAL_ELECTRON_IMPACTS, TOTAL_EMISSIONS
      INTEGER(KIND=8) :: TOTAL_ION_IMPACTS, TOTAL_ION_EMISSIONS
      INTEGER(KIND=8) :: TOTAL_NEUTRAL_IMPACTS, TOTAL_NEUTRAL_EMISSIONS
      INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: MATERIAL_IMPACTS_GLOBAL, MATERIAL_EMISSIONS_GLOBAL
      
      IF (.NOT. BOOL_SEE_ENABLED) RETURN

      CALL GATHER_SEE_STATISTICS(TOTAL_ELECTRON_IMPACTS, TOTAL_EMISSIONS, TOTAL_ION_IMPACTS, TOTAL_ION_EMISSIONS, &
                                 TOTAL_NEUTRAL_IMPACTS, TOTAL_NEUTRAL_EMISSIONS, MATERIAL_IMPACTS_GLOBAL, MATERIAL_EMISSIONS_GLOBAL)

      IF (PROC_ID /= 0) THEN
         IF (ALLOCATED(MATERIAL_IMPACTS_GLOBAL)) DEALLOCATE(MATERIAL_IMPACTS_GLOBAL)
         IF (ALLOCATED(MATERIAL_EMISSIONS_GLOBAL)) DEALLOCATE(MATERIAL_EMISSIONS_GLOBAL)
         RETURN
      END IF
      
      ! Calculate γ coefficient (ion SEE yield)
      IF (TOTAL_ION_IMPACTS > 0) THEN
         gamma_coefficient = REAL(TOTAL_ION_EMISSIONS, KIND=8) / REAL(TOTAL_ION_IMPACTS, KIND=8)
      ELSE
         gamma_coefficient = 0.d0
      END IF
      
      ! Calculate electron SEE yield
      IF (TOTAL_ELECTRON_IMPACTS > 0) THEN
         electron_yield = REAL(TOTAL_EMISSIONS - TOTAL_ION_EMISSIONS - TOTAL_NEUTRAL_EMISSIONS, KIND=8) / &
                          REAL(TOTAL_ELECTRON_IMPACTS, KIND=8)
      ELSE
         electron_yield = 0.d0
      END IF
      
      ! Calculate neutral atom SEE yield
      IF (TOTAL_NEUTRAL_IMPACTS > 0) THEN
         neutral_yield = REAL(TOTAL_NEUTRAL_EMISSIONS, KIND=8) / REAL(TOTAL_NEUTRAL_IMPACTS, KIND=8)
      ELSE
         neutral_yield = 0.d0
      END IF
      
      ! Output diagnostic information
      IF (MOD(TIMESTEP, 100) == 0) THEN
         WRITE(*,*) '=== Glow Discharge SEE Diagnostics (Step ', TIMESTEP, ') ==='
         WRITE(*,*) '  γ coefficient (ion SEE): ', gamma_coefficient
         WRITE(*,*) '  δ coefficient (electron SEE): ', electron_yield
         WRITE(*,*) '  Neutral atom SEE yield: ', neutral_yield
         WRITE(*,*) '  Total ion impacts: ', TOTAL_ION_IMPACTS
         WRITE(*,*) '  Total ion emissions: ', TOTAL_ION_EMISSIONS
         WRITE(*,*) '  Total electron impacts: ', TOTAL_ELECTRON_IMPACTS
         WRITE(*,*) '  Total electron emissions: ', TOTAL_EMISSIONS - TOTAL_ION_EMISSIONS - TOTAL_NEUTRAL_EMISSIONS
         WRITE(*,*) '  Total neutral atom impacts: ', TOTAL_NEUTRAL_IMPACTS
         WRITE(*,*) '  Total neutral emissions: ', TOTAL_NEUTRAL_EMISSIONS
         
         ! Check discharge maintenance condition
         IF (gamma_coefficient > 0.1d0) THEN
            WRITE(*,*) '  Ion SEE sufficient for discharge maintenance'
         ELSE
            WRITE(*,*) '  Ion SEE may be insufficient for discharge maintenance'
         END IF
         WRITE(*,*) '================================================'
      END IF

      IF (ALLOCATED(MATERIAL_IMPACTS_GLOBAL)) DEALLOCATE(MATERIAL_IMPACTS_GLOBAL)
      IF (ALLOCATED(MATERIAL_EMISSIONS_GLOBAL)) DEALLOCATE(MATERIAL_EMISSIONS_GLOBAL)
      
   END SUBROUTINE CHECK_GLOW_DISCHARGE_CONDITIONS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION FIND_BOUNDARY_PG_BY_NAME -> Find physical group by name        !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION FIND_BOUNDARY_PG_BY_NAME(BOUNDARY_NAME) RESULT(PG_ID)

      IMPLICIT NONE
      
      CHARACTER(LEN=*), INTENT(IN) :: BOUNDARY_NAME
      INTEGER :: PG_ID
      INTEGER :: IPG
      
      PG_ID = -1  ! Default: not found
      
      ! Check if GRID_BC is allocated and search through all physical groups
      IF (ALLOCATED(GRID_BC) .AND. N_GRID_BC > 0) THEN
         DO IPG = 1, MIN(N_GRID_BC, SIZE(GRID_BC))  ! Use smaller of the two for safety
            IF (TRIM(GRID_BC(IPG)%PHYSICAL_GROUP_NAME) == TRIM(BOUNDARY_NAME)) THEN
               PG_ID = IPG
               RETURN
            END IF
         END DO
      END IF

   END FUNCTION FIND_BOUNDARY_PG_BY_NAME

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION FIND_SEE_MATERIAL_FOR_BOUNDARY -> Find SEE material for boundary!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION FIND_SEE_MATERIAL_FOR_BOUNDARY(BOUNDARY_ID, WALL_ID) RESULT(MATERIAL_ID)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: BOUNDARY_ID, WALL_ID
      INTEGER :: MATERIAL_ID
      INTEGER :: IB, PG_ID
      
      MATERIAL_ID = -1  ! Default: no SEE material
      
      ! Check if SEE materials are loaded
      IF (N_SEE_MATERIALS <= 0 .OR. .NOT. ALLOCATED(SEE_MATERIALS)) THEN
         RETURN
      END IF
      
      ! Search through boundary mappings
      DO IB = 1, N_SEE_BOUNDARIES
         IF (SEE_BOUNDARY_MAP(IB)%ENABLED) THEN
            ! Check for wall-specific mapping first
            IF (WALL_ID >= 0 .AND. SEE_BOUNDARY_MAP(IB)%WALL_ID == WALL_ID) THEN
               MATERIAL_ID = SEE_BOUNDARY_MAP(IB)%MATERIAL_ID
               ! Validate material ID
               IF (MATERIAL_ID < 1 .OR. MATERIAL_ID > N_SEE_MATERIALS) THEN
                  WRITE(*,*) 'Warning: SEE material ID out of range: ', MATERIAL_ID, ' (max: ', N_SEE_MATERIALS, ')'
                  MATERIAL_ID = -1
               END IF
               RETURN
            ! Check for boundary mapping
            ELSE IF (WALL_ID < 0) THEN
               IF (SEE_BOUNDARY_MAP(IB)%USE_NAME) THEN
                  ! Find physical group by name
                  PG_ID = FIND_BOUNDARY_PG_BY_NAME(SEE_BOUNDARY_MAP(IB)%BOUNDARY_NAME)
                  IF (PG_ID == BOUNDARY_ID) THEN
                     MATERIAL_ID = SEE_BOUNDARY_MAP(IB)%MATERIAL_ID
                     ! Validate material ID
                     IF (MATERIAL_ID < 1 .OR. MATERIAL_ID > N_SEE_MATERIALS) THEN
                        WRITE(*,*) 'Warning: SEE material ID out of range: ', MATERIAL_ID, ' (max: ', N_SEE_MATERIALS, ')'
                        MATERIAL_ID = -1
                     END IF
                     RETURN
                  END IF
               ELSE
                  ! Use numeric ID
                  IF (SEE_BOUNDARY_MAP(IB)%BOUNDARY_ID == BOUNDARY_ID) THEN
                     MATERIAL_ID = SEE_BOUNDARY_MAP(IB)%MATERIAL_ID
                     ! Validate material ID
                     IF (MATERIAL_ID < 1 .OR. MATERIAL_ID > N_SEE_MATERIALS) THEN
                        WRITE(*,*) 'Warning: SEE material ID out of range: ', MATERIAL_ID, ' (max: ', N_SEE_MATERIALS, ')'
                        MATERIAL_ID = -1
                     END IF
                     RETURN
                  END IF
               END IF
            END IF
         END IF
      END DO

   END FUNCTION FIND_SEE_MATERIAL_FOR_BOUNDARY

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION IS_SEE_ION -> Check if species ID is in ion species array    !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION IS_SEE_ION(SPECIES_ID) RESULT(IS_ION)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: SPECIES_ID
      LOGICAL :: IS_ION
      INTEGER :: I
      
      IS_ION = .FALSE.
      
      ! Check if array is allocated
      IF (.NOT. ALLOCATED(SEE_ION_SPECIES_IDS) .OR. N_SEE_ION_SPECIES <= 0) RETURN
      
      ! Search for species ID in array
      DO I = 1, N_SEE_ION_SPECIES
         IF (SPECIES_ID == SEE_ION_SPECIES_IDS(I)) THEN
            IS_ION = .TRUE.
            RETURN
         END IF
      END DO

   END FUNCTION IS_SEE_ION

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION IS_SEE_NEUTRAL -> Check if species ID is in neutral array    !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION IS_SEE_NEUTRAL(SPECIES_ID) RESULT(IS_NEUTRAL)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: SPECIES_ID
      LOGICAL :: IS_NEUTRAL
      INTEGER :: I
      
      IS_NEUTRAL = .FALSE.
      
      ! Check if array is allocated
      IF (.NOT. ALLOCATED(SEE_NEUTRAL_SPECIES_IDS) .OR. N_SEE_NEUTRAL_SPECIES <= 0) RETURN
      
      ! Search for species ID in array
      DO I = 1, N_SEE_NEUTRAL_SPECIES
         IF (SPECIES_ID == SEE_NEUTRAL_SPECIES_IDS(I)) THEN
            IS_NEUTRAL = .TRUE.
            RETURN
         END IF
      END DO

   END FUNCTION IS_SEE_NEUTRAL

END MODULE secondary_electron_emission 
