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

      ! Allocate statistics arrays
      IF (N_SEE_MATERIALS > 0) THEN
         ALLOCATE(SEE_MATERIAL_IMPACTS(N_SEE_MATERIALS))
         ALLOCATE(SEE_MATERIAL_EMISSIONS(N_SEE_MATERIALS))
         SEE_MATERIAL_IMPACTS = 0
         SEE_MATERIAL_EMISSIONS = 0
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
      INTEGER :: ios, IM
      CHARACTER(LEN=512) :: line
      CHARACTER(LEN=64) :: name
      REAL(KIND=8) :: delta_max, e_max, e_th, s_param, w_func
      REAL(KIND=8) :: gamma_max, ion_e_th, ion_alpha, ion_beta, ion_charge_factor
      REAL(KIND=8) :: neutral_gamma_max, neutral_e_th, neutral_alpha, neutral_beta
      REAL(KIND=8) :: p1, p2, e1, e2, sigma
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
         
         ! Parse material definition (extended format with neutral atom SEE parameters)
         READ(line, *, iostat=ios) name, delta_max, e_max, e_th, s_param, &
                                   gamma_max, ion_e_th, ion_alpha, ion_beta, ion_charge_factor, &
                                   w_func, p1, p2, e1, e2, sigma, &
                                   neutral_gamma_max, neutral_e_th, neutral_alpha, neutral_beta
         IF (ios /= 0) THEN
            WRITE(*,*) 'Error parsing SEE material line: ', TRIM(line)
            WRITE(*,*) 'Expected 20 parameters, got error code: ', ios
            CYCLE
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
         ! Common parameters
         SEE_MATERIALS(N_SEE_MATERIALS)%W_WORK_FUNCTION = w_func
         SEE_MATERIALS(N_SEE_MATERIALS)%P1 = p1
         SEE_MATERIALS(N_SEE_MATERIALS)%P2 = p2
         SEE_MATERIALS(N_SEE_MATERIALS)%E1 = e1
         SEE_MATERIALS(N_SEE_MATERIALS)%E2 = e2
         SEE_MATERIALS(N_SEE_MATERIALS)%SIGMA = sigma

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
   ! FUNCTION CALCULATE_SEE_YIELD -> Calculate SEE yield using Vaughan model!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION CALCULATE_SEE_YIELD(MATERIAL_ID, IMPACT_ENERGY_EV, IMPACT_ANGLE_DEGREES) RESULT(DELTA)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: MATERIAL_ID
      REAL(KIND=8), INTENT(IN) :: IMPACT_ENERGY_EV, IMPACT_ANGLE_DEGREES
      REAL(KIND=8) :: DELTA
      
      REAL(KIND=8) :: delta_max, e_max, e_th, s_param
      REAL(KIND=8) :: cos_theta, f_angle, energy_ratio, x
      
      ! Validate inputs
      IF (MATERIAL_ID < 1 .OR. MATERIAL_ID > N_SEE_MATERIALS) THEN
         DELTA = 0.d0
         RETURN
      END IF
      
      ! Check if SEE_MATERIALS array is allocated
      IF (.NOT. ALLOCATED(SEE_MATERIALS)) THEN
         DELTA = 0.d0
         RETURN
      END IF
      
      IF (IMPACT_ENERGY_EV <= 0.d0) THEN
         DELTA = 0.d0
         RETURN
      END IF

      ! Get material properties
      delta_max = SEE_MATERIALS(MATERIAL_ID)%DELTA_MAX
      e_max = SEE_MATERIALS(MATERIAL_ID)%E_MAX
      e_th = SEE_MATERIALS(MATERIAL_ID)%E_TH
      s_param = SEE_MATERIALS(MATERIAL_ID)%S_PARAMETER

      ! Check threshold
      IF (IMPACT_ENERGY_EV < e_th) THEN
         DELTA = 0.d0
         RETURN
      END IF

      ! Angular dependence (Vaughan model)
      cos_theta = COS(IMPACT_ANGLE_DEGREES * PI / 180.d0)
      IF (cos_theta <= 0.d0) THEN
         DELTA = 0.d0
         RETURN
      END IF
      f_angle = 1.d0 / cos_theta

      ! Energy dependence (Vaughan model)
      energy_ratio = IMPACT_ENERGY_EV / e_max
      x = energy_ratio - 1.d0
      
      ! Vaughan formula: δ = δ_max * (s*x / (s-1+x^s)) * exp(1-x)
      IF (energy_ratio > 0.d0) THEN
         DELTA = delta_max * (s_param * energy_ratio) / (s_param - 1.d0 + energy_ratio**s_param) * &
                 EXP(1.d0 - energy_ratio) * f_angle
      ELSE
         DELTA = 0.d0
      END IF

      ! Ensure reasonable bounds
      DELTA = MAX(0.d0, MIN(DELTA, 20.d0))  ! Cap at 20 for numerical stability

   END FUNCTION CALCULATE_SEE_YIELD

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION CALCULATE_ION_SEE_YIELD -> Calculate ion-induced SEE yield (γ coefficient)!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION CALCULATE_ION_SEE_YIELD(MATERIAL_ID, ION_ENERGY_EV, ION_CHARGE, ION_ANGLE_DEGREES) RESULT(GAMMA)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: MATERIAL_ID
      REAL(KIND=8), INTENT(IN) :: ION_ENERGY_EV, ION_CHARGE, ION_ANGLE_DEGREES
      REAL(KIND=8) :: GAMMA
      
      REAL(KIND=8) :: gamma_max, e_threshold, alpha, beta, charge_factor
      REAL(KIND=8) :: energy_ratio, angle_factor, cos_angle
      
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
      
      ! Get material parameters
      gamma_max = SEE_MATERIALS(MATERIAL_ID)%GAMMA_MAX
      e_threshold = SEE_MATERIALS(MATERIAL_ID)%ION_E_THRESHOLD
      alpha = SEE_MATERIALS(MATERIAL_ID)%ION_ALPHA
      beta = SEE_MATERIALS(MATERIAL_ID)%ION_BETA
      charge_factor = SEE_MATERIALS(MATERIAL_ID)%ION_CHARGE_FACTOR
      
      ! Energy threshold check
      IF (ION_ENERGY_EV < e_threshold) THEN
         GAMMA = 0.d0
         RETURN
      END IF
      
      ! Ion SEE yield calculation (based on Hagstrum model)
      energy_ratio = ION_ENERGY_EV / e_threshold
      GAMMA = gamma_max * (energy_ratio**alpha) * EXP(-beta/SQRT(ION_ENERGY_EV))
      
      ! Charge state correction
      GAMMA = GAMMA * (ION_CHARGE * charge_factor)
      
      ! Angle correction (ion SEE has weaker angle dependence)
      cos_angle = COS(ION_ANGLE_DEGREES * PI / 180.d0)
      angle_factor = 1.d0 + 0.1d0 * (1.d0 - cos_angle)  ! Simplified angle correction
      GAMMA = GAMMA * angle_factor
      
      ! Physical limits
      GAMMA = MAX(0.d0, MIN(GAMMA, 10.d0))  ! Limit to reasonable range
      
   END FUNCTION CALCULATE_ION_SEE_YIELD

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION CALCULATE_NEUTRAL_SEE_YIELD -> Calculate neutral atom SEE yield!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   FUNCTION CALCULATE_NEUTRAL_SEE_YIELD(MATERIAL_ID, NEUTRAL_ENERGY_EV, NEUTRAL_ANGLE_DEGREES) RESULT(GAMMA)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: MATERIAL_ID
      REAL(KIND=8), INTENT(IN) :: NEUTRAL_ENERGY_EV, NEUTRAL_ANGLE_DEGREES
      REAL(KIND=8) :: GAMMA
      
      REAL(KIND=8) :: gamma_max, e_threshold, alpha, beta
      REAL(KIND=8) :: energy_ratio, angle_factor, cos_angle
      
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
      
      ! Get material parameters for neutral atoms
      gamma_max = SEE_MATERIALS(MATERIAL_ID)%NEUTRAL_GAMMA_MAX
      e_threshold = SEE_MATERIALS(MATERIAL_ID)%NEUTRAL_E_THRESHOLD
      alpha = SEE_MATERIALS(MATERIAL_ID)%NEUTRAL_ALPHA
      beta = SEE_MATERIALS(MATERIAL_ID)%NEUTRAL_BETA
      
      ! Energy threshold check
      IF (NEUTRAL_ENERGY_EV < e_threshold) THEN
         GAMMA = 0.d0
         RETURN
      END IF
      
      ! Neutral atom SEE yield calculation (pure kinetic emission)
      ! No charge state correction - this is the key difference from ions
      energy_ratio = NEUTRAL_ENERGY_EV / e_threshold
      GAMMA = gamma_max * (energy_ratio**alpha) * EXP(-beta/SQRT(NEUTRAL_ENERGY_EV))
      
      ! Angle correction (weaker than electron, similar to ion)
      cos_angle = COS(NEUTRAL_ANGLE_DEGREES * PI / 180.d0)
      angle_factor = 1.d0 + 0.05d0 * (1.d0 - cos_angle)  ! Even weaker angle dependence
      GAMMA = GAMMA * angle_factor
      
      ! Physical limits (neutral yields should be much lower than ions)
      GAMMA = MAX(0.d0, MIN(GAMMA, 1.0d0))  ! Cap at 1.0 for neutrals
      
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
      impact_energy_ev = 0.5d0 * SPECIES(particles(PARTICLE_INDEX)%S_ID)%MOLECULAR_MASS * velocity_magnitude**2 / QE
      
      ! Calculate impact angle (angle between velocity and surface normal)
      cos_angle = -DOT_PRODUCT(velocity_vec, SURFACE_NORMAL) / velocity_magnitude
      cos_angle = MAX(-1.d0, MIN(1.d0, cos_angle))  ! Clamp to [-1,1]
      impact_angle_degrees = ACOS(ABS(cos_angle)) * 180.d0 / PI
      
      ! Calculate SEE yield based on particle type
      IF (particles(PARTICLE_INDEX)%S_ID == SEE_ELECTRON_SPECIES_ID) THEN
         ! Electron SEE (δ process)
         delta_yield = CALCULATE_SEE_YIELD(MATERIAL_ID, impact_energy_ev, impact_angle_degrees)
         
         ! Update electron statistics
         SEE_TOTAL_IMPACTS = SEE_TOTAL_IMPACTS + 1
         IF (MATERIAL_ID >= 1 .AND. MATERIAL_ID <= N_SEE_MATERIALS .AND. ALLOCATED(SEE_MATERIAL_IMPACTS)) THEN
            SEE_MATERIAL_IMPACTS(MATERIAL_ID) = SEE_MATERIAL_IMPACTS(MATERIAL_ID) + 1
         END IF
         
      ELSE IF (IS_SEE_ION(particles(PARTICLE_INDEX)%S_ID)) THEN
         ! Ion SEE (γ process) - key for glow discharge!
         ! Works for all positive ions (H2+, H3+, etc.)
         ion_charge = ABS(SPECIES(particles(PARTICLE_INDEX)%S_ID)%CHARGE)
         delta_yield = CALCULATE_ION_SEE_YIELD(MATERIAL_ID, impact_energy_ev, ion_charge, impact_angle_degrees)
         
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
   ! SUBROUTINE WRITE_SEE_STATISTICS -> Write SEE statistics to file       !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE WRITE_SEE_STATISTICS(TIMESTEP)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: TIMESTEP
      INTEGER :: IM, out_unit = 30
      CHARACTER(LEN=512) :: filename
      REAL(KIND=8) :: global_yield, electron_yield, ion_yield, neutral_yield
      
      IF (.NOT. BOOL_SEE_ENABLED .OR. PROC_ID /= 0) RETURN
      
      ! Calculate global yield
      IF (SEE_TOTAL_IMPACTS > 0) THEN
         global_yield = REAL(SEE_TOTAL_EMISSIONS, KIND=8) / REAL(SEE_TOTAL_IMPACTS, KIND=8)
      ELSE
         global_yield = 0.d0
      END IF
      
      ! Write to file
      IF (TRIM(SEE_STATS_SAVE_PATH) /= '') THEN
         WRITE(filename, '(A,A,I8.8,A)') TRIM(SEE_STATS_SAVE_PATH), '/see_stats_', TIMESTEP, '.dat'
         OPEN(UNIT=out_unit, FILE=TRIM(filename), STATUS='replace', ACTION='write')
         
         WRITE(out_unit, '(A)') '# SEE Statistics'
         WRITE(out_unit, '(A,I12)') '# Timestep: ', TIMESTEP
         WRITE(out_unit, '(A,I12)') '# Total electron impacts: ', SEE_TOTAL_IMPACTS
         WRITE(out_unit, '(A,I12)') '# Total ion impacts: ', SEE_TOTAL_ION_IMPACTS
         WRITE(out_unit, '(A,I12)') '# Total neutral atom impacts: ', SEE_TOTAL_NEUTRAL_IMPACTS
         WRITE(out_unit, '(A,I12)') '# Total emissions: ', SEE_TOTAL_EMISSIONS
         WRITE(out_unit, '(A,I12)') '# Ion-induced emissions: ', SEE_TOTAL_ION_EMISSIONS
         WRITE(out_unit, '(A,I12)') '# Neutral-induced emissions: ', SEE_TOTAL_NEUTRAL_EMISSIONS
         
         ! Calculate separate yields
         IF (SEE_TOTAL_IMPACTS > 0) THEN
            electron_yield = REAL(SEE_TOTAL_EMISSIONS - SEE_TOTAL_ION_EMISSIONS, KIND=8) / REAL(SEE_TOTAL_IMPACTS, KIND=8)
            WRITE(out_unit, '(A,F8.4)') '# Electron SEE yield (δ): ', electron_yield
         END IF
         
         IF (SEE_TOTAL_ION_IMPACTS > 0) THEN
            ion_yield = REAL(SEE_TOTAL_ION_EMISSIONS, KIND=8) / REAL(SEE_TOTAL_ION_IMPACTS, KIND=8)
            WRITE(out_unit, '(A,F8.4)') '# Ion SEE yield (γ): ', ion_yield
            SEE_TOTAL_ION_YIELD = ion_yield
         END IF
         
         IF (SEE_TOTAL_NEUTRAL_IMPACTS > 0) THEN
            neutral_yield = REAL(SEE_TOTAL_NEUTRAL_EMISSIONS, KIND=8) / REAL(SEE_TOTAL_NEUTRAL_IMPACTS, KIND=8)
            WRITE(out_unit, '(A,F8.4)') '# Neutral atom SEE yield: ', neutral_yield
            SEE_TOTAL_NEUTRAL_YIELD = neutral_yield
         END IF
         
         WRITE(out_unit, '(A,F8.4)') '# Global yield: ', global_yield
         WRITE(out_unit, '(A)') '#'
         WRITE(out_unit, '(A)') '# Material statistics:'
         WRITE(out_unit, '(A)') '# ID  Name                    Impacts      Emissions    Yield'
         
         DO IM = 1, N_SEE_MATERIALS
            IF (SEE_MATERIAL_IMPACTS(IM) > 0) THEN
               WRITE(out_unit, '(I4,2X,A20,2X,I12,2X,I12,2X,F8.4)') IM, &
                     SEE_MATERIALS(IM)%NAME, SEE_MATERIAL_IMPACTS(IM), &
                     SEE_MATERIAL_EMISSIONS(IM), &
                     REAL(SEE_MATERIAL_EMISSIONS(IM), KIND=8) / REAL(SEE_MATERIAL_IMPACTS(IM), KIND=8)
            ELSE
               WRITE(out_unit, '(I4,2X,A20,2X,I12,2X,I12,2X,F8.4)') IM, &
                     SEE_MATERIALS(IM)%NAME, 0, 0, 0.0
            END IF
         END DO
         
         CLOSE(out_unit)
      END IF

   END SUBROUTINE WRITE_SEE_STATISTICS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE CHECK_GLOW_DISCHARGE_CONDITIONS -> Check glow discharge physics!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE CHECK_GLOW_DISCHARGE_CONDITIONS(TIMESTEP)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: TIMESTEP
      REAL(KIND=8) :: gamma_coefficient, electron_yield, neutral_yield
      
      IF (.NOT. BOOL_SEE_ENABLED .OR. PROC_ID /= 0) RETURN
      
      ! Calculate γ coefficient (ion SEE yield)
      IF (SEE_TOTAL_ION_IMPACTS > 0) THEN
         gamma_coefficient = REAL(SEE_TOTAL_ION_EMISSIONS, KIND=8) / REAL(SEE_TOTAL_ION_IMPACTS, KIND=8)
      ELSE
         gamma_coefficient = 0.d0
      END IF
      
      ! Calculate electron SEE yield
      IF (SEE_TOTAL_IMPACTS > 0) THEN
         electron_yield = REAL(SEE_TOTAL_EMISSIONS - SEE_TOTAL_ION_EMISSIONS - &
                               SEE_TOTAL_NEUTRAL_EMISSIONS, KIND=8) / REAL(SEE_TOTAL_IMPACTS, KIND=8)
      ELSE
         electron_yield = 0.d0
      END IF
      
      ! Calculate neutral atom SEE yield
      IF (SEE_TOTAL_NEUTRAL_IMPACTS > 0) THEN
         neutral_yield = REAL(SEE_TOTAL_NEUTRAL_EMISSIONS, KIND=8) / REAL(SEE_TOTAL_NEUTRAL_IMPACTS, KIND=8)
      ELSE
         neutral_yield = 0.d0
      END IF
      
      ! Output diagnostic information
      IF (MOD(TIMESTEP, 100) == 0) THEN
         WRITE(*,*) '=== Glow Discharge SEE Diagnostics (Step ', TIMESTEP, ') ==='
         WRITE(*,*) '  γ coefficient (ion SEE): ', gamma_coefficient
         WRITE(*,*) '  δ coefficient (electron SEE): ', electron_yield
         WRITE(*,*) '  Neutral atom SEE yield: ', neutral_yield
         WRITE(*,*) '  Total ion impacts: ', SEE_TOTAL_ION_IMPACTS
         WRITE(*,*) '  Total ion emissions: ', SEE_TOTAL_ION_EMISSIONS
         WRITE(*,*) '  Total electron impacts: ', SEE_TOTAL_IMPACTS
         WRITE(*,*) '  Total electron emissions: ', &
                    SEE_TOTAL_EMISSIONS - SEE_TOTAL_ION_EMISSIONS - SEE_TOTAL_NEUTRAL_EMISSIONS
         WRITE(*,*) '  Total neutral atom impacts: ', SEE_TOTAL_NEUTRAL_IMPACTS
         WRITE(*,*) '  Total neutral emissions: ', SEE_TOTAL_NEUTRAL_EMISSIONS
         
         ! Check discharge maintenance condition
         IF (gamma_coefficient > 0.1d0) THEN
            WRITE(*,*) '  ✓ Ion SEE sufficient for discharge maintenance'
         ELSE
            WRITE(*,*) '  ⚠ Ion SEE may be insufficient for discharge maintenance'
         END IF
         WRITE(*,*) '================================================'
      END IF
      
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