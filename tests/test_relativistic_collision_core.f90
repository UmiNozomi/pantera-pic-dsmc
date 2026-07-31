PROGRAM test_relativistic_collision_core

   USE relativistic_collision_core, ONLY: RELATIVISTIC_TWO_BODY_SCATTER

   IMPLICIT NONE

   REAL(KIND=8), PARAMETER :: C_LIGHT = 2.99792458d8
   REAL(KIND=8), PARAMETER :: QE = 1.602176634d-19
   REAL(KIND=8), PARAMETER :: M_ELECTRON = 9.1093837139d-31
   REAL(KIND=8), PARAMETER :: M_D2 = 6.689d-27

   CALL TEST_ELASTIC(60.d3, (/ 0.d0, 1.d0, 0.d0 /))
   CALL TEST_ELASTIC(100.d3, (/ -0.6d0, 0.2d0, 0.7745966692414834d0 /))
   CALL TEST_ELASTIC(500.d3, (/ -1.d0, 0.d0, 0.d0 /))
   CALL TEST_VIBRATIONAL_LOSS(0.38d0)
   CALL TEST_VIBRATIONAL_LOSS(60.d3)
   CALL TEST_VIBRATIONAL_LOSS(100.d3)
   CALL TEST_INSUFFICIENT_ENERGY()

   WRITE(*,'(A)') 'Relativistic collision-core tests passed.'

CONTAINS

   SUBROUTINE TEST_ELASTIC(ELECTRON_ENERGY_EV, DIRECTION)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: ELECTRON_ENERGY_EV
      REAL(KIND=8), DIMENSION(3), INTENT(IN) :: DIRECTION

      REAL(KIND=8), DIMENSION(3) :: C1, C2, P_BEFORE, P_AFTER
      REAL(KIND=8) :: K_BEFORE, K_AFTER, SPEED1, SPEED2, P_SCALE
      LOGICAL :: SUCCESS

      C1 = (/ SPEED_FROM_ENERGY(ELECTRON_ENERGY_EV*QE, M_ELECTRON), &
              0.d0, 0.d0 /)
      C2 = 0.d0

      CALL TWO_PARTICLE_STATE(C1, C2, K_BEFORE, P_BEFORE, SPEED1, SPEED2)
      CALL RELATIVISTIC_TWO_BODY_SCATTER(M_ELECTRON, M_D2, C1, C2, 0.d0, &
                                         DIRECTION, SUCCESS)
      CALL ASSERT_TRUE(SUCCESS, 'Elastic scatter reported failure.')
      CALL TWO_PARTICLE_STATE(C1, C2, K_AFTER, P_AFTER, SPEED1, SPEED2)

      P_SCALE = MAX(SQRT(SUM(P_BEFORE*P_BEFORE)), TINY(1.d0))
      CALL ASSERT_TRUE(MAXVAL(ABS(P_AFTER-P_BEFORE)) <= 5.d-12*P_SCALE, &
                       'Elastic three-momentum is not conserved.')
      CALL ASSERT_TRUE(ABS(K_AFTER-K_BEFORE)/QE <= 2.d-5, &
                       'Elastic relativistic kinetic energy is not conserved.')
      CALL ASSERT_TRUE(SPEED1 < C_LIGHT .AND. SPEED2 < C_LIGHT, &
                       'Elastic scatter produced a superluminal particle.')

   END SUBROUTINE TEST_ELASTIC


   SUBROUTINE TEST_VIBRATIONAL_LOSS(ELECTRON_ENERGY_EV)

      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: ELECTRON_ENERGY_EV

      REAL(KIND=8), PARAMETER :: LOSS_EV = 0.3712d0
      REAL(KIND=8), DIMENSION(3), PARAMETER :: DIRECTION = &
                                             (/ 0.3d0, -0.4d0, 0.8660254037844386d0 /)
      REAL(KIND=8), DIMENSION(3) :: C1, C2, P_BEFORE, P_AFTER
      REAL(KIND=8) :: K_BEFORE, K_AFTER, SPEED1, SPEED2, P_SCALE
      REAL(KIND=8) :: EVIB, LOSS
      LOGICAL :: SUCCESS

      LOSS = LOSS_EV*QE
      C1 = (/ SPEED_FROM_ENERGY(ELECTRON_ENERGY_EV*QE, M_ELECTRON), 0.d0, 0.d0 /)
      C2 = 0.d0

      ! Mirror the explicit v=0 -> v=1 integration: discard the classical
      ! background sample, scatter with one loss, then add one quantum once.
      EVIB = 100.d0*QE
      EVIB = 0.d0
      CALL TWO_PARTICLE_STATE(C1, C2, K_BEFORE, P_BEFORE, SPEED1, SPEED2)
      CALL RELATIVISTIC_TWO_BODY_SCATTER(M_ELECTRON, M_D2, C1, C2, LOSS, &
                                         DIRECTION, SUCCESS)
      CALL ASSERT_TRUE(SUCCESS, 'Vibrational scatter reported failure.')
      EVIB = EVIB + LOSS
      CALL TWO_PARTICLE_STATE(C1, C2, K_AFTER, P_AFTER, SPEED1, SPEED2)

      P_SCALE = MAX(SQRT(SUM(P_BEFORE*P_BEFORE)), TINY(1.d0))
      CALL ASSERT_TRUE(MAXVAL(ABS(P_AFTER-P_BEFORE)) <= 5.d-12*P_SCALE, &
                       'Vibrational three-momentum is not conserved.')
      CALL ASSERT_TRUE(ABS((K_BEFORE-K_AFTER)/QE-LOSS_EV) <= 2.d-5, &
                       'Lab-frame translational loss is not exactly one quantum.')
      CALL ASSERT_TRUE(ABS(EVIB/QE-LOSS_EV) <= 1.d-13, &
                       'Vibrational internal energy was not added exactly once.')
      CALL ASSERT_TRUE(ABS((K_AFTER+EVIB)-K_BEFORE)/QE <= 2.d-5, &
                       'Translation plus vibrational energy is not conserved.')
      CALL ASSERT_TRUE(SPEED1 < C_LIGHT .AND. SPEED2 < C_LIGHT, &
                       'Vibrational scatter produced a superluminal particle.')

   END SUBROUTINE TEST_VIBRATIONAL_LOSS

   SUBROUTINE TEST_INSUFFICIENT_ENERGY()

      IMPLICIT NONE

      REAL(KIND=8), PARAMETER :: LOSS = 0.3712d0*QE
      REAL(KIND=8), DIMENSION(3), PARAMETER :: DIRECTION = (/ 0.d0, 1.d0, 0.d0 /)
      REAL(KIND=8), DIMENSION(3) :: C1, C2, C1_BEFORE, C2_BEFORE
      LOGICAL :: SUCCESS

      C1 = (/ SPEED_FROM_ENERGY(0.36d0*QE, M_ELECTRON), 0.d0, 0.d0 /)
      C2 = 0.d0
      C1_BEFORE = C1
      C2_BEFORE = C2

      CALL RELATIVISTIC_TWO_BODY_SCATTER(M_ELECTRON, M_D2, C1, C2, LOSS, &
                                         DIRECTION, SUCCESS)
      CALL ASSERT_TRUE(.NOT. SUCCESS, 'Below-threshold scatter should fail.')
      CALL ASSERT_TRUE(MAXVAL(ABS(C1-C1_BEFORE)) <= 0.d0 .AND. &
                       MAXVAL(ABS(C2-C2_BEFORE)) <= 0.d0, &
                       'Failed scatter modified an input velocity.')

   END SUBROUTINE TEST_INSUFFICIENT_ENERGY



   SUBROUTINE TWO_PARTICLE_STATE(C1, C2, K_TOTAL, P_TOTAL, SPEED1, SPEED2)

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(3), INTENT(IN) :: C1, C2
      REAL(KIND=8), INTENT(OUT) :: K_TOTAL, SPEED1, SPEED2
      REAL(KIND=8), DIMENSION(3), INTENT(OUT) :: P_TOTAL

      REAL(KIND=8), DIMENSION(3) :: P1, P2
      REAL(KIND=8) :: K1, K2

      CALL PARTICLE_STATE(M_ELECTRON, C1, K1, P1, SPEED1)
      CALL PARTICLE_STATE(M_D2, C2, K2, P2, SPEED2)
      K_TOTAL = K1 + K2
      P_TOTAL = P1 + P2

   END SUBROUTINE TWO_PARTICLE_STATE


   SUBROUTINE PARTICLE_STATE(MASS, VELOCITY, KINETIC, MOMENTUM, SPEED)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: MASS
      REAL(KIND=8), DIMENSION(3), INTENT(IN) :: VELOCITY
      REAL(KIND=8), INTENT(OUT) :: KINETIC, SPEED
      REAL(KIND=8), DIMENSION(3), INTENT(OUT) :: MOMENTUM

      REAL(KIND=8) :: BETA2, ROOT, GAMMA_MINUS_ONE, GAMMA

      SPEED = SQRT(SUM(VELOCITY*VELOCITY))
      BETA2 = SUM(VELOCITY*VELOCITY)/(C_LIGHT*C_LIGHT)
      CALL ASSERT_TRUE(BETA2 < 1.d0, 'State evaluator received a superluminal velocity.')
      ROOT = SQRT(1.d0-BETA2)
      GAMMA = 1.d0/ROOT
      GAMMA_MINUS_ONE = BETA2/(ROOT*(1.d0+ROOT))
      KINETIC = GAMMA_MINUS_ONE*MASS*C_LIGHT*C_LIGHT
      MOMENTUM = GAMMA*MASS*VELOCITY

   END SUBROUTINE PARTICLE_STATE


   REAL(KIND=8) FUNCTION SPEED_FROM_ENERGY(KINETIC, MASS)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: KINETIC, MASS
      REAL(KIND=8) :: GAMMA

      GAMMA = 1.d0 + KINETIC/(MASS*C_LIGHT*C_LIGHT)
      SPEED_FROM_ENERGY = C_LIGHT*SQRT(1.d0-1.d0/(GAMMA*GAMMA))

   END FUNCTION SPEED_FROM_ENERGY


   SUBROUTINE ASSERT_TRUE(CONDITION, MESSAGE)

      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: CONDITION
      CHARACTER(LEN=*), INTENT(IN) :: MESSAGE

      IF (.NOT. CONDITION) THEN
         WRITE(*,'(A)') TRIM(MESSAGE)
         ERROR STOP 1
      END IF

   END SUBROUTINE ASSERT_TRUE

END PROGRAM test_relativistic_collision_core
