PROGRAM test_boundary_control_core

   USE boundary_control_core
   IMPLICIT NONE

   CALL CHECK_CHARGE('ion absorption', 1.d0, 1.d0, 1.d0, 1.d0, &
      .FALSE., 0, 1.d0)
   CALL CHECK_CHARGE('electron absorption', -1.d0, 1.d0, -1.d0, 1.d0, &
      .FALSE., 0, -1.d0)
   CALL CHECK_CHARGE('ion reflection', 1.d0, 1.d0, 1.d0, 1.d0, &
      .TRUE., 0, 0.d0)
   CALL CHECK_CHARGE('electron reflection', -1.d0, 1.d0, -1.d0, 1.d0, &
      .TRUE., 0, 0.d0)
   CALL CHECK_CHARGE('ion absorption plus two SEE', 1.d0, 1.d0, &
      1.d0, 1.d0, .FALSE., 2, 3.d0)
   CALL CHECK_CHARGE('electron absorption plus two SEE', -1.d0, 1.d0, &
      -1.d0, 1.d0, .FALSE., 2, 1.d0)
   CALL CHECK_CHARGE('electron reflection plus two SEE', -1.d0, 1.d0, &
      -1.d0, 1.d0, .TRUE., 2, 2.d0)
   CALL CHECK_CHARGE('ion reflection plus two SEE', 1.d0, 1.d0, &
      1.d0, 1.d0, .TRUE., 2, 2.d0)
   CALL CHECK_CHARGE('ion to neutral', 1.d0, 1.d0, 0.d0, 1.d0, &
      .TRUE., 0, 1.d0)
   CALL CHECK_CHARGE('neutral to ion', 0.d0, 1.d0, 1.d0, 1.d0, &
      .TRUE., 0, -1.d0)
   CALL CHECK_CHARGE('neutral reflection plus SEE', 0.d0, 1.d0, &
      0.d0, 1.d0, .TRUE., 1, 1.d0)
   CALL CHECK_CHARGE('neutral absorption', 0.d0, 1.d0, 0.d0, 1.d0, &
      .FALSE., 0, 0.d0)
   CALL CHECK_CHARGE('neutral absorption plus SEE', 0.d0, 1.d0, &
      0.d0, 1.d0, .FALSE., 1, 1.d0)
   CALL CHECK_CHARGE('ion weight conversion', 1.d0, 2.d0, 1.d0, 1.d0, &
      .TRUE., 0, 1.d0)
   CALL CHECK_CHARGE('electron weight conversion', -1.d0, 2.d0, &
      -1.d0, 1.d0, .TRUE., 0, -1.d0)

   CALL CHECK_COMPONENTS('ion absorption plus SEE components', &
      1.d0, 1.d0, 1.d0, 1.d0, .FALSE., 2, 1.d0, 0.d0, 2.d0)
   CALL CHECK_COMPONENTS('electron absorption plus SEE components', &
      -1.d0, 1.d0, -1.d0, 1.d0, .FALSE., 2, 0.d0, -1.d0, 2.d0)
   CALL CHECK_COMPONENTS('electron reflection plus SEE components', &
      -1.d0, 1.d0, -1.d0, 1.d0, .TRUE., 2, 0.d0, 0.d0, 2.d0)
   CALL CHECK_COMPONENTS('ion to neutral components', &
      1.d0, 1.d0, 0.d0, 1.d0, .TRUE., 0, 1.d0, 0.d0, 0.d0)
   CALL CHECK_COMPONENTS('neutral to ion components', &
      0.d0, 1.d0, 1.d0, 1.d0, .TRUE., 0, -1.d0, 0.d0, 0.d0)

   CALL CHECK_HIGH_ENERGY_REFLECTION_AVERAGE
   CALL CHECK_SMOOTH_STARTUP_RAMP
   CALL CHECK_ERROR_DRIVEN_PID

   WRITE(*,'(A)') 'PASS: boundary charge, startup ramp and error-driven PID tests'

CONTAINS

   SUBROUTINE CHECK_CHARGE(NAME, QIN, WIN, QOUT, WOUT, SURVIVES, &
                           NSEE, EXPECTED)

      CHARACTER(LEN=*), INTENT(IN) :: NAME
      REAL(KIND=8), INTENT(IN) :: QIN, WIN, QOUT, WOUT, EXPECTED
      LOGICAL, INTENT(IN) :: SURVIVES
      INTEGER, INTENT(IN) :: NSEE
      REAL(KIND=8) :: QION, QELEC, QSEE

      CALL SURFACE_EVENT_CHARGE_UNITS(QIN, WIN, QOUT, WOUT, SURVIVES, &
         NSEE, -1.d0, 1.d0, QION, QELEC, QSEE)
      CALL ASSERT_CLOSE(NAME, QION + QELEC + QSEE, EXPECTED)

   END SUBROUTINE CHECK_CHARGE

   SUBROUTINE CHECK_COMPONENTS(NAME, QIN, WIN, QOUT, WOUT, SURVIVES, &
                               NSEE, EXPECTED_ION, EXPECTED_ELEC, &
                               EXPECTED_SEE)

      CHARACTER(LEN=*), INTENT(IN) :: NAME
      REAL(KIND=8), INTENT(IN) :: QIN, WIN, QOUT, WOUT
      REAL(KIND=8), INTENT(IN) :: EXPECTED_ION, EXPECTED_ELEC, EXPECTED_SEE
      LOGICAL, INTENT(IN) :: SURVIVES
      INTEGER, INTENT(IN) :: NSEE
      REAL(KIND=8) :: QION, QELEC, QSEE

      CALL SURFACE_EVENT_CHARGE_UNITS(QIN, WIN, QOUT, WOUT, SURVIVES, &
         NSEE, -1.d0, 1.d0, QION, QELEC, QSEE)
      CALL ASSERT_CLOSE(TRIM(NAME)//' ion', QION, EXPECTED_ION)
      CALL ASSERT_CLOSE(TRIM(NAME)//' electron', QELEC, EXPECTED_ELEC)
      CALL ASSERT_CLOSE(TRIM(NAME)//' SEE', QSEE, EXPECTED_SEE)

   END SUBROUTINE CHECK_COMPONENTS

   SUBROUTINE CHECK_HIGH_ENERGY_REFLECTION_AVERAGE

      REAL(KIND=8) :: QION, QELEC, QSEE
      REAL(KIND=8) :: QABS, QREFL

      CALL SURFACE_EVENT_CHARGE_UNITS(-1.d0, 1.d0, -1.d0, 1.d0, &
         .FALSE., 0, -1.d0, 1.d0, QION, QELEC, QSEE)
      QABS = QION + QELEC + QSEE
      CALL SURFACE_EVENT_CHARGE_UNITS(-1.d0, 1.d0, -1.d0, 1.d0, &
         .TRUE., 0, -1.d0, 1.d0, QION, QELEC, QSEE)
      QREFL = QION + QELEC + QSEE
      CALL ASSERT_CLOSE('25 percent absorption above energy threshold', &
         0.25d0*QABS + 0.75d0*QREFL, -0.25d0)

      CALL SURFACE_EVENT_CHARGE_UNITS(-1.d0, 1.d0, -1.d0, 1.d0, &
         .FALSE., 2, -1.d0, 1.d0, QION, QELEC, QSEE)
      QABS = QION + QELEC + QSEE
      CALL SURFACE_EVENT_CHARGE_UNITS(-1.d0, 1.d0, -1.d0, 1.d0, &
         .TRUE., 2, -1.d0, 1.d0, QION, QELEC, QSEE)
      QREFL = QION + QELEC + QSEE
      CALL ASSERT_CLOSE('high-energy absorption plus two SEE', &
         0.25d0*QABS + 0.75d0*QREFL, 1.75d0)

   END SUBROUTINE CHECK_HIGH_ENERGY_REFLECTION_AVERAGE

   SUBROUTINE CHECK_SMOOTH_STARTUP_RAMP

      REAL(KIND=8), PARAMETER :: START_VOLTAGE = 0.d0
      REAL(KIND=8), PARAMETER :: TARGET_VOLTAGE = -3.d4
      REAL(KIND=8), PARAMETER :: DURATION = 1.d-7

      CALL ASSERT_CLOSE('startup ramp before zero', &
         SMOOTH_STARTUP_VOLTAGE(-1.d-9, DURATION, START_VOLTAGE, &
                                TARGET_VOLTAGE), START_VOLTAGE)
      CALL ASSERT_CLOSE('startup ramp begins at configured voltage', &
         SMOOTH_STARTUP_VOLTAGE(0.d0, DURATION, START_VOLTAGE, &
                                TARGET_VOLTAGE), START_VOLTAGE)
      CALL ASSERT_CLOSE('startup ramp smoothstep midpoint', &
         SMOOTH_STARTUP_VOLTAGE(0.5d0*DURATION, DURATION, START_VOLTAGE, &
                                TARGET_VOLTAGE), -1.5d4)
      CALL ASSERT_CLOSE('startup ramp reaches CV voltage', &
         SMOOTH_STARTUP_VOLTAGE(DURATION, DURATION, START_VOLTAGE, &
                                TARGET_VOLTAGE), TARGET_VOLTAGE)
      CALL ASSERT_CLOSE('startup ramp remains at CV voltage', &
         SMOOTH_STARTUP_VOLTAGE(2.d0*DURATION, DURATION, START_VOLTAGE, &
                                TARGET_VOLTAGE), TARGET_VOLTAGE)
      CALL ASSERT_CLOSE('disabled ramp immediately reaches target', &
         SMOOTH_STARTUP_VOLTAGE(0.d0, 0.d0, START_VOLTAGE, &
                                TARGET_VOLTAGE), TARGET_VOLTAGE)

   END SUBROUTINE CHECK_SMOOTH_STARTUP_RAMP

   SUBROUTINE CHECK_ERROR_DRIVEN_PID

      REAL(KIND=8) :: INTEGRAL, PREV_ERROR, PREVIOUS_VOLTAGE, VOLTAGE
      LOGICAL :: LIMITED

      INTEGRAL = 2.d0
      PREV_ERROR = 0.d0
      CALL ERROR_DRIVEN_PID_STEP(0.05d0, 0.d0, 3.d5, 0.d0, 0.d0, &
         1.d-8, -6.d4, -1.d0, HUGE(1.d0), .FALSE., &
         -9.d5, -5.d3, .TRUE., INTEGRAL, PREV_ERROR, VOLTAGE, LIMITED)
      CALL ASSERT_CLOSE('positional P first update', VOLTAGE, -7.5d4)
      CALL ASSERT_CLOSE('KI=0 leaves integral unchanged', INTEGRAL, 2.d0)
      IF (LIMITED) ERROR STOP 'unexpected positional P limit'

      PREVIOUS_VOLTAGE = VOLTAGE
      CALL ERROR_DRIVEN_PID_STEP(0.05d0, 0.d0, 3.d5, 0.d0, 0.d0, &
         1.d-8, PREVIOUS_VOLTAGE, -1.d0, HUGE(1.d0), .FALSE., &
         -9.d5, -5.d3, .FALSE., INTEGRAL, PREV_ERROR, VOLTAGE, LIMITED)
      CALL ASSERT_CLOSE('same-sign error keeps correcting', VOLTAGE, -9.0d4)
      CALL ASSERT_CLOSE('repeated KI=0 leaves integral unchanged', INTEGRAL, 2.d0)

      INTEGRAL = 0.d0
      PREV_ERROR = 0.d0
      CALL ERROR_DRIVEN_PID_STEP(0.05d0, 0.d0, 3.d5, 2.d12, 0.d0, &
         1.d-8, -6.d4, -1.d0, 1.d2, .TRUE., &
         -9.d5, -5.d3, .TRUE., INTEGRAL, PREV_ERROR, VOLTAGE, LIMITED)
      CALL ASSERT_CLOSE('100 V cathode slew step one', VOLTAGE, -6.01d4)
      CALL ASSERT_CLOSE('slew anti-windup', INTEGRAL, 0.d0)
      IF (.NOT. LIMITED) ERROR STOP 'slew limit was not reported'

      PREVIOUS_VOLTAGE = VOLTAGE
      CALL ERROR_DRIVEN_PID_STEP(0.05d0, 0.d0, 3.d5, 2.d12, 0.d0, &
         1.d-8, PREVIOUS_VOLTAGE, -1.d0, 1.d2, .TRUE., &
         -9.d5, -5.d3, .FALSE., INTEGRAL, PREV_ERROR, VOLTAGE, LIMITED)
      CALL ASSERT_CLOSE('100 V cathode slew step two', VOLTAGE, -6.02d4)
      CALL ASSERT_CLOSE('repeated slew anti-windup', INTEGRAL, 0.d0)

      INTEGRAL = 0.d0
      PREV_ERROR = 0.d0
      CALL ERROR_DRIVEN_PID_STEP(0.05d0, 0.d0, 1.d6, 1.d0, 0.d0, &
         1.d0, -6.d4, -1.d0, HUGE(1.d0), .TRUE., &
         -9.d4, -5.d3, .TRUE., INTEGRAL, PREV_ERROR, VOLTAGE, LIMITED)
      CALL ASSERT_CLOSE('absolute voltage clamp', VOLTAGE, -9.d4)
      CALL ASSERT_CLOSE('clamp anti-windup', INTEGRAL, 0.d0)
      IF (.NOT. LIMITED) ERROR STOP 'absolute limit was not reported'

      INTEGRAL = 0.d0
      PREV_ERROR = 0.d0
      CALL ERROR_DRIVEN_PID_STEP(0.05d0, 0.d0, 1.d2, 0.d0, 0.d0, &
         1.d0, -60.d0, -1.d0, HUGE(1.d0), .FALSE., &
         -1.d3, 1.d3, .TRUE., INTEGRAL, PREV_ERROR, VOLTAGE, LIMITED)
      CALL ASSERT_CLOSE('negative cathode control direction', VOLTAGE, -65.d0)

      INTEGRAL = 0.d0
      PREV_ERROR = 0.d0
      PREVIOUS_VOLTAGE = -8.d4
      CALL ERROR_DRIVEN_PID_STEP(0.01d0, 0.05d0, 5.d5, 0.d0, 0.d0, &
         1.d-8, PREVIOUS_VOLTAGE, -1.d0, HUGE(1.d0), .FALSE., &
         -1.d5, -5.d3, .TRUE., INTEGRAL, PREV_ERROR, VOLTAGE, LIMITED)
      CALL ASSERT_CLOSE('over-current first correction is less negative', VOLTAGE, -6.d4)

      PREVIOUS_VOLTAGE = VOLTAGE
      CALL ERROR_DRIVEN_PID_STEP(0.01d0, 0.03d0, 5.d5, 0.d0, 0.d0, &
         1.d-8, PREVIOUS_VOLTAGE, -1.d0, HUGE(1.d0), .FALSE., &
         -1.d5, -5.d3, .FALSE., INTEGRAL, PREV_ERROR, VOLTAGE, LIMITED)
      CALL ASSERT_CLOSE('over-current remains less negative as error shrinks', VOLTAGE, -5.d4)

      CALL ERROR_DRIVEN_PID_STEP(0.05d0, 0.05d0, 1.d2, 0.d0, 1.d0, &
         0.d0, -60.d0, -1.d0, HUGE(1.d0), .FALSE., &
         -1.d3, 1.d3, .FALSE., INTEGRAL, PREV_ERROR, VOLTAGE, LIMITED)
      CALL ASSERT_CLOSE('zero dt remains finite', VOLTAGE, -60.d0)

   END SUBROUTINE CHECK_ERROR_DRIVEN_PID

   SUBROUTINE ASSERT_CLOSE(NAME, VALUE, EXPECTED)

      CHARACTER(LEN=*), INTENT(IN) :: NAME
      REAL(KIND=8), INTENT(IN) :: VALUE, EXPECTED
      REAL(KIND=8) :: TOL

      TOL = 1.d-11*MAX(1.d0, ABS(EXPECTED))
      IF (ABS(VALUE - EXPECTED) > TOL) THEN
         WRITE(*,*) 'FAIL: ', TRIM(NAME), VALUE, EXPECTED
         ERROR STOP 1
      END IF

   END SUBROUTINE ASSERT_CLOSE

END PROGRAM test_boundary_control_core
