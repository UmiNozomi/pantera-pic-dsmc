MODULE boundary_control_core

   IMPLICIT NONE

   CONTAINS

   ! Return charge transferred to a material by one surface event in units
   ! of the elementary charge.  Charge state and statistical weight are kept
   ! separate so reflected particles with a changed species/weight are handled
   ! without assuming complete absorption.
   PURE SUBROUTINE SURFACE_EVENT_CHARGE_UNITS(IN_CHARGE, IN_WEIGHT, &
                                               OUT_CHARGE, OUT_WEIGHT, &
                                               PRIMARY_SURVIVES, N_SECONDARY, &
                                               SEE_CHARGE, SEE_WEIGHT, &
                                               Q_ION, Q_ELEC, Q_SEE)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: IN_CHARGE, IN_WEIGHT
      REAL(KIND=8), INTENT(IN) :: OUT_CHARGE, OUT_WEIGHT
      REAL(KIND=8), INTENT(IN) :: SEE_CHARGE, SEE_WEIGHT
      LOGICAL, INTENT(IN) :: PRIMARY_SURVIVES
      INTEGER, INTENT(IN) :: N_SECONDARY
      REAL(KIND=8), INTENT(OUT) :: Q_ION, Q_ELEC, Q_SEE

      Q_ION = 0.d0
      Q_ELEC = 0.d0
      Q_SEE = 0.d0

      IF (IN_CHARGE > 0.d0) THEN
         Q_ION = Q_ION + IN_CHARGE*IN_WEIGHT
      ELSE IF (IN_CHARGE < 0.d0) THEN
         Q_ELEC = Q_ELEC + IN_CHARGE*IN_WEIGHT
      END IF

      IF (PRIMARY_SURVIVES) THEN
         IF (OUT_CHARGE > 0.d0) THEN
            Q_ION = Q_ION - OUT_CHARGE*OUT_WEIGHT
         ELSE IF (OUT_CHARGE < 0.d0) THEN
            Q_ELEC = Q_ELEC - OUT_CHARGE*OUT_WEIGHT
         END IF
      END IF

      IF (N_SECONDARY > 0) THEN
         Q_SEE = -REAL(N_SECONDARY, KIND=8)*SEE_CHARGE*SEE_WEIGHT
      END IF

   END SUBROUTINE SURFACE_EVENT_CHARGE_UNITS

   ! Cubic smoothstep startup command.  The zero slopes at both ends avoid
   ! applying an impulsive electric-field derivative to the initial seed.
   PURE FUNCTION SMOOTH_STARTUP_VOLTAGE(TIME, DURATION, START_VOLTAGE, &
                                        TARGET_VOLTAGE) RESULT(VOLTAGE)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: TIME, DURATION
      REAL(KIND=8), INTENT(IN) :: START_VOLTAGE, TARGET_VOLTAGE
      REAL(KIND=8) :: VOLTAGE, FRACTION

      IF (DURATION <= 0.d0) THEN
         VOLTAGE = TARGET_VOLTAGE
         RETURN
      END IF

      FRACTION = MAX(0.d0, MIN(1.d0, TIME/DURATION))
      FRACTION = FRACTION*FRACTION*(3.d0 - 2.d0*FRACTION)
      VOLTAGE = START_VOLTAGE + &
                (TARGET_VOLTAGE - START_VOLTAGE)*FRACTION

   END FUNCTION SMOOTH_STARTUP_VOLTAGE

   ! Error-driven incremental PID update.  Each correction is referenced to the
   ! last applied voltage, so the command cannot return toward the startup
   ! voltage while the target-current error still has the same sign.
   PURE SUBROUTINE ERROR_DRIVEN_PID_STEP(TARGET_CURRENT, MEASURED_CURRENT, &
                                       KP, KI, KD, CONTROL_DT, &
                                       PREVIOUS_VOLTAGE, &
                                       PID_SIGN, MAX_VOLTAGE_STEP, &
                                       APPLY_VOLTAGE_LIMITS, VOLTAGE_MIN, &
                                       VOLTAGE_MAX, FIRST_UPDATE, &
                                       ERROR_INTEGRAL, ERROR_PREV, &
                                       NEW_VOLTAGE, CONTROL_LIMITED)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: TARGET_CURRENT, MEASURED_CURRENT
      REAL(KIND=8), INTENT(IN) :: KP, KI, KD, CONTROL_DT
      REAL(KIND=8), INTENT(IN) :: PREVIOUS_VOLTAGE
      REAL(KIND=8), INTENT(IN) :: PID_SIGN, MAX_VOLTAGE_STEP
      REAL(KIND=8), INTENT(IN) :: VOLTAGE_MIN, VOLTAGE_MAX
      LOGICAL, INTENT(IN) :: APPLY_VOLTAGE_LIMITS, FIRST_UPDATE
      REAL(KIND=8), INTENT(INOUT) :: ERROR_INTEGRAL, ERROR_PREV
      REAL(KIND=8), INTENT(OUT) :: NEW_VOLTAGE
      LOGICAL, INTENT(OUT) :: CONTROL_LIMITED

      REAL(KIND=8) :: CURRENT_ERROR, ERROR_DERIVATIVE
      REAL(KIND=8) :: INTEGRAL_OLD, INTEGRAL_TRIAL
      REAL(KIND=8) :: VOLTAGE_UNSAT

      CURRENT_ERROR = TARGET_CURRENT - MEASURED_CURRENT
      INTEGRAL_OLD = ERROR_INTEGRAL

      IF (KI > 0.d0 .AND. CONTROL_DT > 0.d0) THEN
         INTEGRAL_TRIAL = INTEGRAL_OLD + CURRENT_ERROR*CONTROL_DT
      ELSE
         INTEGRAL_TRIAL = INTEGRAL_OLD
      END IF

      IF (FIRST_UPDATE .OR. CONTROL_DT <= 0.d0) THEN
         ERROR_DERIVATIVE = 0.d0
      ELSE
         ERROR_DERIVATIVE = (CURRENT_ERROR - ERROR_PREV)/CONTROL_DT
      END IF

      VOLTAGE_UNSAT = PREVIOUS_VOLTAGE + PID_SIGN*( &
         KP*CURRENT_ERROR + KI*INTEGRAL_TRIAL + KD*ERROR_DERIVATIVE)
      CALL APPLY_VOLTAGE_CONSTRAINTS(VOLTAGE_UNSAT, PREVIOUS_VOLTAGE, &
         MAX_VOLTAGE_STEP, APPLY_VOLTAGE_LIMITS, VOLTAGE_MIN, &
         VOLTAGE_MAX, NEW_VOLTAGE)

      CONTROL_LIMITED = ABS(NEW_VOLTAGE - VOLTAGE_UNSAT) > &
         16.d0*EPSILON(1.d0)*MAX(1.d0, ABS(NEW_VOLTAGE), &
                                 ABS(VOLTAGE_UNSAT))

      ! Conditional integration: reject a trial integral only when a voltage
      ! constraint blocks motion in the direction that would reduce the error.
      IF (CONTROL_LIMITED .AND. &
          (VOLTAGE_UNSAT - NEW_VOLTAGE)*(PID_SIGN*CURRENT_ERROR) > 0.d0) THEN
         ERROR_INTEGRAL = INTEGRAL_OLD
         VOLTAGE_UNSAT = PREVIOUS_VOLTAGE + PID_SIGN*( &
            KP*CURRENT_ERROR + KI*ERROR_INTEGRAL + KD*ERROR_DERIVATIVE)
         CALL APPLY_VOLTAGE_CONSTRAINTS(VOLTAGE_UNSAT, PREVIOUS_VOLTAGE, &
            MAX_VOLTAGE_STEP, APPLY_VOLTAGE_LIMITS, VOLTAGE_MIN, &
            VOLTAGE_MAX, NEW_VOLTAGE)
         CONTROL_LIMITED = ABS(NEW_VOLTAGE - VOLTAGE_UNSAT) > &
            16.d0*EPSILON(1.d0)*MAX(1.d0, ABS(NEW_VOLTAGE), &
                                    ABS(VOLTAGE_UNSAT))
      ELSE
         ERROR_INTEGRAL = INTEGRAL_TRIAL
      END IF

      ERROR_PREV = CURRENT_ERROR

   END SUBROUTINE ERROR_DRIVEN_PID_STEP

   PURE SUBROUTINE APPLY_VOLTAGE_CONSTRAINTS(VOLTAGE_COMMAND, &
                                              PREVIOUS_VOLTAGE, &
                                              MAX_VOLTAGE_STEP, &
                                              APPLY_VOLTAGE_LIMITS, &
                                              VOLTAGE_MIN, VOLTAGE_MAX, &
                                              VOLTAGE_APPLIED)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: VOLTAGE_COMMAND, PREVIOUS_VOLTAGE
      REAL(KIND=8), INTENT(IN) :: MAX_VOLTAGE_STEP
      REAL(KIND=8), INTENT(IN) :: VOLTAGE_MIN, VOLTAGE_MAX
      LOGICAL, INTENT(IN) :: APPLY_VOLTAGE_LIMITS
      REAL(KIND=8), INTENT(OUT) :: VOLTAGE_APPLIED

      VOLTAGE_APPLIED = VOLTAGE_COMMAND
      IF (MAX_VOLTAGE_STEP > 0.d0 .AND. &
          MAX_VOLTAGE_STEP < HUGE(1.d0)) THEN
         VOLTAGE_APPLIED = MAX(PREVIOUS_VOLTAGE - MAX_VOLTAGE_STEP, &
            MIN(PREVIOUS_VOLTAGE + MAX_VOLTAGE_STEP, VOLTAGE_APPLIED))
      END IF
      IF (APPLY_VOLTAGE_LIMITS) THEN
         VOLTAGE_APPLIED = MAX(VOLTAGE_MIN, MIN(VOLTAGE_MAX, VOLTAGE_APPLIED))
      END IF

   END SUBROUTINE APPLY_VOLTAGE_CONSTRAINTS

END MODULE boundary_control_core
