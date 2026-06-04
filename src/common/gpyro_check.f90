MODULE GPYRO_CHECK

USE PREC
USE GPYRO_VARS
USE GPYRO_IO

IMPLICIT NONE

CONTAINS


!******************************************************************************	
SUBROUTINE CHECK_GPYRO()
!******************************************************************************	

CALL CHECK_SPROPS
CALL CHECK_GPROPS
CALL CHECK_RXNS
CALL CHECK_HGRXNS
CALL CHECK_IC
CALL CHECK_ALLBC
CALL CHECK_GEOM
CALL CHECK_CASES
CALL CHECK_GENERAL
CALL CHECK_OUTPUT

!******************************************************************************
END SUBROUTINE CHECK_GPYRO
!******************************************************************************

!******************************************************************************  
SUBROUTINE CHECK_SPROPS  
!******************************************************************************  
INTEGER :: ISPEC
CHARACTER(3) :: THREE
CHARACTER(1000) :: MESSAGE

GPG%SOLVE_POROSITY = .FALSE.
GPG%IS_DENSITY_DEPENDENT_ON_TEMPERATURE = .FALSE.

DO ISPEC = 1, SPROP%NSSPEC  
   WRITE(THREE, '(I0)') ISPEC  

   CALL CHECK_DENSITY 
   IF (ALL(GPG%ZEROD(1:GPG%NCASES))) RETURN  ! If all ATG the density is the only needed parameter
 
   CALL CHECK_THERMAL_CONDUCTIVITY
   CALL CHECK_HEAT_CAPACITY
   CALL CHECK_EMISSIVITY
   CALL CHECK_BULK_DENSITY
ENDDO

IF (GPG%SOLVE_PRESSURE) THEN
   DO ISPEC = 1, SPROP%NSSPEC
      WRITE(THREE, '(I0)') ISPEC
      CALL CHECK_PERMEABILITY
   ENDDO
ENDIF

CONTAINS  
!******************************************!
SUBROUTINE CHECK_DENSITY()  

   IF (SPROP%R0(ISPEC) .LT. 0 ) THEN  
      MESSAGE = 'Error: The density of species ISPEC=' // THREE // ' is not defined. ' // &  
                  'Please set R0(' // TRIM(THREE) // ') in &GPYRO_SPORPS. ' // &  
                  '(Note that this quantity is linked to the reaction stoichiometry.)'  
      CALL SHUTDOWN_GPYRO(MESSAGE)  
   ENDIF  

   IF (SPROP%NR(ISPEC) .NE. 0) THEN  
      MESSAGE = 'Warning: The density of species ISPEC=' // THREE // ' depends on temperature. ' // &  
                  'This may cause the stoichiometric coefficients to vary over temperature. ' // &  
                  'If this is not desired, set NR(' // TRIM(THREE) // ')=0 in &GPYRO_SPORPS.'  
      WRITE(0,*) MESSAGE  
      GPG%IS_DENSITY_DEPENDENT_ON_TEMPERATURE = .TRUE.  
   ENDIF  
END SUBROUTINE CHECK_DENSITY 

!******************************************!
SUBROUTINE CHECK_THERMAL_CONDUCTIVITY()  

   IF (SPROP%K0Z(ISPEC) .LT. 0 ) THEN  
      MESSAGE='Error: The thermal conductivity of species ISPEC='// THREE // ' is not defined. ' // &  
               'Please set K0Z(' // TRIM(THREE) // ') in &GPYRO_SPORPS.'  
      CALL SHUTDOWN_GPYRO(MESSAGE)  
   END IF  

   ! Default isotropic values  
   IF (SPROP%K0X(ISPEC) .LT. 0     ) SPROP%K0X(ISPEC)= SPROP%K0Z(ISPEC)
   IF (SPROP%K0Y(ISPEC) .LT. 0     ) SPROP%K0Y(ISPEC)= SPROP%K0Z(ISPEC)
   IF (SPROP%NKX(ISPEC) .LT. -1000 ) SPROP%NKX(ISPEC)= SPROP%NKZ(ISPEC)
   IF (SPROP%NKY(ISPEC) .LT. -1000 ) SPROP%NKY(ISPEC)= SPROP%NKZ(ISPEC)
END SUBROUTINE CHECK_THERMAL_CONDUCTIVITY  

!******************************************!
SUBROUTINE CHECK_HEAT_CAPACITY()  
   IF (SPROP%C0(ISPEC) .LT. 0) THEN  
      MESSAGE='Error: The specific heat capacity of species ISPEC='// THREE // ' is not defined. ' // &  
               'Please set C0(' // TRIM(THREE) // ') in &GPYRO_SPORPS.'  
      CALL SHUTDOWN_GPYRO(MESSAGE)  
   ENDIF  
END SUBROUTINE CHECK_HEAT_CAPACITY  

!******************************************!
SUBROUTINE CHECK_EMISSIVITY()
   IF (SPROP%EMIS(ISPEC) .LT. 0 ) THEN  
      MESSAGE='Error: The emissivity of species ISPEC='// THREE // ' is not defined. ' // &  
               'Please set EMIS(' // TRIM(THREE) // ') in &GPYRO_SPORPS.'  
      CALL SHUTDOWN_GPYRO(MESSAGE)  
   ENDIF  
END SUBROUTINE CHECK_EMISSIVITY 

!******************************************!
SUBROUTINE CHECK_BULK_DENSITY()

   IF (SPROP%RS0(ISPEC) .LT. 0D0) SPROP%RS0(ISPEC) = SPROP%R0(ISPEC)  ! By defaut if RS0 not specified

   IF (SPROP%R0(ISPEC) .GT. SPROP%RS0(ISPEC)) THEN  
      MESSAGE = 'Error: The bulk density R0 of species ISPEC=' // THREE // ' cannot be greater than the pure solid density RS0. ' // &  
                  'Please change RS0(' // TRIM(THREE) // ') or R0(' // TRIM(THREE) // ') in &GPYRO_SPORPS.'  
      CALL SHUTDOWN_GPYRO(MESSAGE)  
   END IF  

   IF (SPROP%R0(ISPEC) .NE. SPROP%RS0(ISPEC)) GPG%SOLVE_POROSITY = .TRUE.  !solve porosity only if is needed

END SUBROUTINE CHECK_BULK_DENSITY  

!******************************************!
SUBROUTINE CHECK_PERMEABILITY()

   IF (SPROP%PERMZ(ISPEC) .LT. 0D0) THEN  
      MESSAGE = 'Error: The permeabilty PERMZ of species ISPEC=' // THREE // 'sould be defined when SOLVE_PRESSURE=.TRUE.. '
      CALL SHUTDOWN_GPYRO(MESSAGE)
   END IF

   ! Default isotropic values  
   IF (SPROP%PERMX(ISPEC) .LT. 0 ) SPROP%PERMX(ISPEC)= SPROP%PERMZ(ISPEC)
   IF (SPROP%PERMY(ISPEC) .LT. 0 ) SPROP%PERMY(ISPEC)= SPROP%PERMZ(ISPEC)

END SUBROUTINE CHECK_PERMEABILITY



END SUBROUTINE CHECK_SPROPS  
!******************************************************************************  

!******************************************************************************
SUBROUTINE CHECK_RXNS
!******************************************************************************
INTEGER  J,IRXN
CHARACTER(2) :: TWO
CHARACTER(1000) :: MESSAGE
REAL(EB) :: SUMVAL

DO IRXN = 1, SPROP%NRXN
   SUMVAL = 0D0
   DO J = 1, GPROP%NGSPEC
      SUMVAL = SUMVAL + GPROP%YIELDS(J,IRXN)
   ENDDO
   IF (SUMVAL .LT. 0.999999 .OR. SUMVAL .GT. 1.000001) THEN
      WRITE(TWO,'(I2.2)') IRXN
      MESSAGE='Be careful. Gaseous yields for heterogeneous reaction ' // TWO // ' do not sum to unity.'
      CALL SHUTDOWN_GPYRO(MESSAGE) ! This can be commented out to circumvent this
      IF (RXN(IRXN)%CHI .LT. SPROP%R0(1) / (SPROP%R0(1) - SPROP%R0(2) ) ) THEN
         WRITE(TWO,'(I2.2)') IRXN
         MESSAGE='Error: CHI value for heterogeneous reaction ' // TWO // ' is lower than possible for a condensation reaction.'
         CALL SHUTDOWN_GPYRO(MESSAGE) !This can be commented out to circumvent this 
      ENDIF
   ENDIF
ENDDO

!******************************************************************************
END SUBROUTINE CHECK_RXNS
!******************************************************************************

!******************************************************************************
SUBROUTINE CHECK_HGRXNS
!******************************************************************************
INTEGER  J,IRXN
CHARACTER(2) :: TWO
CHARACTER(1000) :: MESSAGE
REAL(EB) :: SUMVAL


DO IRXN = 1, GPROP%NHGRXN
   SUMVAL = 0D0
   DO J = 1, GPROP%NGSPEC
      SUMVAL = SUMVAL + GPROP%HGYIELDS(J,IRXN)
   ENDDO
   IF (SUMVAL .LT. -0.000001 .OR. SUMVAL .GT. 0.000001) THEN
      WRITE(TWO,'(I2.2)') IRXN
      MESSAGE='Problem with gaseous yields for homogeneous gaseous reaction ' // TWO // '. Yields do not sum to 0'
      CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF
ENDDO
!******************************************************************************
END SUBROUTINE CHECK_HGRXNS
!******************************************************************************

!******************************************************************************  
SUBROUTINE CHECK_GEOM  
!******************************************************************************  
INTEGER :: IOBST
INTEGER :: IMESH, NCELLZ, NCELLX, NCELLY, NMESH
CHARACTER(1000) :: MESSAGE
CHARACTER(LEN=12) :: FACENAME(6)
INTEGER :: IFACE, ISURF, BC_IDX, ICNUM
LOGICAL :: FOUND
CHARACTER(LEN=20) ::  STR_IMESH, STR_BCIDX, STR_FACEIDX

! Face order and names
FACENAME(1) = 'west (-x)'  
FACENAME(2) = 'east (+x)'  
FACENAME(3) = 'south (-y)'  
FACENAME(4) = 'north (+y)'  
FACENAME(5) = 'top (+z)'  
FACENAME(6) = 'bottom (-z)'  

! First, check dimension consistency (0D/1D/2D/3D)
CALL CHECK_MESH_DIMENSION_VALIDITY
! Check all mesh BC and IC
CALL CHECK_MESH_BC_VALIDITY  
CALL CHECK_MESH_IC_VALIDITY 
! Check all OBST BC and IC
CALL CHECK_OBST_BC_VALIDITY
CALL CHECK_OBST_IC_VALIDITY

CONTAINS
!******************************************!
SUBROUTINE CHECK_MESH_DIMENSION_VALIDITY()  
   NMESH = GPG%NUM_GPYRO_MESHES  
   DO IMESH = 1, NMESH  

      NCELLZ = GPG%NCELLZ(IMESH)  
      NCELLX = GPG%NCELLX(IMESH)  
      NCELLY = GPG%NCELLY(IMESH)  

      ! Check 1D validity
      IF (NCELLZ .EQ. 1 .AND. NCELLX .GT. 1 .AND. NCELLY .EQ. 1) THEN  
         MESSAGE = 'Error: For 1D simulations set NCELLY = 1, NCELLX = 1, and NCELLZ > 1'  
         CALL SHUTDOWN_GPYRO(MESSAGE)  
      ENDIF  

      IF (NCELLZ .EQ. 1 .AND. NCELLY .GT. 1 .AND. NCELLX .EQ. 1) THEN  
         MESSAGE = 'Error: For 1D simulations set NCELLY = 1, NCELLX = 1, and NCELLZ > 1'  
         CALL SHUTDOWN_GPYRO(MESSAGE)  
      ENDIF  

      ! Check 2D validity
      IF (NCELLY .GT. 1 .AND. NCELLX .EQ. 1) THEN  
         MESSAGE = 'Error: For 2D simulations set NCELLY = 1 and NCELLX > 1'  
         CALL SHUTDOWN_GPYRO(MESSAGE)  
      ENDIF  
   ENDDO
END SUBROUTINE CHECK_MESH_DIMENSION_VALIDITY  

!******************************************!
SUBROUTINE CHECK_MESH_BC_VALIDITY()
   NMESH = GPG%NUM_GPYRO_MESHES  
   DO IMESH = 1, NMESH  
      DO IFACE = 1, 6 
         ! Skip unused directions depending on mesh dimension
         IF ((GPG%NCELLX(IMESH) .LE. 1) .AND. (IFACE == 1 .OR. IFACE == 2)) CYCLE  
         IF ((GPG%NCELLY(IMESH) .LE. 1) .AND. (IFACE == 3 .OR. IFACE == 4)) CYCLE  
         IF ((GPG%NCELLZ(IMESH) .LE. 1) .AND. (IFACE == 5 .OR. IFACE == 6)) CYCLE  

         BC_IDX = GPG%DEFAULT_SURF_IDX(IMESH, IFACE)

         ! Convert integers to strings
         WRITE(STR_BCIDX,'(I0)') BC_IDX
         WRITE(STR_IMESH,'(I0)') IMESH 
         WRITE(STR_FACEIDX,'(I0)') IFACE

         ! Check if BC index is defined (positive)
         IF (BC_IDX .LE. 0) THEN  
               MESSAGE = 'Error: Boundary condition not defined for ' // TRIM(FACENAME(IFACE)) // &  
                        ' face of mesh ' // TRIM(STR_IMESH) // &  
                        '. Please set DEFAULT_SURF_IDX('// TRIM(STR_IMESH)//','// TRIM(STR_FACEIDX)// ').'  
               CALL SHUTDOWN_GPYRO(MESSAGE)  
         ENDIF  

         ! Check if BC index exists in SURF_IDX list
         FOUND = .FALSE.  
         DO ISURF = 1, GPG%NSURF_IDX  
               IF (BC_IDX == GPG%ALLBC(ISURF)%SURF_IDX) THEN  
                  FOUND = .TRUE.  
                  EXIT  
               ENDIF  
         END DO  

         IF (.NOT. FOUND) THEN  
               MESSAGE = 'Error: Invalid boundary condition for ' // TRIM(FACENAME(IFACE)) // &  
                        ' face of mesh ' // TRIM(STR_IMESH) // &  
                        '. The boundary index ' // TRIM(STR_BCIDX) // ' is not defined.' // &  
                        ' Please change DEFAULT_SURF_IDX('// TRIM(STR_IMESH)//','// TRIM(STR_FACEIDX)// ').'  
               CALL SHUTDOWN_GPYRO(MESSAGE)  
         ENDIF  
      ENDDO  
   ENDDO

   DO IOBST = 1, GPG%NOBST
      IMESH = GPG%GEOM(IOBST)%IMESH 
      DO IFACE = 1, 6 
         ! Skip unused directions depending on mesh dimension
         IF ((GPG%NCELLX(IMESH) .LE. 1) .AND. (IFACE == 1 .OR. IFACE == 2)) CYCLE  
         IF ((GPG%NCELLY(IMESH) .LE. 1) .AND. (IFACE == 3 .OR. IFACE == 4)) CYCLE  
         IF ((GPG%NCELLZ(IMESH) .LE. 1) .AND. (IFACE == 5 .OR. IFACE == 6)) CYCLE  
         
         BC_IDX = GPG%GEOM(IOBST)%SURF_IDX(IFACE)

         !If not provided take default mesh value
         IF (BC_IDX .LE. 0 ) THEN
            GPG%GEOM(IOBST)%SURF_IDX(IFACE) =  GPG%DEFAULT_SURF_IDX(IMESH, IFACE)
         ENDIF
      ENDDO
   ENDDO

END SUBROUTINE CHECK_MESH_BC_VALIDITY  

!******************************************!
SUBROUTINE CHECK_MESH_IC_VALIDITY()
   NMESH = GPG%NUM_GPYRO_MESHES  

   DO IMESH = 1, NMESH  
      ICNUM = GPG%DEFAULT_IC(IMESH)  
      ! Convert integers to strings
      WRITE(STR_BCIDX,'(I0)') ICNUM  
      WRITE(STR_IMESH,'(I0)') IMESH  

      ! Check if IC index is within valid range
      IF (ICNUM .GT. GPG%NIC) THEN  
         MESSAGE = 'Error: Invalid initial condition for mesh ' // TRIM(STR_IMESH) // &  
                     '. The initial condition index ' // TRIM(STR_BCIDX) // ' is not defined.' // &  
                     ' Please change DEFAULT_IC('// TRIM(STR_IMESH)//').'  
         CALL SHUTDOWN_GPYRO(MESSAGE)  
      ENDIF
   ENDDO
END SUBROUTINE CHECK_MESH_IC_VALIDITY  

!******************************************!
SUBROUTINE CHECK_OBST_BC_VALIDITY()

   DO IOBST = 1, GPG%NOBST 
      IMESH = GPG%GEOM(IOBST)%IMESH
      DO IFACE = 1, 6
         ! Skip unused directions depending on mesh dimension
         IF ((GPG%NCELLX(IMESH) .LE. 1) .AND. (IFACE == 1 .OR. IFACE == 2)) CYCLE  
         IF ((GPG%NCELLY(IMESH) .LE. 1) .AND. (IFACE == 3 .OR. IFACE == 4)) CYCLE  
         IF ((GPG%NCELLZ(IMESH) .LE. 1) .AND. (IFACE == 5 .OR. IFACE == 6)) CYCLE  

         BC_IDX = GPG%GEOM(IOBST)%SURF_IDX(IFACE)  
         IF (BC_IDX .LE. 0) CYCLE  

         ! Convert integers to strings
         WRITE(STR_BCIDX,'(I0)') BC_IDX  
         WRITE(STR_IMESH,'(I0)') IOBST  
         WRITE(STR_FACEIDX,'(I0)') IFACE  

         ! Check if BC index exists in SURF_IDX list
         FOUND = .FALSE.  
         DO ISURF = 1, GPG%NSURF_IDX  
               IF (BC_IDX == GPG%ALLBC(ISURF)%SURF_IDX) THEN  
                  FOUND = .TRUE.  
                  EXIT  
               ENDIF  
         END DO  

         IF (.NOT. FOUND) THEN  
               MESSAGE = 'Error: Invalid boundary condition for ' // TRIM(FACENAME(IFACE)) // &  
                        ' face of OBST ' // TRIM(STR_IMESH) // &  
                        '. The boundary index ' // TRIM(STR_BCIDX) // ' is not defined.' // &  
                        ' Please change SURF_IDX2D('// TRIM(STR_IMESH)//','// TRIM(STR_FACEIDX)// ').'  
               CALL SHUTDOWN_GPYRO(MESSAGE)  
         ENDIF  
      ENDDO  
   ENDDO
END SUBROUTINE CHECK_OBST_BC_VALIDITY  

!******************************************!
SUBROUTINE CHECK_OBST_IC_VALIDITY()
   DO IOBST = 1, GPG%NOBST
      ICNUM = GPG%GEOM(IOBST)%ICNUM  
      ! Convert integers to strings
      WRITE(STR_BCIDX,'(I0)') ICNUM  
      WRITE(STR_IMESH,'(I0)') IOBST  
      ! Check if IC index is within valid range
      IF (ICNUM .GT. GPG%NIC) THEN  
         MESSAGE = 'Error: Invalid initial condition for OBST ' // TRIM(STR_IMESH) // &  
                     '. The initial condition index ' // TRIM(STR_BCIDX) // ' is not defined.' // &  
                     ' Please change ICNUM('// TRIM(STR_IMESH)//').'  
         CALL SHUTDOWN_GPYRO(MESSAGE)  
      ENDIF  
   ENDDO
END SUBROUTINE CHECK_OBST_IC_VALIDITY  

END SUBROUTINE CHECK_GEOM  
!******************************************************************************  


!******************************************************************************
SUBROUTINE CHECK_GPROPS
!******************************************************************************
INTEGER ::IGSPEC
CHARACTER(1000) :: MESSAGE
CHARACTER(3) :: THREE
CHARACTER(20) :: STATUS


IF (GPG%SOLVE_PRESSURE .OR. GPG%SOLVE_GAS_YJ .OR. GPG%SOLVE_GAS_ENERGY) THEN
   CALL CHECK_BACKGROUND_SPECIES
   CALL CHECK_OXIDATIVE_SPECIES

   IGSPEC = GPROP%IBG
   STATUS = "background"
   WRITE(THREE, '(I0)') IGSPEC
   CALL CHECK_SIGMA
   CALL CHECK_EPSOK

   IGSPEC = GPROP%IO2
   STATUS = "oxidative"

   CALL CHECK_SIGMA
   CALL CHECK_EPSOK

   DO IGSPEC = 1, SPROP%NSSPEC
      WRITE(THREE, '(I0)') IGSPEC
      CALL CHECK_MOLAR_MASS
   ENDDO
ENDIF

CONTAINS

!******************************************!
SUBROUTINE CHECK_SIGMA
   IF (GPROP%SIGMA(IGSPEC) .LT. 0) THEN
       MESSAGE = 'Error: The σ parameter of Lennard-Jones model for the ' // TRIM(STATUS) // &
                 ' species IGSPEC=' // TRIM(THREE) // ' is not defined, and it is needed for the gas phase solvers. ' // &
                 'Please set SIGMA(' // TRIM(THREE) // ') in &GPYRO_GPROPS.'
       CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF
END SUBROUTINE CHECK_SIGMA

!******************************************!
SUBROUTINE CHECK_EPSOK()  
   IF (GPROP%EPSOK(IGSPEC) .LT. 0) THEN  
      MESSAGE = 'Error: The ε/k parameter of Lennard-Jones model for the ' // TRIM(STATUS) // &
                  ' species IGSPEC=' // TRIM(THREE) // ' is not defined, and it is needed for the gas phase solvers. ' // &
                  'Please set EPSOK(' // TRIM(THREE) // ') in &GPYRO_GPROPS.'
      CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF  
END SUBROUTINE CHECK_EPSOK

!******************************************!
SUBROUTINE CHECK_MOLAR_MASS()  
   IF (GPROP%M(IGSPEC) .LT. 0) THEN  
      MESSAGE='Error: The molar mass of gas species IGSPEC='// THREE // ' is not defined, and it is needed for the gas phase solvers. ' // &  
               'Please set M(' // TRIM(THREE) // ') in &GPYRO_GPORPS.'
      CALL SHUTDOWN_GPYRO(MESSAGE)  
   ENDIF  
END SUBROUTINE CHECK_MOLAR_MASS

!******************************************!
SUBROUTINE CHECK_BACKGROUND_SPECIES()
   IF (GPROP%IBG .GT. GPROP%NGSPEC) THEN
      MESSAGE = 'Error: The index of the background gas species is larger than the number of gas species NGSPEC. ' // &
                'Please change IBG in &GPYRO_GPROPS.'
      CALL SHUTDOWN_GPYRO(MESSAGE)

  ELSEIF (GPROP%IBG .LT. 0) THEN
      MESSAGE = 'Error: The background species is not defined, and it is needed for the gas phase solvers. ' // &
                'Please set IBG in &GPYRO_GPROPS.'
      CALL SHUTDOWN_GPYRO(MESSAGE)
  ENDIF

END SUBROUTINE CHECK_BACKGROUND_SPECIES

!******************************************!

SUBROUTINE CHECK_OXIDATIVE_SPECIES()
   IF (GPROP%IO2 .GT. GPROP%NGSPEC) THEN
      MESSAGE = 'Error: The index of the oxidative gas species is larger than the number of gas species NGSPEC. ' // &
                'Please change IO2 in &GPYRO_GPROPS.'
      CALL SHUTDOWN_GPYRO(MESSAGE)

  ELSEIF (GPROP%IO2 .LT. 0) THEN
      MESSAGE = 'Error: The oxidative species is not defined, and it is needed for the gas phase solvers. ' // &
                'Please set IO2 in &GPYRO_GPROPS.'
      CALL SHUTDOWN_GPYRO(MESSAGE)
  ENDIF

END SUBROUTINE CHECK_OXIDATIVE_SPECIES

!******************************************************************************
END SUBROUTINE CHECK_GPROPS
!******************************************************************************

!******************************************************************************
SUBROUTINE CHECK_GENERAL
!******************************************************************************
INTEGER :: ISPEC
CHARACTER(1000) :: MESSAGE

IF (GPROP%NHGRXN .GT. 0) THEN
   IF ( (.NOT. GPG%SOLVE_GAS_ENERGY) ) THEN
      MESSAGE='Set SOLVE_GAS_ENERGY = .TRUE. when using homogeneous gaseous reactions (NHGRXN > 0).'
      CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF
ENDIF

!SOLVE_POROSITY is determined in CHECK_SPROPS
IF (GPG%SOLVE_GAS_ENERGY .AND. (.NOT. GPG%SOLVE_POROSITY)) THEN
   MESSAGE = 'Error: SOLVE_GAS_ENERGY requires porosity to be resolved.' // NEW_LINE('A') // &
             'Porosity can only be computed if the pure solid density RS0 of at least one ' // &
             'condensed species is defined and different from its bulk density R0.'
   CALL SHUTDOWN_GPYRO(MESSAGE)
ENDIF

IF (GPG%SOLVE_GAS_YJ .AND. (.NOT. GPG%SOLVE_POROSITY)) THEN
   MESSAGE = 'Error: SOLVE_GAS_YJ requires porosity to be resolved.' // NEW_LINE('A') // &
             'Porosity can only be computed if the pure solid density RS0 of at least one ' // &
             'condensed species is defined and different from its bulk density R0.'
   CALL SHUTDOWN_GPYRO(MESSAGE)
ENDIF

IF (GPG%SOLVE_PRESSURE .AND. (.NOT. GPG%SOLVE_POROSITY)) THEN
   MESSAGE = 'Error: SOLVE_PRESSURE requires porosity to be resolved.' // NEW_LINE('A') // &
             'Porosity can only be computed if the pure solid density RS0 of at least one ' // &
             'condensed species is defined and different from its bulk density R0.'
   CALL SHUTDOWN_GPYRO(MESSAGE)
ENDIF

IF (GPG%THERMAL_EQUILIBRIUM .AND. GPG%FULL_QSG .AND. (.NOT. GPG%SOLVE_POROSITY)) THEN
   MESSAGE = 'Error: THERMAL_EQUILIBRIUM with FULL_QSG requires porosity to be resolved.' // NEW_LINE('A') // &
             'Porosity can only be computed if the pure solid density RS0 of at least one ' // &
             'condensed species is defined and different from its bulk density R0.'
   CALL SHUTDOWN_GPYRO(MESSAGE)
ENDIF



IF (GPG%FDSMODE) THEN
   IF (GPG%FDS_MATL_VER .NE. 5 .AND. GPG%FDS_MATL_VER .NE. 6) THEN
      MESSAGE='Error, when FDSMODE=.TRUE., set FDS_MATL_VER to either 5 or 6'
      CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF   
ENDIF

! Warnings
IF (.NOT. GPG%SHYI_CORRECTION) THEN
   DO ISPEC = 1, SPROP%NSSPEC
      IF (SPROP%C0(ISPEC) .NE. SPROP%C0(1) .OR. SPROP%NC(ISPEC) .NE. SPROP%NC(1) ) THEN
         WRITE(*,*) '***** WARNING ***** '
         WRITE(*,*) 'If different solid species have different specific heat capacities '
         WRITE(*,*) 'Then Gpyro should be run with GPG%SHYI_CORRECTION = .TRUE. '
         WRITE(*,*)
      ENDIF   
   ENDDO
ENDIF

IF (GPG%EXPLICIT_T) THEN
   WRITE(*,*) '***** WARNING ***** '
   WRITE(*,*) 'Gpyro is an implicit code.' 
   WRITE(*,*) 'Setting EXPLICIT_T = .TRUE. should be used with caution.'
   WRITE(*,*)
ENDIF

IF (GPG%FDSMODE) THEN
   WRITE(*,*) '***** WARNING ***** '
   WRITE(*,*) 'Gpyro is running in FDS mode.' 
   WRITE(*,*) "Be sure to check Gpyro's .out file as some defaults are different in FDS mode."
   WRITE(*,*)
ENDIF

!******************************************************************************
END SUBROUTINE CHECK_GENERAL
!******************************************************************************

!******************************************************************************  
SUBROUTINE CHECK_IC  
!******************************************************************************  
INTEGER :: I, J, IRXN, ICINDEX  
CHARACTER(2) :: TWO  
CHARACTER(3) :: THREE  
CHARACTER(1000) :: MESSAGE  
REAL(EB) :: SUMVAL  

CALL CHECK_NEED_GAS_MASS_FRACTION  ! Determine if gas-phase mass fractions is need

DO ICINDEX = 1, GPG%NIC  
   CALL CHECK_IC_CONDENSED_PHASE
   IF (GPG%NEED_GAS_YJ) CALL CHECK_IC_GAS_PHASE
END DO

CONTAINS

!******************************************!
SUBROUTINE CHECK_NEED_GAS_MASS_FRACTION  

   GPG%NEED_GAS_YJ = GPG%SOLVE_GAS_YJ  !If the user specified solving gas YJ, then it is of course needed
   ! It can be needed even if SOLVE_GAS_YJ = False; in this case, the gas mass fraction will be the initial mass fraction
   IF (GPG%SOLVE_PRESSURE) GPG%NEED_GAS_YJ = .TRUE. !Needed to get gas density
   DO IRXN = 1, SPROP%NRXN
      DO J = 1, GPROP%NGSPEC
         ! That means a gaseous species is consumed, so the gas mass fraction is needed  
         IF (GPROP%YIELDS(J, IRXN) .LT. 0D0) GPG%NEED_GAS_YJ = .TRUE.
      ENDDO  
   ENDDO  
END SUBROUTINE CHECK_NEED_GAS_MASS_FRACTION  

!******************************************!
SUBROUTINE CHECK_IC_CONDENSED_PHASE()

   SUMVAL = 0D0  
   DO I = 1, SPROP%NSSPEC  
      IF (GPG%INITIAL_CONDITIONS(ICINDEX)%YI0(I) .GT. 1 .OR. &  
            GPG%INITIAL_CONDITIONS(ICINDEX)%YI0(I) .LT. 0) THEN  
            WRITE(THREE, '(I0)') ICINDEX  
            WRITE(TWO , '(I0)') I  
            MESSAGE = "Error: The initial condensed-phase mass fractions (kg/kg) must be between 0 and 1. " // &  
                     "Please check YI0(" // TRIM(THREE) // "," // TRIM(TWO) // ") in &GPYRO_IC."  
            CALL SHUTDOWN_GPYRO(MESSAGE)  
      ENDIF  
      SUMVAL = SUMVAL + GPG%INITIAL_CONDITIONS(ICINDEX)%YI0(I)  
   ENDDO  

   IF (SUMVAL .LT. 0.999999 .OR. SUMVAL .GT. 1.000001) THEN  
      WRITE(THREE,'(I0)') ICINDEX  
      MESSAGE='Error: Initial condensed-phase mass fractions for the initial condition ' // THREE // &  
               ' do not sum to 1.0. Please Check YI0(' // THREE // ', 1:NSPEC) in &GPYRO_IC.'  
      CALL SHUTDOWN_GPYRO(MESSAGE)  
   ENDIF  
END SUBROUTINE CHECK_IC_CONDENSED_PHASE  

!******************************************!
SUBROUTINE CHECK_IC_GAS_PHASE()
SUMVAL = 0D0  
DO J = 1, GPROP%NGSPEC
   SUMVAL = SUMVAL + GPG%INITIAL_CONDITIONS(ICINDEX)%YJ0(J)
ENDDO

IF (SUMVAL .LT. 0.999999 .OR. SUMVAL .GT. 1.000001) THEN
   WRITE(THREE,'(I3.3)') ICINDEX
   MESSAGE='Error: Initial gas-phase mass fractions do not sum to 1.0 for ICNUM ' // THREE // &
            '. Check IC worksheet and make sure you have specified initial conditions for all gas-phase species.'
   CALL SHUTDOWN_GPYRO(MESSAGE)  
ENDIF  
END SUBROUTINE CHECK_IC_GAS_PHASE

!******************************************************************************
END SUBROUTINE CHECK_IC
!******************************************************************************


!******************************************************************************
SUBROUTINE CHECK_ALLBC  
!******************************************************************************
INTEGER :: J, IOBST
INTEGER :: IMESH, IFACE, ISURF
CHARACTER(3) :: THREE
CHARACTER(1000) :: MESSAGE
REAL(EB) :: SUMVAL
CHARACTER(LEN=12) :: FACENAME(6)
CHARACTER(LEN=20) :: STR_IMESH, STR_BCIDX, STR_FACEIDX

! Face order and names
FACENAME(1) = 'west (-x)'
FACENAME(2) = 'east (+x)'
FACENAME(3) = 'south (-y)'
FACENAME(4) = 'north (+y)'
FACENAME(5) = 'top (+z)'
FACENAME(6) = 'bottom (-z)'
! Check gas transport boundary conditions for meshes and obstacles if instantaneous transport is used
IF (.NOT. GPG%SOLVE_PRESSURE) THEN
   CALL CHECK_MESH_BC_INSTANTANEOUS_GAS_RELEASE
   CALL CHECK_OBST_BC_INSTANTANEOUS_GAS_RELEASE
ENDIF

! Check boundary conditions for gas species mass fractions
IF (GPG%SOLVE_GAS_YJ) CALL CHECK_GAS_SPECIES_BC

CONTAINS 

!******************************************!
SUBROUTINE CHECK_MESH_BC_INSTANTANEOUS_GAS_RELEASE
   INTEGER :: IBC
   DO IMESH = 1, GPG%NUM_GPYRO_MESHES
      DO IFACE = 1, 6  
         IF (IFACE .EQ. 6) CYCLE ! Bottom face is allowed for gas mass flux  
         ! Skip unused directions depending on mesh dimension  
         IF ( (GPG%NCELLX(IMESH) .LE. 1) .AND. (IFACE == 1 .OR. IFACE == 2) ) CYCLE  
         IF ( (GPG%NCELLY(IMESH) .LE. 1) .AND. (IFACE == 3 .OR. IFACE == 4) ) CYCLE  
         IF ( (GPG%NCELLZ(IMESH) .LE. 1) .AND. (IFACE == 5 .OR. IFACE == 6) ) CYCLE  

         IBC = GPG%DEFAULT_SURF_IDX(IMESH, IFACE)
         DO ISURF = 1, GPG%NSURF_IDX
         IF  (GPG%ALLBC(ISURF)%SURF_IDX .NE. IBC) CYCLE
         IF (ABS(GPG%ALLBC(ISURF)%MDOTPP) .GT. EPSILON_EB) THEN  
               WRITE(STR_IMESH,'(I0)') IMESH  
               WRITE(STR_BCIDX,'(I0)') ISURF  
               WRITE(STR_FACEIDX,'(I0)') IFACE  

               MESSAGE = 'Error: Invalid gas transport boundary condition for ' // TRIM(FACENAME(IFACE)) // &
                        ' face of mesh ' // TRIM(STR_IMESH) // '.' // NEW_LINE('A') // &
                        'For instantaneous gas transport (SOLVE_PRESSURE = .FALSE.),' // &
                        ' a gas mass flux can be imposed only on the bottom surface (index 6).' // NEW_LINE('A') // &
                        'However, a mass flux is defined for DEFAULT_SURF_IDX(' // TRIM(STR_IMESH) // ',' // TRIM(STR_FACEIDX) // ')' // &
                        ' related to boundary condition number ' // TRIM(STR_BCIDX) // '.' // NEW_LINE('A') // &
                        'Please set MDOTPP(' // TRIM(STR_BCIDX) // ')=0D0' // &
                        ' or modify the boundary condition associated with the ' // TRIM(FACENAME(IFACE)) // &
                        ' surface of mesh ' // TRIM(STR_IMESH) // ' via DEFAULT_SURF_IDX(' // TRIM(STR_IMESH) // ',' // TRIM(STR_FACEIDX) // ').'
               CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
         ENDDO
      ENDDO
   ENDDO
END SUBROUTINE CHECK_MESH_BC_INSTANTANEOUS_GAS_RELEASE

!******************************************!
SUBROUTINE CHECK_OBST_BC_INSTANTANEOUS_GAS_RELEASE
   INTEGER :: IBC
   DO IOBST = 1, GPG%NOBST 
      IMESH = GPG%GEOM(IOBST)%IMESH  
      DO IFACE = 1, 6  
         IF (IFACE .EQ. 6) CYCLE ! Bottom face is allowed for gas mass flux  
         ! Skip unused directions depending on mesh dimension  
         IF ( (GPG%NCELLX(IMESH) .LE. 1) .AND. (IFACE == 1 .OR. IFACE == 2) ) CYCLE
         IF ( (GPG%NCELLY(IMESH) .LE. 1) .AND. (IFACE == 3 .OR. IFACE == 4) ) CYCLE
         IF ( (GPG%NCELLZ(IMESH) .LE. 1) .AND. (IFACE == 5 .OR. IFACE == 6) ) CYCLE
         IBC = GPG%GEOM(IOBST)%SURF_IDX(IFACE)
         DO ISURF = 1, GPG%NSURF_IDX
         IF  (GPG%ALLBC(ISURF)%SURF_IDX .NE. IBC) CYCLE
         IF (ABS(GPG%ALLBC(ISURF)%MDOTPP) .GT. EPSILON_EB) THEN
               WRITE(STR_IMESH,'(I0)') IOBST
               WRITE(STR_BCIDX,'(I0)') ISURF
               WRITE(STR_FACEIDX,'(I0)') IFACE

               MESSAGE = 'Error: Invalid gas transport boundary condition for ' // TRIM(FACENAME(IFACE)) // &
                        ' face of obstruction ' // TRIM(STR_IMESH) // '.' // NEW_LINE('A') // &
                        'For instantaneous gas transport (SOLVE_PRESSURE = .FALSE.),' // &
                        ' a gas mass flux can be imposed only on the bottom surface (index 6).' // NEW_LINE('A') // &
                        'However, a mass flux is defined for SURF_IDX2D(' // TRIM(STR_IMESH) // ',' // TRIM(STR_FACEIDX) // ')' // &
                        ' related to boundary condition number ' // TRIM(STR_BCIDX) // '.' // NEW_LINE('A') // &
                        'Please set MDOTPP(' // TRIM(STR_BCIDX) // ')=0D0' // &
                        ' or modify the boundary condition associated with the ' // TRIM(FACENAME(IFACE)) // &
                        ' surface of obstruction ' // TRIM(STR_IMESH) // ' via SURF_IDX2D(' // TRIM(STR_IMESH) // ',' // TRIM(STR_FACEIDX) // ').'  
               CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
         ENDDO
      ENDDO
   ENDDO
END SUBROUTINE CHECK_OBST_BC_INSTANTANEOUS_GAS_RELEASE

!******************************************!
SUBROUTINE CHECK_GAS_SPECIES_BC
   DO ISURF = 1, GPG%NSURF_IDX
      WRITE(THREE,'(I0)') ISURF
      SUMVAL = 0D0
      DO J = 1, GPROP%NGSPEC
         SUMVAL = SUMVAL + GPG%ALLBC(ISURF)%YJINF(J)
      ENDDO

      IF (SUMVAL .LT. 0.999999 .OR. SUMVAL .GT. 1.000001) THEN
         MESSAGE = 'Error: Problem with gaseous species boundary condition # ' // THREE // &  
                     '. Mass fractions do not sum to 1.0' 
         CALL SHUTDOWN_GPYRO(MESSAGE)
      ENDIF  

      IF (ABS(GPG%ALLBC(ISURF)%MDOTPP) .GT. EPSILON_FB .AND. GPG%ALLBC(ISURF)%HM .GT. EPSILON_FB) THEN  
         WRITE(THREE,'(I0)') ISURF
         MESSAGE = 'Error: Problem with gaseous species boundary condition # ' // THREE // &  
                     '. Only one of HM or MDOTPP can be nonzero.'
         CALL SHUTDOWN_GPYRO(MESSAGE)  
      ENDIF
   ENDDO
END SUBROUTINE CHECK_GAS_SPECIES_BC
!******************************************************************************
END SUBROUTINE CHECK_ALLBC
!******************************************************************************

!******************************************************************************
SUBROUTINE CHECK_CASES
!******************************************************************************
INTEGER :: ICASE,IMESH, NCELLZ, NCELLX, NCELLY, NMESH
CHARACTER(2) :: TWO
CHARACTER(3) :: THREE
CHARACTER(1000) :: MESSAGE

DO ICASE = 1, GPG%NCASES
   IMESH  = GPG%IMESH (ICASE)
   NCELLZ = GPG%NCELLZ(IMESH)
   NCELLX = GPG%NCELLX(IMESH)
   NCELLY = GPG%NCELLY(IMESH)

   IF (IMESH .GT. GPG%NUM_GPYRO_MESHES) THEN
      WRITE(THREE, '(I0)') ICASE
      WRITE(TWO, '(I0)') NMESH
      MESSAGE = 'Error: For case ' // TRIM(THREE) // ', the mesh index ' // TRIM(TWO) // ' does not exist. ' // &
      'Please change IMESH(' // TRIM(THREE) // ') in &GPYRO_CASES.'
   CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF

   IF (GPG%TSTOP(ICASE) .LT. 0) THEN
      WRITE(THREE, '(I0)') ICASE
       MESSAGE = 'Error: For case ' // TRIM(THREE) // ', the total time of the simulation is not defined '&
      'Please change TSTOP(' // TRIM(THREE) // ') in &GPYRO_CASES.'
      CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF

   IF(GPG%ZEROD(ICASE) .AND. (NCELLX .GT. 1 .OR. NCELLY .GT. 1 .OR. NCELLZ .GT. 1)) THEN
      MESSAGE='Error:  For meshes corresponding to cases where ZEROD=.TRUE., ensure that NCELLZ=NCELLX=NCELLY=1.'
      CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF
ENDDO

!******************************************************************************
END SUBROUTINE CHECK_CASES
!******************************************************************************

!******************************************************************************
SUBROUTINE CHECK_OUTPUT
!******************************************************************************

CALL CHECK_QUANTITY_VALIDITY
CALL CHECK_POINT_COORDINATE
CALL CHECK_PROFILE_COORDINATE
CALL CHECK_SMOKEVIEW_COORDINATE

CONTAINS

!******************************************!
SUBROUTINE CHECK_QUANTITY_VALIDITY
   INTEGER I, IHIBOUND,ITYPE
   INTEGER :: IMESH
   CHARACTER(60), ALLOCATABLE, DIMENSION(:) :: QUANTITY_NAME
   INTEGER, ALLOCATABLE, DIMENSION(:) :: QUANTITY_INDEX, QUANTITY_IMESH
   CHARACTER(3) :: THREE
   CHARACTER(1000) :: MESSAGE

   IHIBOUND = MAX(GPG%N_POINT_QUANTITIES, GPG%N_PROFILE_QUANTITIES, GPG%N_SMOKEVIEW_QUANTITIES)
   ALLOCATE(QUANTITY_NAME (1:IHIBOUND))
   ALLOCATE(QUANTITY_INDEX(1:IHIBOUND))
   ALLOCATE(QUANTITY_IMESH(1:IHIBOUND))

   ! Check for requested output that could cause segmentation faults 
   DO ITYPE = 1, 3
      IF (ITYPE .EQ. 1) THEN
         IHIBOUND = GPG%N_POINT_QUANTITIES
         QUANTITY_NAME (1:IHIBOUND) = GPG%POINT_QUANTITY      (1:IHIBOUND)
         QUANTITY_INDEX(1:IHIBOUND) = GPG%POINT_QUANTITY_INDEX(1:IHIBOUND)
         QUANTITY_IMESH(1:IHIBOUND) = GPG%POINT_IMESH         (1:IHIBOUND)
      ENDIF

      IF (ITYPE .EQ. 2) THEN
         IHIBOUND = GPG%N_PROFILE_QUANTITIES
         QUANTITY_NAME (1:IHIBOUND) = GPG%PROFILE_QUANTITY      (1:IHIBOUND)
         QUANTITY_INDEX(1:IHIBOUND) = GPG%PROFILE_QUANTITY_INDEX(1:IHIBOUND)
         QUANTITY_IMESH(1:IHIBOUND) = GPG%PROFILE_IMESH         (1:IHIBOUND)
      ENDIF

      IF (ITYPE .EQ. 3) THEN
         IHIBOUND = GPG%N_SMOKEVIEW_QUANTITIES
         QUANTITY_NAME (1:IHIBOUND) = GPG%SMOKEVIEW_QUANTITY      (1:IHIBOUND)
         QUANTITY_INDEX(1:IHIBOUND) = GPG%SMOKEVIEW_QUANTITY_INDEX(1:IHIBOUND)
         QUANTITY_IMESH(1:IHIBOUND) = GPG%SMOKEVIEW_IMESH         (1:IHIBOUND)
      ENDIF

      ! First check to make sure quantity name is valid
      DO I = 1, IHIBOUND
         SELECT CASE(QUANTITY_NAME(I))
            CASE ('TEMPERATURE')
            CASE ('ENTHALPY')
            CASE ('YI')
               IF (QUANTITY_INDEX(I) .GT. SPROP%NSSPEC) THEN 
                  MESSAGE='For output quantity YI, species index cannot be greater than the number of solid species.' 
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
               IF (QUANTITY_INDEX(I) .LT. 0) THEN 
                  MESSAGE='Error : For output quantity YI, solid species index are not specified.' 
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            CASE ('XI')
               IF (QUANTITY_INDEX(I) .GT. SPROP%NSSPEC) THEN 
                  MESSAGE='For output quantity XI, species index cannot begreater than the number of solid species.' 
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
               IF (QUANTITY_INDEX(I) .LT. 0) THEN 
                  MESSAGE='Error : For output quantity XI, solid species index are not specified.' 
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            CASE ('CI')
               IF (QUANTITY_INDEX(I) .GT. SPROP%NSSPEC) THEN 
                  MESSAGE='For output quantity CI, species index cannot be greater than the number of solid species.' 
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
               IF (QUANTITY_INDEX(I) .LT. 0) THEN 
                  MESSAGE='Error : For output quantity CI, solid species index are not specified.' 
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            CASE ('YISUM')
            CASE ('REACTION_RATE_K')
               IF (QUANTITY_INDEX(I) .GT. SPROP%NRXN) THEN 
                  MESSAGE='For output quantity REACTION_RATE_K, reaction index cannot be greater than the number of reactions.' 
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
               IF (QUANTITY_INDEX(I) .LT. 0) THEN 
                  MESSAGE='Error : For output quantity REACTION_RATE_K, reaction index are not specified.' 
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            CASE ('YJ')
               IF (QUANTITY_INDEX(I) .GT. GPROP%NGSPEC) THEN 
                  MESSAGE='For output quantity YJ, species index cannot be greater than the number of gaseous species.' 
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
               IF (QUANTITY_INDEX(I) .LT. 0) THEN 
                  MESSAGE='Error : For output quantity YJ, gas species index are not specified.' 
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            CASE ('CJ')
               IF (QUANTITY_INDEX(I) .GT. GPROP%NGSPEC) THEN 
                  MESSAGE='For output quantity CJ, species index cannot be greater than the number of gaseous species.' 
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
               IF (QUANTITY_INDEX(I) .LT. 0) THEN 
                  MESSAGE='Error : For output quantity CJ, gas species index are not specified.' 
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            CASE ('YJSUM')
            CASE('GAS_DENSITY')
            CASE ('REACTION_RATE_L')
               IF (QUANTITY_INDEX(I) .GT. GPROP%NHGRXN) THEN 
                  MESSAGE='For output quantity REACTION_RATE_L, reaction index cannot be greater than the number of homogeneous gaseous reactions.' 
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
               IF (QUANTITY_INDEX(I) .LT. 0) THEN 
                  MESSAGE='Error : For output quantity REACTION_RATE_L, gas-gas reaction index are not specified.' 
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            CASE ('S')
            CASE ('THERMAL_CONDUCTIVITY_Z')
            CASE ('THERMAL_CONDUCTIVITY_X')
            CASE ('THERMAL_CONDUCTIVITY_Y')
            CASE ('BULK_DENSITY')
            CASE ('SOLID_DENSITY')
               IF (.NOT. GPG%SOLVE_POROSITY ) THEN
                  MESSAGE = 'Error: SOLID_DENSITY is not a valid output quantity unless the pure solid density RS0 of ' // &
                  'condensed species is defined, and at least one of them differs from the bulk density R0.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            CASE ('SPECIFIC_HEAT_CAPACITY')
            CASE ('PRESSURE')
               IF (.NOT. GPG%SOLVE_PRESSURE) THEN 
                  MESSAGE='PRESSURE is not a valid output quantity unless SOLVE_PRESSURE = .TRUE.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            CASE ('MASS_FLUX_TOTAL_Z')
            CASE ('MASS_FLUX_TOTAL_X')
               IF (.NOT. GPG%SOLVE_PRESSURE) THEN 
                  MESSAGE='MASS_FLUX_TOTAL_X is not a valid output quantity unless SOLVE_PRESSURE = .TRUE.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            CASE ('MASS_FLUX_TOTAL_Y')
               IF (.NOT. GPG%SOLVE_PRESSURE) THEN 
                  MESSAGE='MASS_FLUX_TOTAL_Y is not a valid output quantity unless SOLVE_PRESSURE = .TRUE.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            CASE ('GAS_TEMPERATURE')
                  IF ( (.NOT. GPG%SOLVE_GAS_ENERGY) .OR. (GPG%THERMAL_EQUILIBRIUM .AND. GPG%SOLVE_GAS_ENERGY)) THEN 
                     MESSAGE='GAS_TEMPERATURE is not a valid output quantity unless SOLVE_GAS_ENERGY = .TRUE. or in thermal equilibrium mode.'
                     CALL SHUTDOWN_GPYRO(MESSAGE)
                  ENDIF
            CASE ('GAS_ENTHALPY')
               IF ((.NOT. GPG%SOLVE_GAS_ENERGY) .OR. (GPG%THERMAL_EQUILIBRIUM .AND. GPG%SOLVE_GAS_ENERGY)) THEN 
                  MESSAGE='GAS_ENTHALPY is not a valid output quantity unless SOLVE_GAS_ENERGY = .TRUE.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            CASE ('TG-T')
               IF ( (.NOT. GPG%SOLVE_GAS_ENERGY) .OR. (GPG%THERMAL_EQUILIBRIUM .AND. GPG%SOLVE_GAS_ENERGY)) THEN 
                  MESSAGE='TG-T is not a valid output quantity unless SOLVE_GAS_ENERGY = .TRUE.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            CASE ('SHYI')
            CASE ('SHP')
            CASE ('SHM')
            CASE ('POROSITY')
               IF (.NOT. GPG%SOLVE_POROSITY ) THEN
                  MESSAGE = 'Error: POROSITY is not a valid output quantity unless the pure solid density RS0 of ' // &
                  'condensed species is defined, and at least one of them differs from the bulk density R0.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            CASE ('D12')
               IF (.NOT. (GPG%SOLVE_PRESSURE .OR. GPG%SOLVE_GAS_YJ) ) THEN 
                  MESSAGE='D12 is not a valid output quantity unless SOLVE_PRESSURE = .TRUE. .OR. SOLVE_GAS_YJ = .TRUE.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF     
            CASE ('SHGP')
               IF ( (.NOT. GPG%SOLVE_GAS_ENERGY) .OR. (GPG%THERMAL_EQUILIBRIUM .AND. GPG%SOLVE_GAS_ENERGY)) THEN 
                  MESSAGE='SHGP is not a valid output quantity unless SOLVE_GAS_ENERGY = .TRUE.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            CASE ('SHGM')
               IF ( (.NOT. GPG%SOLVE_GAS_ENERGY) .OR. (GPG%THERMAL_EQUILIBRIUM .AND. GPG%SOLVE_GAS_ENERGY)) THEN 
                  MESSAGE='SHGP is not a valid output quantity unless SOLVE_GAS_ENERGY = .TRUE.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            CASE ('QSG')
            CASE ('QSC')               
            CASE ('RYIDZSIGMA')
            CASE ('UNREACTEDNESS')
            CASE ('GOMEGA3')
               IF (QUANTITY_INDEX(I) .GT. SPROP%NSSPEC) THEN 
                  MESSAGE='For output quantity GOMEGA3, species index cannot be greater than the number of solid species.' 
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF   
            CASE ('RE')
               IF (GPG%THERMAL_EQUILIBRIUM .OR. GPG%HCV .GT. 0D0) THEN
                  MESSAGE='RE is not a valid output quantity unless THERMAL_EQUILIBRIUM = .FALSE. and HCV < 0'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            CASE ('NU')
               IF (GPG%THERMAL_EQUILIBRIUM .OR. GPG%HCV .GT. 0D0) THEN
                  MESSAGE='NU is not a valid output quantity unless THERMAL_EQUILIBRIUM = .FALSE. and HCV < 0'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            CASE ('HCV')
               IF (GPG%THERMAL_EQUILIBRIUM .OR. GPG%HCV .GT. 0D0) THEN
                  MESSAGE='HCV is not a valid output quantity unless THERMAL_EQUILIBRIUM = .FALSE. and HCV < 0'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            CASE ('NEEDSBCT')
            CASE ('NEEDSBCB')
            CASE ('NEEDSBCE')
            CASE ('NEEDSBCW')
            CASE ('NEEDSBCN')
            CASE ('NEEDSBCS')
            CASE ('DLTZN')
            CASE ('DLTXN')
            CASE ('DLTYN')
            CASE ('PERMEABILITY_Z')
               IF ( .NOT. GPG%SOLVE_PRESSURE) THEN
                  MESSAGE='Error:  cannot dump PERMEABILITY_Z unless SOLVE_PRESSURE=.TRUE.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            CASE ('PERMEABILITY_X')
               IF ( .NOT. GPG%SOLVE_PRESSURE) THEN
                  MESSAGE='Error:  cannot dump PERMEABILITY_X unless SOLVE_PRESSURE=.TRUE.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF      
            CASE ('PERMEABILITY_Y')
               IF ( .NOT. GPG%SOLVE_PRESSURE) THEN
                  MESSAGE='Error:  cannot dump PERMEABILITY_Y unless SOLVE_PRESSURE=.TRUE.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF 

            CASE ('M/M0')
               IF (ITYPE .NE. 1) THEN
                  MESSAGE='Error, M/M0 is only a point dump quantity.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF

            CASE ('CML')
               IF (ITYPE .NE. 1) THEN
                  MESSAGE='Error, CML is only a point dump quantity.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            
            CASE ('MLR')
               IF (ITYPE .NE. 1) THEN
                  MESSAGE='Error, MLR is only a point dump quantity.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF

               IF (QUANTITY_INDEX(I) .GT. GPROP%NGSPEC) THEN
                  MESSAGE='For MLR, INDEX cannot be greater than # of gas species.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF

               IF (QUANTITY_INDEX(I) .LT. 0) THEN
                  MESSAGE='For MLR, INDEX cannot be less than 0.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF

            CASE ('GGR') ! Gas generation rate
               IF (ITYPE .NE. 1) THEN
                  MESSAGE='Error, GGR is only a point dump quantity.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF

               IF (QUANTITY_INDEX(I) .GT. GPROP%NGSPEC) THEN
                  MESSAGE='For GGR, INDEX cannot be greater than # of gas species.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF

               IF (QUANTITY_INDEX(I) .LT. 0) THEN
                  MESSAGE='For GGR, INDEX cannot be less than 0.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF

            CASE ('HRR')
               IF (ITYPE .NE. 1) THEN
                  MESSAGE='Error, HRR is only a point dump quantity.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF

            CASE ('MDOTPPZ')
               IF (QUANTITY_INDEX(I) .GT. GPROP%NGSPEC) THEN
                  MESSAGE='For MDOTPPZ, INDEX cannot be greater than # of gas species.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF

               IF (QUANTITY_INDEX(I) .LT. 0) THEN
                  MESSAGE='For MDOTPPZ, INDEX cannot be less than 0.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
               IMESH  = QUANTITY_IMESH(I)           
               IF (GPG%NCELLX(IMESH) .GT. 1 .OR. GPG%NCELLY(IMESH) .GT. 1) THEN
                  MESSAGE='MDOTPPZ can only be used for 0D/1D simulations.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF

               IF (ITYPE .EQ. 2) THEN
                  IF (GPG%PROFILE_DIRECTION(I) .NE. 'z' .AND. GPG%PROFILE_DIRECTION(I) .NE. 'Z') THEN
                     MESSAGE='For MDOTPPZ, only valid PROFILE_DIRECTION is z'
                     CALL SHUTDOWN_GPYRO(MESSAGE)
                  ENDIF
               ENDIF

               IF (ITYPE .EQ. 3) THEN
                  MESSAGE='Cannot dump Smokeview plane for MDOTPPZ'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            CASE ('MPPI')
               IF (ITYPE .NE. 1) THEN
                  MESSAGE='Error, MPPI is only a point dump quantity.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF

            CASE ('THICKNESS')
               IF (ITYPE .NE. 1) THEN
                  MESSAGE='Error, THICKNESS is only a point dump quantity.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF

            CASE ('DT')
               IF (ITYPE .NE. 1) THEN
                  MESSAGE='Error, DT is only a point dump quantity.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF

            CASE ('N_ITERATIONS')
               IF (ITYPE .NE. 1) THEN
                  MESSAGE='Error, N_ITERATIONS is only a point dump quantity.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF

            CASE ('TOTAL_MASS')
               IF (ITYPE .NE. 1) THEN
                  MESSAGE='Error, TOTAL_MASS is only a point dump quantity.'
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
               IF (QUANTITY_INDEX(I) .GT. SPROP%NSSPEC) THEN 
                  MESSAGE='For output quantity TOTAL_MASS, species index cannot be greater than the number of solid species.' 
                  CALL SHUTDOWN_GPYRO(MESSAGE)
               ENDIF
            !!! DEAD OPTION !!!
            CASE ('MDOTPPDARCY')
               MESSAGE='MDOTPPDARCY is no longer a valid output quantity. Try MASS_FLUX_TOTAL_Z or MASS_FLUX_TOTAL_X or MASS_FLUX_TOTAL_Y.'
               CALL SHUTDOWN_GPYRO(MESSAGE)
            CASE('THERMAL_CONDUCTIVITY')
               MESSAGE='THERMAL_CONDUCTIVITY is not a valid output quantity. Specify THERMAL_CONDUCTIVITY_Z, THERMAL_CONDUCTIVITY_Y, or THERMAL_CONDUCTIVITY_X'
               CALL SHUTDOWN_GPYRO(MESSAGE)
            CASE('PERMEABILITY')
               MESSAGE='PERMEABILITY is not a valid output quantity. Specify PERMEABILITY_Z, PERMEABILITY_Y, or PERMEABILITY_X'
               CALL SHUTDOWN_GPYRO(MESSAGE)
            CASE DEFAULT
               WRITE(THREE,'(I3.3)') I
               IF (ITYPE .EQ. 1) MESSAGE = 'Error, quantity ' // TRIM(QUANTITY_NAME(I)) // ' not valid for Point dump. INDEX = ' // THREE
               IF (ITYPE .EQ. 2) MESSAGE = 'Error, quantity ' // TRIM(QUANTITY_NAME(I)) // ' not valid for profile dump. INDEX = ' // THREE
               IF (ITYPE .EQ. 3) MESSAGE = 'Error, quantity ' // TRIM(QUANTITY_NAME(I)) // ' not valid for Smokeview dump. INDEX = ' // THREE
               CALL SHUTDOWN_GPYRO(MESSAGE)
         END SELECT
      ENDDO
   ENDDO !ITYPE

   DEALLOCATE(QUANTITY_NAME)
   DEALLOCATE(QUANTITY_INDEX)

END SUBROUTINE CHECK_QUANTITY_VALIDITY

!******************************************!
SUBROUTINE CHECK_POINT_COORDINATE
   INTEGER I, IMESH, NCELLZ, NCELLX, NCELLY
   REAL(EB):: ZDIM, XDIM, YDIM
   CHARACTER(3) :: THREE
   CHARACTER(1000) :: MESSAGE
   ! Check point dumps:
   DO I = 1, GPG%N_POINT_QUANTITIES
      IMESH  = GPG%POINT_IMESH(I)
      IF (IMESH .EQ. 0) IMESH=1
      NCELLZ = GPG%NCELLZ(IMESH)
      NCELLX = GPG%NCELLX(IMESH)
      NCELLY = GPG%NCELLY(IMESH)
      ZDIM   = GPG%ZDIM(IMESH) 
      XDIM   = GPG%XDIM(IMESH) 
      YDIM   = GPG%YDIM(IMESH) 
   
      IF (NCELLZ .GT. 1) THEN
         IF (GPG%POINT_Z(I) .LT. 0. .OR. GPG%POINT_Z(I) .GT. ZDIM) THEN
            WRITE(THREE,'(I3.3)') I
            MESSAGE='Error:  Point dump ' // THREE // ' z coordinate is not between 0 and zdim'
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF
   
      IF (NCELLX .GT. 1) THEN
         IF (GPG%POINT_X(I) .LT. 0. .OR. GPG%POINT_X(I) .GT. XDIM) THEN
            WRITE(THREE,'(I3.3)') I
            MESSAGE='Error:  Point dump ' // THREE // ' x coordinate is not between 0 and xdim'
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF
   
      IF (NCELLY .GT. 1) THEN
         IF (GPG%POINT_Y(I) .LT. 0. .OR. GPG%POINT_Y(I) .GT. YDIM) THEN
            WRITE(THREE,'(I3.3)') I
            MESSAGE='Error:  Point dump ' // THREE // ' y coordinate is not between 0 and ydim'
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF
   
   ENDDO
   
END SUBROUTINE CHECK_POINT_COORDINATE

!******************************************!
SUBROUTINE CHECK_PROFILE_COORDINATE
   INTEGER I, IMESH, NCELLZ, NCELLX, NCELLY
   REAL(EB):: ZDIM, XDIM, YDIM
   CHARACTER(3) :: THREE
   CHARACTER(1000) :: MESSAGE

   ! Check profile dumps:
   DO I = 1, GPG%N_PROFILE_QUANTITIES
      IMESH  = GPG%PROFILE_IMESH(I)
      IF (IMESH .EQ. 0) IMESH=1
      NCELLZ = GPG%NCELLZ(IMESH)
      NCELLX = GPG%NCELLX(IMESH)
      NCELLY = GPG%NCELLY(IMESH)
      ZDIM   = GPG%ZDIM(IMESH) 
      XDIM   = GPG%XDIM(IMESH) 
      YDIM   = GPG%YDIM(IMESH) 
   
      IF (GPG%PROFILE_DIRECTION(I) .NE. 'z' .AND. GPG%PROFILE_DIRECTION(I) .NE. 'Z' .AND. & 
         GPG%PROFILE_DIRECTION(I) .NE. 'x' .AND. GPG%PROFILE_DIRECTION(I) .NE. 'X' .AND. &
         GPG%PROFILE_DIRECTION(I) .NE. 'y' .AND. GPG%PROFILE_DIRECTION(I) .NE. 'Y') THEN
         MESSAGE='Error:  For profile dumps, set PROFILE_DIRECTION to one of z, x, or y'
         CALL SHUTDOWN_GPYRO(MESSAGE)
      ENDIF
      IF (GPG%PROFILE_DIRECTION(I) .EQ. 'Z') GPG%PROFILE_DIRECTION(I)='z' 
      IF (GPG%PROFILE_DIRECTION(I) .EQ. 'X') GPG%PROFILE_DIRECTION(I)='x' 
      IF (GPG%PROFILE_DIRECTION(I) .EQ. 'Y') GPG%PROFILE_DIRECTION(I)='y'
   
      IF (GPG%PROFILE_DIRECTION(I) .EQ. 'x' .AND. NCELLX .EQ. 1) THEN
         MESSAGE='Error:  cannot have x-direction profile dump with one cell in x-direction'
         CALL SHUTDOWN_GPYRO(MESSAGE)
      ENDIF
   
      IF (GPG%PROFILE_DIRECTION(I) .EQ. 'y' .AND. NCELLY .EQ. 1) THEN
         MESSAGE='Error:  cannot have y-direction profile dump with one cell in y-direction'
         CALL SHUTDOWN_GPYRO(MESSAGE)
      ENDIF
      
      IF (GPG%PROFILE_DIRECTION(I) .EQ. 'z') THEN
         IF (NCELLX .GT. 1) THEN
            IF (GPG%PROFILE_COORD1(I) .LT. 0. .OR. GPG%PROFILE_COORD1(I) .GT. XDIM) THEN
               WRITE(THREE,'(I3.3)') I
               MESSAGE='Error:  Profile dump ' // THREE // ' x coordinate is not between 0 and xdim'
               CALL SHUTDOWN_GPYRO(MESSAGE)
            ENDIF
         ENDIF
         IF (NCELLY .GT. 1) THEN
            IF (GPG%PROFILE_COORD2(I) .LT. 0. .OR. GPG%PROFILE_COORD2(I) .GT. YDIM) THEN
               WRITE(THREE,'(I3.3)') I
               MESSAGE='Error:  Profile dump ' // THREE // ' y coordinate is not between 0 and ydim'
               CALL SHUTDOWN_GPYRO(MESSAGE)
            ENDIF
         ENDIF
      ENDIF
   
      IF (GPG%PROFILE_DIRECTION(I) .EQ. 'x') THEN
         IF (NCELLY .GT. 1) THEN
            IF (GPG%PROFILE_COORD1(I) .LT. 0. .OR. GPG%PROFILE_COORD1(I) .GT. YDIM) THEN
               WRITE(THREE,'(I3.3)') I
               MESSAGE='Error:  Profile dump ' // THREE // ' y coordinate is not between 0 and ydim'
               CALL SHUTDOWN_GPYRO(MESSAGE)
            ENDIF
         ENDIF
         IF (NCELLZ .GT. 1) THEN
            IF (GPG%PROFILE_COORD2(I) .LT. 0. .OR. GPG%PROFILE_COORD2(I) .GT. ZDIM) THEN
               WRITE(THREE,'(I3.3)') I
               MESSAGE='Error:  Profile dump ' // THREE // ' z coordinate is not between 0 and zdim'
               CALL SHUTDOWN_GPYRO(MESSAGE)
            ENDIF
         ENDIF
      ENDIF
   
      IF (GPG%PROFILE_DIRECTION(I) .EQ. 'y') THEN
         IF (NCELLX .GT. 1) THEN
            IF (GPG%PROFILE_COORD1(I) .LT. 0. .OR. GPG%PROFILE_COORD1(I) .GT. XDIM) THEN
               WRITE(THREE,'(I3.3)') I
               MESSAGE='Error:  Profile dump ' // THREE // ' x coordinate is not between 0 and xdim'
               CALL SHUTDOWN_GPYRO(MESSAGE)
            ENDIF
         ENDIF
         IF (NCELLZ .GT. 1) THEN
            IF (GPG%PROFILE_COORD2(I) .LT. 0. .OR. GPG%PROFILE_COORD2(I) .GT. ZDIM) THEN
               WRITE(THREE,'(I3.3)') I
               MESSAGE='Error:  Profile dump ' // THREE // ' z coordinate is not between 0 and zdim'
               CALL SHUTDOWN_GPYRO(MESSAGE)
            ENDIF
         ENDIF
      ENDIF
   ENDDO
END SUBROUTINE CHECK_PROFILE_COORDINATE

!******************************************!
SUBROUTINE CHECK_SMOKEVIEW_COORDINATE
   INTEGER N, IMESH, NCELLZ, NCELLX, NCELLY
   REAL(EB):: ZDIM, XDIM, YDIM
   CHARACTER(1000) :: MESSAGE

   !Check Smokeview dumps:
   DO N = 1, GPG%N_SMOKEVIEW_QUANTITIES
      IMESH  = GPG%SMOKEVIEW_IMESH(N)
      IF (IMESH .EQ. 0) IMESH=1
      NCELLZ = GPG%NCELLZ(IMESH)
      NCELLX = GPG%NCELLX(IMESH)
      NCELLY = GPG%NCELLY(IMESH)
      ZDIM   = GPG%ZDIM(IMESH) 
      XDIM   = GPG%XDIM(IMESH) 
      YDIM   = GPG%YDIM(IMESH) 
      IF (GPG%SMOKEVIEW_PLANE(N) .EQ. 'xz' .OR. GPG%SMOKEVIEW_PLANE(N) .EQ. 'XZ' .OR. &
            GPG%SMOKEVIEW_PLANE(N) .EQ. 'zx' .OR. GPG%SMOKEVIEW_PLANE(N) .EQ. 'ZX') THEN !y=const plane
         GPG%SMOKEVIEW_PLANE(N)='xz'
      ENDIF
   
      IF (GPG%SMOKEVIEW_PLANE(N) .EQ. 'yz' .OR. GPG%SMOKEVIEW_PLANE(N) .EQ. 'YZ' .OR. &
            GPG%SMOKEVIEW_PLANE(N) .EQ. 'zy' .OR. GPG%SMOKEVIEW_PLANE(N) .EQ. 'ZY') THEN !x=const plane
         GPG%SMOKEVIEW_PLANE(N)='yz'   
      ENDIF
   
      IF (GPG%SMOKEVIEW_PLANE(N) .EQ. 'xy' .OR. GPG%SMOKEVIEW_PLANE(N) .EQ. 'XY' .OR. &
         GPG%SMOKEVIEW_PLANE(N) .EQ. 'yx' .OR. GPG%SMOKEVIEW_PLANE(N) .EQ. 'YX') THEN !z=const plane
         GPG%SMOKEVIEW_PLANE(N)='xy'
      ENDIF
   
      IF (GPG%SMOKEVIEW_PLANE(N) .NE. 'xz' .AND. GPG%SMOKEVIEW_PLANE(N) .NE. 'yz' .AND. GPG%SMOKEVIEW_PLANE(N) .NE. 'xy') THEN
         MESSAGE='Error:  For Smokeview output, set SMOKEVIEW_PLANE to one of xy, yz, or xz'
         CALL SHUTDOWN_GPYRO(MESSAGE)
      ENDIF
   
      IF (GPG%SMOKEVIEW_PLANE(N) .EQ. 'xz') THEN
         IF (NCELLZ .EQ. 1 .OR. NCELLX .EQ. 1) THEN
            MESSAGE='Error:  Cannot have xz SMOKEVIEW_PLANE with only one cell in x or z direction'
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF
   
      IF (GPG%SMOKEVIEW_PLANE(N) .EQ. 'xy') THEN
         IF (NCELLX .EQ. 1 .OR. NCELLY .EQ. 1) THEN
            MESSAGE='Error:  Cannot have xy SMOKEVIEW_PLANE with only one cell in x or y direction'
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF
   
      IF (GPG%SMOKEVIEW_PLANE(N) .EQ. 'yz') THEN
         IF (NCELLY .EQ. 1 .OR. NCELLZ .EQ. 1) THEN
            MESSAGE='Error:  Cannot have yz SMOKEVIEW_PLANE with only one cell in y or z direction'
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF   
   ENDDO
END SUBROUTINE CHECK_SMOKEVIEW_COORDINATE

!******************************************************************************
END SUBROUTINE CHECK_OUTPUT
!******************************************************************************


!******************************************************************************	
END MODULE GPYRO_CHECK
!******************************************************************************