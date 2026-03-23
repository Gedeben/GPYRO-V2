MODULE GPYRO_INIT

USE PREC
USE GPYRO_VARS
USE GPYRO_FUNCS
USE GPYRO_IO
USE GPYRO_MASS, ONLY : RHOOFT

IMPLICIT NONE

CONTAINS

!******************************************************************************	
SUBROUTINE CHECK_GPYRO(IRANK)
!******************************************************************************	

INTEGER ,INTENT(IN) :: IRANK
INTEGER I, J, N, IHIBOUND,ITYPE,IRXN,ICINDEX,ISPEC,IOBST
INTEGER :: ICASE,IMESH, NCELLZ, NCELLX, NCELLY, NMESH
REAL(EB):: ZDIM, XDIM, YDIM
CHARACTER(60), ALLOCATABLE, DIMENSION(:) :: QUANTITY_NAME
CHARACTER(2) :: TWO
CHARACTER(3) :: THREE
CHARACTER(300) :: MESSAGE
INTEGER, ALLOCATABLE, DIMENSION(:) :: QUANTITY_INDEX, QUANTITY_IMESH
REAL(EB) :: SUMVAL

CHARACTER(LEN=12) :: FACENAME(6)
INTEGER :: IFACE, ISURF, BC_IDX, ICNUM
LOGICAL :: FOUND
CHARACTER(LEN=20) ::  STR_IMESH, STR_BCIDX, STR_FACEIDX


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
         CASE ('XI')
         CASE ('CI')
         CASE ('YISUM')
         CASE ('REACTION_RATE_K')
         CASE ('YJ')
         CASE ('CJ') 
         CASE ('YJSUM') 
         CASE ('REACTION_RATE_L')
         CASE ('S')
         CASE ('THERMAL_CONDUCTIVITY_Z')
         CASE ('THERMAL_CONDUCTIVITY_X')
         CASE ('THERMAL_CONDUCTIVITY_Y')
         CASE ('BULK_DENSITY')
         CASE ('SOLID_DENSITY')
         CASE ('SPECIFIC_HEAT_CAPACITY')
         CASE ('PRESSURE')
         CASE ('MASS_FLUX_TOTAL_Z')
         CASE ('MASS_FLUX_TOTAL_X')
         CASE ('MASS_FLUX_TOTAL_Y')
         CASE ('GAS_TEMPERATURE')
         CASE ('GAS_ENTHALPY')
         CASE ('TG-T')
         CASE ('SHYI')
         CASE ('SHP')
         CASE ('SHM')
         CASE ('POROSITY')
         CASE ('D12')
         CASE ('SHGP')
         CASE ('SHGM')
         CASE ('QSG')
         CASE ('RYIDZSIGMA')
         CASE ('UNREACTEDNESS')
         CASE ('GOMEGA3')
         CASE ('RE')
         CASE ('NU')
         CASE ('HCV')
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
            NCELLX = GPG%NCELLX(IMESH)
            NCELLY = GPG%NCELLY(IMESH)
            
            IF (NCELLX .GT. 1 .OR. NCELLY .GT. 1) THEN
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
            
         !CASE ('FRONT_FACE_DIFFUSIVE_FLUX')
         !   IF (ITYPE .NE. 1) THEN
         !      MESSAGE='Error, FRONT_FACE_DIFFUSIVE_FLUX is only a point dump quantity.'
         !      CALL SHUTDOWN_GPYRO(MESSAGE)
         !   ENDIF

         !CASE ('FRONT_FACE_CONVECTIVE_FLUX')
         !   IF (ITYPE .NE. 1) THEN
         !      MESSAGE='Error, FRONT_FACE_CONVECTIVE_FLUX is only a point dump quantity.'
         !      CALL SHUTDOWN_GPYRO(MESSAGE)
         !   ENDIF

         !CASE ('FRONT_FACE_MASS_FLUX')
         !   IF (ITYPE .NE. 1) THEN
         !      MESSAGE='Error, FRONT_FACE_MASS_FLUX is only a point dump quantity.'
         !      CALL SHUTDOWN_GPYRO(MESSAGE)
         !   ENDIF

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

         CASE DEFAULT
            WRITE(THREE,'(I3.3)') I
            IF (ITYPE .EQ. 1) MESSAGE = 'Error, quantity ' // TRIM(QUANTITY_NAME(I)) // ' not valid for Point dump. INDEX = ' // THREE
            IF (ITYPE .EQ. 2) MESSAGE = 'Error, quantity ' // TRIM(QUANTITY_NAME(I)) // ' not valid for profile dump. INDEX = ' // THREE
            IF (ITYPE .EQ. 3) MESSAGE = 'Error, quantity ' // TRIM(QUANTITY_NAME(I)) // ' not valid for Smokeview dump. INDEX = ' // THREE
            CALL SHUTDOWN_GPYRO(MESSAGE)
      END SELECT
   ENDDO

   DO I = 1, IHIBOUND

      IF (QUANTITY_NAME(I) .EQ. 'YI') THEN
         IF (QUANTITY_INDEX(I) .GT. SPROP%NSSPEC) THEN 
            MESSAGE='For output quantity YI, species index cannot be greater than the number of solid species.' 
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF

      IF (QUANTITY_NAME(I) .EQ. 'XI') THEN
         IF (QUANTITY_INDEX(I) .GT. SPROP%NSSPEC) THEN 
            MESSAGE='For output quantity XI, species index cannot begreater than the number of solid species.' 
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF

      IF (QUANTITY_NAME(I) .EQ. 'CI') THEN
         IF (QUANTITY_INDEX(I) .GT. SPROP%NSSPEC) THEN 
            MESSAGE='For output quantity CI, species index cannot be greater than the number of solid species.' 
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF

      IF (QUANTITY_NAME(I) .EQ. 'TOTAL_MASS') THEN
         IF (QUANTITY_INDEX(I) .GT. SPROP%NSSPEC) THEN 
            MESSAGE='For output quantity TOTAL_MASS, species index cannot be greater than the number of solid species.' 
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF

      IF (QUANTITY_NAME(I) .EQ. 'REACTION_RATE_K') THEN
         IF (QUANTITY_INDEX(I) .GT. SPROP%NRXN) THEN 
            MESSAGE='For output quantity REACTION_RATE_K, reaction index cannot be greater than the number of reactions.' 
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF

      IF (QUANTITY_NAME(I) .EQ. 'YJ') THEN
         IF (QUANTITY_INDEX(I) .GT. GPROP%NGSPEC) THEN 
            MESSAGE='For output quantity YJ, species index cannot be greater than the number of gaseous species.' 
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF

      IF (QUANTITY_NAME(I) .EQ. 'CJ') THEN
         IF (QUANTITY_INDEX(I) .GT. GPROP%NGSPEC) THEN 
            MESSAGE='For output quantity CJ, species index cannot be greater than the number of gaseous species.' 
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF

      IF (QUANTITY_NAME(I) .EQ. 'REACTION_RATE_L') THEN
         IF (QUANTITY_INDEX(I) .GT. GPROP%NHGRXN) THEN 
            MESSAGE='For output quantity REACTION_RATE_L, reaction index cannot be greater than the number of homogeneous gaseous reactions.' 
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF

      IF (QUANTITY_NAME(I) .EQ. 'PRESSURE') THEN
         IF (.NOT. GPG%SOLVE_PRESSURE) THEN 
            MESSAGE='PRESSURE is not a valid output quantity unless SOLVE_PRESSURE = .TRUE.'
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF

      IF (QUANTITY_NAME(I) .EQ. 'MDOTPPDARCY') THEN
         MESSAGE='MDOTPPDARCY is no longer a valid output quantity. Try MASS_FLUX_TOTAL_Z or MASS_FLUX_TOTAL_X or MASS_FLUX_TOTAL_Y.'
         CALL SHUTDOWN_GPYRO(MESSAGE)
      ENDIF

      IF (QUANTITY_NAME(I) .EQ. 'THERMAL_CONDUCTIVITY') THEN
         MESSAGE='THERMAL_CONDUCTIVITY is not a valid output quantity. Specify THERMAL_CONDUCTIVITY_Z, THERMAL_CONDUCTIVITY_Y, or THERMAL_CONDUCTIVITY_X'
         CALL SHUTDOWN_GPYRO(MESSAGE)
      ENDIF

      IF (QUANTITY_NAME(I) .EQ. 'PERMEABILITY') THEN
         MESSAGE='PERMEABILITY is not a valid output quantity. Specify PERMEABILITY_Z, PERMEABILITY_Y, or PERMEABILITY_X'
         CALL SHUTDOWN_GPYRO(MESSAGE)
      ENDIF


      IF (QUANTITY_NAME(I) .EQ. 'MASS_FLUX_TOTAL_X') THEN
         IF (.NOT. GPG%SOLVE_PRESSURE) THEN 
            MESSAGE='MASS_FLUX_TOTAL_X is not a valid output quantity unless SOLVE_PRESSURE = .TRUE.'
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF

      IF (QUANTITY_NAME(I) .EQ. 'MASS_FLUX_TOTAL_Y') THEN
         IF (.NOT. GPG%SOLVE_PRESSURE) THEN 
            MESSAGE='MASS_FLUX_TOTAL_Y is not a valid output quantity unless SOLVE_PRESSURE = .TRUE.'
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF

      IF (QUANTITY_NAME(I) .EQ. 'GAS_TEMPERATURE') THEN
         IF ( (.NOT. GPG%SOLVE_GAS_ENERGY) .OR. (GPG%THERMAL_EQUILIBRIUM .AND. GPG%SOLVE_GAS_ENERGY)) THEN 
            MESSAGE='GAS_TEMPERATURE is not a valid output quantity unless SOLVE_GAS_ENERGY = .TRUE. or in thermal equilibrium mode.'
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF

      IF (QUANTITY_NAME(I) .EQ. 'GAS_ENTHALPY') THEN
         IF ( (.NOT. GPG%SOLVE_GAS_ENERGY) .OR. (GPG%THERMAL_EQUILIBRIUM .AND. GPG%SOLVE_GAS_ENERGY)) THEN 
            MESSAGE='GAS_ENTHALPY is not a valid output quantity unless SOLVE_GAS_ENERGY = .TRUE.'
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF

      IF (QUANTITY_NAME(I) .EQ. 'TG-T') THEN
         IF ( (.NOT. GPG%SOLVE_GAS_ENERGY) .OR. (GPG%THERMAL_EQUILIBRIUM .AND. GPG%SOLVE_GAS_ENERGY)) THEN 
            MESSAGE='TG-T is not a valid output quantity unless SOLVE_GAS_ENERGY = .TRUE.'
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF

      IF (QUANTITY_NAME(I) .EQ. 'PERMEABILITY') THEN
         IF (.NOT. GPG%SOLVE_PRESSURE) THEN 
            MESSAGE='PERMEABILITY is not a valid output quantity unless SOLVE_PRESSURE = .TRUE.'
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF

      IF (QUANTITY_NAME(I) .EQ. 'POROSITY') THEN
         IF (.NOT. GPG%SOLVE_POROSITY ) THEN
            MESSAGE = 'Error: POROSITY is not a valid output quantity unless the pure solid density RS0 of ' // &
            'condensed species is defined, and at least one of them differs from the bulk density R0.'
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF

      IF (QUANTITY_NAME(I) .EQ. 'SOLID_DENSITY') THEN
         IF (.NOT. GPG%SOLVE_POROSITY ) THEN
            MESSAGE = 'Error: SOLID_DENSITY is not a valid output quantity unless the pure solid density RS0 of ' // &
            'condensed species is defined, and at least one of them differs from the bulk density R0.'
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF

      IF (QUANTITY_NAME(I) .EQ. 'D12') THEN
         IF (.NOT. (GPG%SOLVE_PRESSURE .OR. GPG%SOLVE_GAS_YJ) ) THEN 
            MESSAGE='D12 is not a valid output quantity unless SOLVE_PRESSURE = .TRUE. .OR. SOLVE_GAS_YJ = .TRUE.'
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF

      IF (QUANTITY_NAME(I) .EQ. 'SHGP') THEN
         IF ( (.NOT. GPG%SOLVE_GAS_ENERGY) .OR. (GPG%THERMAL_EQUILIBRIUM .AND. GPG%SOLVE_GAS_ENERGY)) THEN 
            MESSAGE='SHGP is not a valid output quantity unless SOLVE_GAS_ENERGY = .TRUE.'
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF

      IF (QUANTITY_NAME(I) .EQ. 'SHGM') THEN
         IF ( (.NOT. GPG%SOLVE_GAS_ENERGY) .OR. (GPG%THERMAL_EQUILIBRIUM .AND. GPG%SOLVE_GAS_ENERGY)) THEN 
            MESSAGE='SHGM is not a valid output quantity unless SOLVE_GAS_ENERGY = .TRUE.'
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF

      IF (QUANTITY_NAME(I) .EQ. 'GOMEGA3') THEN
         IF (QUANTITY_INDEX(I) .GT. SPROP%NSSPEC) THEN 
            MESSAGE='For output quantity GOMEGA3, species index cannot be greater than the number of solid species.' 
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF

      IF (GPG%THERMAL_EQUILIBRIUM .OR. GPG%HCV .GT. 0D0) THEN
         IF (QUANTITY_NAME(I) .EQ. 'RE') THEN
            MESSAGE='RE is not a valid output quantity unless THERMAL_EQUILIBRIUM = .FALSE. and HCV < 0'
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF

         IF (QUANTITY_NAME(I) .EQ. 'NU') THEN
            MESSAGE='NU is not a valid output quantity unless THERMAL_EQUILIBRIUM = .FALSE. and HCV < 0'
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
         
         IF (QUANTITY_NAME(I) .EQ. 'HCV') THEN
            MESSAGE='HCV is not a valid output quantity unless THERMAL_EQUILIBRIUM = .FALSE. and HCV < 0'
            CALL SHUTDOWN_GPYRO(MESSAGE)
         ENDIF
      ENDIF

   ENDDO !IHIBOUND
ENDDO !ITYPE

DEALLOCATE(QUANTITY_NAME)
DEALLOCATE(QUANTITY_INDEX)


NMESH = GPG%NUM_GPYRO_MESHES

! Check 0D/1D/2D/3D initialization
DO IMESH = 1, NMESH

   NCELLZ = GPG%NCELLZ(IMESH)
   NCELLX = GPG%NCELLX(IMESH)
   NCELLY = GPG%NCELLY(IMESH)

   IF (NCELLZ .EQ. 1 .AND. NCELLX .GT. 1 .AND. NCELLY .EQ. 1) THEN
      MESSAGE='Error:  For 1D simulations set NCELLY = 1, NCELLX = 1, and NCELLZ > 1'
      CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF
   
   IF (NCELLZ .EQ. 1 .AND. NCELLY .GT. 1 .AND. NCELLX .EQ. 1) THEN
      MESSAGE='Error:  For 1D simulations set NCELLY = 1, NCELLX = 1, and NCELLZ > 1'
      CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF
   
   IF (NCELLY .GT. 1 .AND. NCELLX .EQ. 1) THEN
      MESSAGE='Error:  for 2D simulations set NCELLY = 1 and NCELLX > 1'
      CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF
ENDDO





! Check mesh boundary condition definitions 
! Face order and names
FACENAME(1) = 'west (-x)'
FACENAME(2) = 'east (+x)'
FACENAME(3) = 'south (-y)'
FACENAME(4) = 'north (+y)'
FACENAME(5) = 'top (+z)'
FACENAME(6) = 'bottom (-z)'

DO IMESH = 1, GPG%NUM_GPYRO_MESHES

   DO IFACE = 1, 6

      ! Skip unused directions depending on mesh dimension
      IF ( (GPG%NCELLX(IMESH) .LE. 1) .AND. (IFACE == 1 .OR. IFACE == 2) ) CYCLE
      IF ( (GPG%NCELLY(IMESH) .LE. 1) .AND. (IFACE == 3 .OR. IFACE == 4) ) CYCLE
      IF ( (GPG%NCELLZ(IMESH) .LE. 1) .AND. (IFACE == 5 .OR. IFACE == 6) ) CYCLE

      BC_IDX= GPG%DEFAULT_SURF_IDX(IMESH, IFACE)

      ! Convert integers to strings
      WRITE(STR_BCIDX,'(I0)') BC_IDX
      WRITE(STR_IMESH,'(I0)') IMESH
      WRITE(STR_FACEIDX, '(I0)') IFACE
      
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
      ENDDO
      
      IF (.NOT. FOUND) THEN
         MESSAGE = 'Error: Invalid boundary condition for ' // TRIM(FACENAME(IFACE)) //&
          ' face of mesh ' // TRIM(STR_IMESH) // &
          '. The Boundary index ' // TRIM(STR_BCIDX) // ' is not defined.' //&
          ' Please change DEFAULT_SURF_IDX('// TRIM(STR_IMESH)//','// TRIM(STR_FACEIDX)// ').'

         CALL SHUTDOWN_GPYRO(MESSAGE)
      ENDIF
      
   ENDDO

   ICNUM = GPG%DEFAULT_IC(IMESH)

   ! Convert integers to strings
   WRITE(STR_BCIDX,'(I0)') ICNUM
   WRITE(STR_IMESH,'(I0)') IMESH
   
   ! Check if BC index is defined (positive)
   IF (ICNUM .GT. GPG%NIC) THEN
      MESSAGE = 'Error: Invalid initial condition for mesh ' // TRIM(STR_IMESH) // &
      '. The initial condition index ' // TRIM(STR_BCIDX) // ' is not defined.' //&
      ' Please change DEFAULT_IC('// TRIM(STR_IMESH)//').'

     CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF

ENDDO


DO IOBST = 1, GPG%NOBST
   IMESH= GPG%GEOM(IOBST)%IMESH 

   DO IFACE = 1, 6

      ! Skip unused directions depending on mesh dimension
      IF ( (GPG%NCELLX(IMESH) .LE. 1) .AND. (IFACE == 1 .OR. IFACE == 2) ) CYCLE
      IF ( (GPG%NCELLY(IMESH) .LE. 1) .AND. (IFACE == 3 .OR. IFACE == 4) ) CYCLE
      IF ( (GPG%NCELLZ(IMESH) .LE. 1) .AND. (IFACE == 5 .OR. IFACE == 6) ) CYCLE

      BC_IDX= GPG%GEOM(IOBST)%SURF_IDX(IFACE)
      IF (BC_IDX .LE. 0) CYCLE


      ! Convert integers to strings
      WRITE(STR_BCIDX,'(I0)') BC_IDX
      WRITE(STR_IMESH,'(I0)') IOBST
      WRITE(STR_FACEIDX, '(I0)') IFACE
      ! Check if BC index exists in SURF_IDX list
      FOUND = .FALSE.
      DO ISURF = 1, GPG%NSURF_IDX
         IF (BC_IDX == GPG%ALLBC(ISURF)%SURF_IDX) THEN
            FOUND = .TRUE.
            EXIT
         ENDIF
      ENDDO
      
      IF (.NOT. FOUND) THEN
         MESSAGE = 'Error: Invalid boundary condition for ' // TRIM(FACENAME(IFACE)) //&
          ' face of OBST ' // TRIM(STR_IMESH) // &
          '. The Boundary index ' // TRIM(STR_BCIDX) // ' is not defined.' //&
          ' Please change SURF_IDX2D('// TRIM(STR_IMESH)//','// TRIM(STR_FACEIDX)// ').'

         CALL SHUTDOWN_GPYRO(MESSAGE)
      ENDIF
      
   ENDDO



   ICNUM = GPG%GEOM(IOBST)%ICNUM
   ! Convert integers to strings
   WRITE(STR_BCIDX,'(I0)') ICNUM
   WRITE(STR_IMESH,'(I0)') IOBST
   
   ! Check if BC index is defined (positive)
   IF (ICNUM .GT. GPG%NIC) THEN
      MESSAGE = 'Error: Invalid initial condition for OBST ' // TRIM(STR_IMESH) // &
      '. The initial condition index ' // TRIM(STR_BCIDX) // ' is not defined.' //&
      ' Please change ICNUM('// TRIM(STR_IMESH)//').'

     CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF

ENDDO







! Check quantities to dump

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



DO ICASE = 1, GPG%NCASES
   IMESH  = GPG%IMESH (ICASE)
   NCELLZ = GPG%NCELLZ(IMESH)
   NCELLX = GPG%NCELLX(IMESH)
   NCELLY = GPG%NCELLY(IMESH)
   IF(GPG%ZEROD(ICASE) .AND. (NCELLX .GT. 1 .OR. NCELLY .GT. 1 .OR. NCELLZ .GT. 1)) THEN
      MESSAGE='Error:  For meshes corresponding to cases where ZEROD=.TRUE., ensure that NCELLZ=NCELLX=NCELLY=1.'
      CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF
ENDDO

DO ISPEC = 1, GPROP%NGSPEC
   IF (TRIM(GPROP%NAME(ISPEC)) .EQ. 'null') THEN
      MESSAGE='Error:  Ensure all gaseous species are defined (and check that # of defined species matches NGSPEC)'
      CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF
ENDDO

IF (GPROP%NHGRXN .GT. 0) THEN
   IF ( (.NOT. GPG%SOLVE_GAS_ENERGY) ) THEN
      MESSAGE='Set SOLVE_GAS_ENERGY = .TRUE. when using homogeneous gaseous reactions (NHGRXN > 0).'
      CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF
ENDIF

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

! Check initial condensed-phase species mass fractions:
DO ICINDEX = 1, GPG%NIC
   SUMVAL = 0D0
   DO I = 1, SPROP%NSSPEC
      SUMVAL = SUMVAL + GPG%INITIAL_CONDITIONS(ICINDEX)%YI0(I)
   ENDDO
   IF (SUMVAL .LT. 0.999999 .OR. SUMVAL .GT. 1.000001) THEN
      WRITE(THREE,'(I3.3)') ICINDEX
      MESSAGE='Initial condensed-phase mass fractions do not sum to 1.0 for IC ' // THREE // &
              '. Check IC worksheet and make sure you have specified initial conditions for all condensed-phase species.'
      CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF

ENDDO

DO I = 1, GPG%NSURF_IDX
   SUMVAL = 0D0
   DO J = 1, GPROP%NGSPEC   
      SUMVAL = SUMVAL + GPG%ALLBC(I)%YJINF(J)
   ENDDO
   IF (SUMVAL .LT. 0.999999 .OR. SUMVAL .GT. 1.000001) THEN
      WRITE(THREE,'(I3.3)') I
      MESSAGE='Problem with gaseous species boundary condition # ' // THREE // '. Mass fractions do not sum to 1.0'
      CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF
ENDDO




GPG%NEED_GAS_YJ = GPG%SOLVE_GAS_YJ  ! If the user specified to solve gas YJ, then it is of course needed

DO IRXN = 1, SPROP%NRXN
   DO J = 1, GPROP%NGSPEC
      IF (GPROP%YIELDS(J, IRXN) .LT. 0D0) THEN
         ! That means a gaseous species is consumed, so the gas mass fraction is needed
         GPG%NEED_GAS_YJ = .TRUE.
      END IF
   ENDDO
ENDDO


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

DO IRXN = 1, GPROP%NHGRXN
   SUMVAL = 0D0
   DO J = 1, GPROP%NGSPEC
      SUMVAL = SUMVAL + GPROP%HGYIELDS(J,IRXN)
   ENDDO
   IF (SUMVAL .LT. -0.000001 .OR. SUMVAL .GT. 0.000001) THEN
      WRITE(TWO,'(I2.2)') IRXN
      MESSAGE='Problem with gaseous yields for homogeneous gaseous reaction ' // TWO // '. Yields do not sum to 1.0'
      CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF
ENDDO


! Check initial gas-phase species mass fractions:
IF (GPG%NEED_GAS_YJ) THEN
   DO ICINDEX = 1, GPG%NIC
      SUMVAL = 0D0
      DO J = 1, GPROP%NGSPEC
         SUMVAL = SUMVAL + GPG%INITIAL_CONDITIONS(ICINDEX)%YJ0(J)
      ENDDO
      IF (SUMVAL .LT. 0.999999 .OR. SUMVAL .GT. 1.000001) THEN
         WRITE(THREE,'(I3.3)') ICINDEX
         MESSAGE='Initial gas-phase mass fractions do not sum to 1.0 for ICNUM ' // THREE // &
               '. Check IC worksheet and make sure you have specified initial conditions for all gas-phase species.'
         CALL SHUTDOWN_GPYRO(MESSAGE)
      ENDIF
   ENDDO
ENDIF

IF (GPG%FDSMODE) THEN
   IF (GPG%FDS_MATL_VER .NE. 5 .AND. GPG%FDS_MATL_VER .NE. 6) THEN
      MESSAGE='Error, when FDSMODE=.TRUE., set FDS_MATL_VER to either 5 or 6'
      CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF   
ENDIF

! Warnings
IF (IRANK .EQ. 0 .AND. (.NOT. GPG%SHYI_CORRECTION)) THEN
   DO I = 1, SPROP%NSSPEC
      IF (SPROP%C0(I) .NE. SPROP%C0(1) .OR. SPROP%NC(I) .NE. SPROP%NC(1) ) THEN
         WRITE(*,*) '***** WARNING ***** '
         WRITE(*,*) 'If different solid species have different specific heat capacities '
         WRITE(*,*) 'Then Gpyro should be run with GPG%SHYI_CORRECTION = .TRUE. '
         WRITE(*,*)
      ENDIF   
   ENDDO
ENDIF

IF (IRANK .EQ. 0 .AND. GPG%EXPLICIT_T) THEN
   WRITE(*,*) '***** WARNING ***** '
   WRITE(*,*) 'Gpyro is an implicit code.' 
   WRITE(*,*) 'Setting EXPLICIT_T = .TRUE. should be used with caution.'
   WRITE(*,*)
ENDIF

IF (IRANK .EQ. 0 .AND. GPG%FDSMODE) THEN
   WRITE(*,*) '***** WARNING ***** '
   WRITE(*,*) 'Gpyro is running in FDS mode.' 
   WRITE(*,*) "Be sure to check Gpyro's .out file as some defaults are different in FDS mode."
   WRITE(*,*)
ENDIF

!******************************************************************************
END SUBROUTINE CHECK_GPYRO
!******************************************************************************	

!******************************************************************************	
SUBROUTINE DEALLOCATE_GPYRO
!******************************************************************************	

IF (ALLOCATED(GPM   )) DEALLOCATE(GPM   )

!******************************************************************************	
END SUBROUTINE DEALLOCATE_GPYRO
!******************************************************************************	

!******************************************************************************	
SUBROUTINE ALLOCATE_GPYRO(IMESH)
!******************************************************************************	

INTEGER, INTENT(IN) :: IMESH
INTEGER :: NGSPEC,NSSPEC,NRXN,NHGRXN
INTEGER :: IZ, IX, IY, IOR ,I, STR_SIZE

INTEGER :: NCELLZ 
INTEGER :: NCELLX 
INTEGER :: NCELLY

IF (.NOT. ALLOCATED (GPM)) ALLOCATE(GPM(GPG%NUM_GPYRO_MESHES))

! Initialize system clock
IF (IGPYRO_TYPE .EQ. 2) THEN !If this is an FDS run initialize the system clock here
   CALL SYSTEM_CLOCK(COUNT_RATE=CLOCK_COUNT_RATE)
   CALL SYSTEM_CLOCK(COUNT_MAX=CLOCK_COUNT_MAX)
ENDIF

G => GPM(IMESH)

! This is ALLOCATABLE (not POINTER) so if it's already allocated then we've already done
! this mesh so return
IF (ALLOCATED(G%CONV_INFO%CONVERGED_YIS)) RETURN
G%NCELLZ = GPG%NCELLZ(IMESH)
G%NCELLX = GPG%NCELLX(IMESH)
G%NCELLY = GPG%NCELLY(IMESH)
G%HALF_CELLS_AT_BC = GPG%HALF_CELLS_AT_BC(IMESH)
IF (.NOT. G%HALF_CELLS_AT_BC ) THEN
   ! If the cell at the boundary is a complete cell, add 2 ghost cells at the edge.
   IF (G%NCELLZ .NE. 1) G%NCELLZ = G%NCELLZ +2
   IF (G%NCELLX .NE. 1) G%NCELLX = G%NCELLX +2
   IF (G%NCELLY .NE. 1) G%NCELLY = G%NCELLY +2
ENDIF

G%ZDIM = GPG%ZDIM(IMESH); IF (G%NCELLZ .EQ. 1) G%ZDIM = 1.
G%XDIM = GPG%XDIM(IMESH); IF (G%NCELLX .EQ. 1) G%XDIM = 1.
G%YDIM = GPG%YDIM(IMESH); IF (G%NCELLY .EQ. 1) G%YDIM = 1.

!Allocate NGPYRO_FACES_NEEDING_BCS
IF (.NOT. ALLOCATED (NGPYRO_FACES_NEEDING_BCS) ) THEN 
   ALLOCATE(NGPYRO_FACES_NEEDING_BCS(0:GPG%NUM_GPYRO_MESHES))
   NGPYRO_FACES_NEEDING_BCS(:) = 0
ENDIF

!Set default gravity and position (then update for FDS)
G%GX = GPG%GX
G%GY = GPG%GY
G%GZ = GPG%GZ

IF (IGPYRO_TYPE .EQ. 2) THEN !FDS
   G%MESH_CENTROIDS(1:3) = GPG%MESH_CENTROIDS (IMESH,1:3)
   G%MESH_EXTENTS  (1:3) = GPG%MESH_EXTENTS   (IMESH,1:3)

   G%MESH_LL(1:3) = G%MESH_CENTROIDS(1:3) - 0.5D0 * G%MESH_EXTENTS(1:3)
   
   !Set x0, y0, z0:
   G%X0 = G%MESH_LL(1)
   G%Y0 = G%MESH_LL(2)
   G%Z0 = G%MESH_LL(3)
ELSE
   !Set x0, y0, z0:
   G%X0 = 0D0
   G%Y0 = 0D0
   G%Z0 = 0D0
ENDIF

NCELLZ = G%NCELLZ
NCELLX = G%NCELLX
NCELLY = G%NCELLY

NGSPEC = GPROP%NGSPEC
NSSPEC = SPROP%NSSPEC
NRXN   = SPROP%NRXN
NHGRXN = GPROP%NHGRXN

G%DIMENSION = 0
IF (NCELLZ .GT.1 ) G%DIMENSION = G%DIMENSION +1
IF (NCELLX .GT.1 ) G%DIMENSION = G%DIMENSION +1
IF (NCELLY .GT.1 ) G%DIMENSION = G%DIMENSION +1

!Allow deformation only in 1D
IF (G%DIMENSION .LE. 1) G%DEFORMATION = .TRUE.
IF (G%DIMENSION .GT. 1) G%DEFORMATION = .FALSE.

!For extension to 3D, consolidate all directional (x,y,z) quantities here:

!Coordinates in z, x, and y directions:
ALLOCATE (G%Z(1:NCELLZ) ); G%Z=0D0
ALLOCATE (G%X(1:NCELLX) ); G%X=0D0
ALLOCATE (G%Y(1:NCELLY) ); G%Y=0D0

ALLOCATE (G%FT(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%FT=0D0
ALLOCATE (G%FB(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%FB=0D0
ALLOCATE (G%FE(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%FE=0D0
ALLOCATE (G%FW(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%FW=0D0
ALLOCATE (G%FN(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%FN=0D0
ALLOCATE (G%FS(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%FS=0D0

ALLOCATE(G%NEEDSBCT(1:NCELLZ,1:NCELLX,1:NCELLY)); G%NEEDSBCT(:,:,:) = .FALSE.
ALLOCATE(G%NEEDSBCB(1:NCELLZ,1:NCELLX,1:NCELLY)); G%NEEDSBCB(:,:,:) = .FALSE.
ALLOCATE(G%NEEDSBCW(1:NCELLZ,1:NCELLX,1:NCELLY)); G%NEEDSBCW(:,:,:) = .FALSE.
ALLOCATE(G%NEEDSBCE(1:NCELLZ,1:NCELLX,1:NCELLY)); G%NEEDSBCE(:,:,:) = .FALSE.
ALLOCATE(G%NEEDSBCS(1:NCELLZ,1:NCELLX,1:NCELLY)); G%NEEDSBCS(:,:,:) = .FALSE.
ALLOCATE(G%NEEDSBCN(1:NCELLZ,1:NCELLX,1:NCELLY)); G%NEEDSBCN(:,:,:) = .FALSE.

ALLOCATE(G%SURF_IDX_BCT(1:NCELLZ,1:NCELLX,1:NCELLY)); G%SURF_IDX_BCT(:,:,:) = -1
ALLOCATE(G%SURF_IDX_BCB(1:NCELLZ,1:NCELLX,1:NCELLY)); G%SURF_IDX_BCB(:,:,:) = -1
ALLOCATE(G%SURF_IDX_BCW(1:NCELLZ,1:NCELLX,1:NCELLY)); G%SURF_IDX_BCW(:,:,:) = -1
ALLOCATE(G%SURF_IDX_BCE(1:NCELLZ,1:NCELLX,1:NCELLY)); G%SURF_IDX_BCE(:,:,:) = -1
ALLOCATE(G%SURF_IDX_BCS(1:NCELLZ,1:NCELLX,1:NCELLY)); G%SURF_IDX_BCS(:,:,:) = -1
ALLOCATE(G%SURF_IDX_BCN(1:NCELLZ,1:NCELLX,1:NCELLY)); G%SURF_IDX_BCN(:,:,:) = -1


CALL ALLOCATE_MEMORY_FOR_FIELD(IMESH, NCELLZ, NCELLX, NCELLY, NSSPEC, NGSPEC)
CALL ASSIGN_FIELD_TO_MEMORY

!Grid size in z, x, and y directions (current and "next"):

ALLOCATE (G%DLTX(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%DLTX=1D0
ALLOCATE (G%DLTY(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%DLTY=1D0

!ALLOCATE (G%DLTXN(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%DLTXN=1D0
!ALLOCATE (G%DLTYN(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%DLTYN=1D0

ALLOCATE (G%DXDY  (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%DXDY  =1D0
ALLOCATE (G%DXDZ  (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%DXDZ  =1D0
ALLOCATE (G%DYDZ  (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%DYDZ  =1D0
ALLOCATE (G%DV(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%DV=1D0

!Distance between cell centers in z, x, and y directions
ALLOCATE (G%DZT  (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%DZT=9D9 !Delta-z (top)
ALLOCATE (G%DZB  (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%DZB=9D9 !Delta-z (bottom)
ALLOCATE (G%DXE  (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%DXE=9D9 !Delta-x (east)
ALLOCATE (G%DXW  (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%DXW=9D9 !Delta-x (west)
ALLOCATE (G%DYN  (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%DYN=9D9 !Delta-y (north)
ALLOCATE (G%DYS  (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%DYS=9D9 !Delta-y (south)

ALLOCATE (G%KZ    (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%KZ=0D0
ALLOCATE (G%KX    (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%KX=0D0
ALLOCATE (G%KY    (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%KY=0D0


!k/c in x, y, and z directions:
ALLOCATE (G%KOCT  (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%KOCT=0D0 !k/c (top)
ALLOCATE (G%KOCB  (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%KOCB=0D0 !k/c (bottom)
ALLOCATE (G%KOCE  (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%KOCE=0D0 !k/c (east)
ALLOCATE (G%KOCW  (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%KOCW=0D0 !k/c (west)
ALLOCATE (G%KOCN  (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%KOCN=0D0 !k/c (north)
ALLOCATE (G%KOCS  (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%KOCS=0D0 !k/c (south)

ALLOCATE (G%KAPPA(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%KAPPA=0D0 
ALLOCATE(G%RESIDUAL_TMP(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RESIDUAL_TMP = 0D0 

IF (GPG%SOLVE_PRESSURE) THEN 
   ALLOCATE(G%RESIDUAL_P(1:NCELLZ,1:NCELLX,1:NCELLY))
   G%RESIDUAL_P = 0D0 
ENDIF

IF (GPG%SOLVE_GAS_ENERGY) THEN 
   ALLOCATE(G%RESIDUAL_HG(1:NCELLZ,1:NCELLX,1:NCELLY))
   G%RESIDUAL_HG = 0D0 
ENDIF

IF (NSSPEC .GT. 0) THEN 
   ALLOCATE(G%RESIDUAL_YIS(1:NSSPEC,1:NCELLZ,1:NCELLX,1:NCELLY))
   G%RESIDUAL_YIS = 0D0 
ENDIF

IF (GPG%SOLVE_GAS_YJ) THEN 
   ALLOCATE(G%RESIDUAL_YJG(1:NGSPEC,1:NCELLZ,1:NCELLX,1:NCELLY))
   G%RESIDUAL_YJG = 0D0 
ENDIF

IF (GPG%SHYI_CORRECTION) THEN
   ALLOCATE (G%KOCHIT  (1:NSSPEC,1:NCELLZ,1:NCELLX,1:NCELLY) ); G%KOCHIT=0D0 !(k/c) * hi (top)
   ALLOCATE (G%KOCHIB  (1:NSSPEC,1:NCELLZ,1:NCELLX,1:NCELLY) ); G%KOCHIB=0D0 !(k/c) * hi (bottom) 
   ALLOCATE (G%KOCHIE  (1:NSSPEC,1:NCELLZ,1:NCELLX,1:NCELLY) ); G%KOCHIE=0D0 !(k/c) * hi (east)
   ALLOCATE (G%KOCHIW  (1:NSSPEC,1:NCELLZ,1:NCELLX,1:NCELLY) ); G%KOCHIW=0D0 !(k/c) * hi (west)
   ALLOCATE (G%KOCHIN  (1:NSSPEC,1:NCELLZ,1:NCELLX,1:NCELLY) ); G%KOCHIN=0D0 !(k/c) * hi (north)
   ALLOCATE (G%KOCHIS  (1:NSSPEC,1:NCELLZ,1:NCELLX,1:NCELLY) ); G%KOCHIS=0D0 !(k/c) * hi (south)

   ALLOCATE (G%SHYI     (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%SHYI     =0D0

ENDIF

IF (GPG%SOLVE_PRESSURE) THEN
   ALLOCATE (G%BUOYANCYZ(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%BUOYANCYZ=0D0
   ALLOCATE (G%BUOYANCYX(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%BUOYANCYX=0D0
   ALLOCATE (G%BUOYANCYY(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%BUOYANCYY=0D0

   ALLOCATE (G%PERMZ      (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%PERMZ=0D0
   ALLOCATE (G%PERMX      (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%PERMX=0D0
   ALLOCATE (G%PERMY      (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%PERMY=0D0

   ALLOCATE (G%PERMONUT   (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%PERMONUT=0D0
   ALLOCATE (G%PERMONUB   (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%PERMONUB=0D0
   ALLOCATE (G%PERMONUE   (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%PERMONUE=0D0
   ALLOCATE (G%PERMONUW   (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%PERMONUW=0D0
   ALLOCATE (G%PERMONUN   (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%PERMONUN=0D0
   ALLOCATE (G%PERMONUS   (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%PERMONUS=0D0

   ALLOCATE (G%MDOTPPDARCYT(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%MDOTPPDARCYT=0D0
   ALLOCATE (G%MDOTPPDARCYB(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%MDOTPPDARCYB=0D0
   ALLOCATE (G%MDOTPPDARCYE(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%MDOTPPDARCYE=0D0
   ALLOCATE (G%MDOTPPDARCYW(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%MDOTPPDARCYW=0D0
   ALLOCATE (G%MDOTPPDARCYN(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%MDOTPPDARCYN=0D0
   ALLOCATE (G%MDOTPPDARCYS(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%MDOTPPDARCYS=0D0
ENDIF

ALLOCATE (G%INDICE_OF_BC_GAS_OUT(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%INDICE_OF_BC_GAS_OUT=0D0

!z-direction mass flux, calculated from continuity:
ALLOCATE (G%MDOTPPZ(0:NGSPEC,1:NCELLZ,1:NCELLX,1:NCELLY) ); G%MDOTPPZ=0D0

IF (GPG%SOLVE_GAS_YJ .OR. GPG%SOLVE_PRESSURE) THEN
   ALLOCATE (G%PSIRGDT(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%PSIRGDT=0D0
   ALLOCATE (G%PSIRGDB(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%PSIRGDB=0D0
   ALLOCATE (G%PSIRGDE(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%PSIRGDE=0D0
   ALLOCATE (G%PSIRGDW(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%PSIRGDW=0D0
   ALLOCATE (G%PSIRGDN(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%PSIRGDN=0D0
   ALLOCATE (G%PSIRGDS(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%PSIRGDS=0D0

   ALLOCATE (G%RGNT(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%RGNT=0D0
   ALLOCATE (G%RGNB(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%RGNB=0D0
   ALLOCATE (G%RGNE(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%RGNE=0D0
   ALLOCATE (G%RGNW(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%RGNW=0D0
   ALLOCATE (G%RGNN(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%RGNN=0D0
   ALLOCATE (G%RGNS(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%RGNS=0D0
ENDIF

!End of directional quantities

!Now allocate major variables:

ALLOCATE (G%CPS   (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%CPS=1D0
ALLOCATE (G%SHP   (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%SHP=0D0
ALLOCATE (G%SHM   (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%SHM=0D0
ALLOCATE (G%SHNET (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%SHNET=0D0

ALLOCATE (G%OMEGASFG(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%OMEGASFG=0D0
   
ALLOCATE (G%OMSOLIDFRAC(1:NRXN,1:NCELLZ,1:NCELLX,1:NCELLY) ); G%OMSOLIDFRAC=0D0

ALLOCATE (G%HI      (1:NSSPEC,1:NCELLZ,1:NCELLX,1:NCELLY) ); G%HI=0D0
ALLOCATE (G%RYIDZ0  (0:NSSPEC,1:NCELLZ,1:NCELLX,1:NCELLY) ); G%RYIDZ0=0D0
   
ALLOCATE (G%UNREACTEDNESS(1:NSSPEC,1:NCELLZ,1:NCELLX,1:NCELLY) ); G%UNREACTEDNESS=1D0
            
ALLOCATE (G%OMEGASFBK  (1:NRXN,1:NCELLZ,1:NCELLX,1:NCELLY) ); G%OMEGASFBK=0D0
ALLOCATE (G%OMEGASDAK  (1:NRXN,1:NCELLZ,1:NCELLX,1:NCELLY) ); G%OMEGASDAK=0D0
ALLOCATE (G%OMEGASFGK  (1:NRXN,1:NCELLZ,1:NCELLX,1:NCELLY) ); G%OMEGASFGK=0D0
            
ALLOCATE (G%RRT(1:NRXN,1:NCELLZ,1:NCELLX,1:NCELLY) ); G%RRT=0D0
ALLOCATE (G%RRY(1:NRXN,1:NCELLZ,1:NCELLX,1:NCELLY) ); G%RRY=0D0

ALLOCATE (G%SOMEGA   (1:3,1:NSSPEC,1:NCELLZ,1:NCELLX,1:NCELLY)); G%SOMEGA=0D0
ALLOCATE (G%GOMEGA   (1:3,0:NGSPEC,1:NCELLZ,1:NCELLX,1:NCELLY)); G%GOMEGA=0D0
ALLOCATE (G%OMEGASFJK(1:NGSPEC,1:NRXN,1:NCELLZ,1:NCELLX,1:NCELLY)); G%OMEGASFJK=0D0
ALLOCATE (G%OMEGASDJK(1:NGSPEC,1:NRXN,1:NCELLZ,1:NCELLX,1:NCELLY)); G%OMEGASDJK=0D0

ALLOCATE (G%IS_REACTING(0:NRXN,1:NCELLZ,1:NCELLX,1:NCELLY));G%IS_REACTING=.FALSE.


ALLOCATE (G%OMEGAGDJL(1:NGSPEC,1:NHGRXN,1:NCELLZ,1:NCELLX,1:NCELLY)); G%OMEGAGDJL=0D0
                        
ALLOCATE (G%QSG       (1:NCELLZ,1:NCELLX,1:NCELLY));G%QSG=0D0
ALLOCATE (G%CONSUMED  (1:NCELLZ,1:NCELLX,1:NCELLY));G%CONSUMED=.FALSE.

ALLOCATE(G%IMASK(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%IMASK= .FALSE.
            

IF (GPG%NEED_GAS_YJ) THEN               
   ALLOCATE (G%YJG     (1:NGSPEC,1:NCELLZ,1:NCELLX,1:NCELLY) ); G%YJG=0D0
   ALLOCATE (G%YJGN    (1:NGSPEC,1:NCELLZ,1:NCELLX,1:NCELLY) ); G%YJGN=0D0
ENDIF

IF (GPG%SOLVE_GAS_ENERGY) THEN
   ALLOCATE (G%SHGP   (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%SHGP=0D0
   ALLOCATE (G%SHGM   (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%SHGM=0D0
   ALLOCATE (G%HGRR   (1:NHGRXN,1:NCELLZ,1:NCELLX,1:NCELLY)); G%HGRR=0D0
   IF (GPROP%NHGRXN .GT. 0) THEN
      ALLOCATE (G%HGOMEGA(1:3,0:NGSPEC,1:NCELLZ,1:NCELLX,1:NCELLY))
      G%HGOMEGA=0D0
   ENDIF
ENDIF

            
IF (GPG%SOLVE_PRESSURE .OR. GPG%SOLVE_GAS_YJ) THEN
   ALLOCATE (G%D12(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%D12=0D0
ENDIF
            

IF ( (.NOT. GPG%THERMAL_EQUILIBRIUM) .AND. GPG%HCV .LT. 0D0) THEN 
   ALLOCATE(G%HCV(1:NCELLZ,1:NCELLX,1:NCELLY) ); G%HCV=0D0
   ALLOCATE(G%RE (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%RE=0D0
   ALLOCATE(G%NU (1:NCELLZ,1:NCELLX,1:NCELLY) ); G%NU=0D0
ENDIF

ALLOCATE (G%THICKNESS (1:NCELLX,1:NCELLY) ); G%THICKNESS = 0D0

! If this is a TGA experiment...
IF (NCELLZ .EQ. 1 .AND. (.NOT. GPG%NEED_GAS_YJ) ) THEN
   ALLOCATE (G%YJG     (1:NGSPEC,1:NCELLZ,1:NCELLX,1:NCELLY) ); G%YJG=0D0
   ALLOCATE (G%YJGN    (1:NGSPEC,1:NCELLZ,1:NCELLX,1:NCELLY) ); G%YJGN=0D0
ENDIF

!Zero arrays for applying boundary conditions:
ALLOCATE(G%GPYRO_BOUNDARY_CONDITION(1:NCELLZ,1:NCELLX,1:NCELLY,-3:3))

GPBCP(1:,1:,1:,-3:)=>G%GPYRO_BOUNDARY_CONDITION(1:NCELLZ,1:NCELLX,1:NCELLY,-3:3)

GPBCP(:,:,:,:)%QE           = 0D0
GPBCP(:,:,:,:)%QENET        = 0D0
GPBCP(:,:,:,:)%HC0          = 12D0
GPBCP(:,:,:,:)%NHC          = 0D0
GPBCP(:,:,:,:)%TINF         = GPG%TAMB
GPBCP(:,:,:,:)%TFIXED       = -1000D0
GPBCP(:,:,:,:)%HFIXED       = 0D0
GPBCP(:,:,:,:)%PRES         = GPG%P0
GPBCP(:,:,:,:)%MFLUX        = 0D0
GPBCP(:,:,:,:)%HM0          = 0.012D0
GPBCP(:,:,:,:)%QEG          = 0D0
GPBCP(:,:,:,:)%HC0G         = 12D0
GPBCP(:,:,:,:)%TINFG        = GPG%TAMB
GPBCP(:,:,:,:)%TFIXEDG      = -1000D0
GPBCP(:,:,:,:)%HFIXEDG      = 0D0
GPBCP(:,:,:,:)%RERAD        = .FALSE.
GPBCP(:,:,:,:)%EMISSIVITY   = 0D0
GPBCP(:,:,:,:)%T_SURFACE    = GPG%TAMB
GPBCP(:,:,:,:)%T_SURFACE_OLD= GPG%TAMB
GPBCP(:,:,:,:)%QRADOUT      = 0D0
GPBCP(:,:,:,:)%QCONF        = 0D0

DO IZ = 1, NCELLZ
DO IX = 1, NCELLX
DO IY = 1, NCELLY
DO IOR =-3, 3
   ALLOCATE(G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%YJINF   (1:GPROP%NGSPEC))
   G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%YJINF   (1:GPROP%NGSPEC) = 0D0
   ALLOCATE(G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%MASSFLUX(1:GPROP%NGSPEC))
   G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%MASSFLUX(1:GPROP%NGSPEC) = 0D0
   ALLOCATE(G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%MLR     (1:GPROP%NGSPEC))
   G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%MLR     (1:GPROP%NGSPEC) = 0D0
ENDDO
ENDDO
ENDDO
ENDDO

ALLOCATE(G%RWORK01(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK01=0D0
ALLOCATE(G%RWORK02(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK02=0D0
ALLOCATE(G%RWORK03(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK03=0D0
ALLOCATE(G%RWORK04(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK04=0D0
ALLOCATE(G%RWORK05(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK05=0D0
ALLOCATE(G%RWORK06(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK06=0D0
ALLOCATE(G%RWORK07(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK07=0D0
ALLOCATE(G%RWORK08(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK08=0D0
ALLOCATE(G%RWORK09(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK09=0D0
ALLOCATE(G%RWORK10(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK10=0D0
ALLOCATE(G%RWORK11(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK11=0D0
ALLOCATE(G%RWORK12(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK12=0D0
ALLOCATE(G%RWORK13(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK13=0D0
ALLOCATE(G%RWORK14(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK14=0D0
ALLOCATE(G%RWORK15(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK15=0D0
ALLOCATE(G%RWORK16(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK16=0D0
ALLOCATE(G%RWORK17(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK17=0D0
ALLOCATE(G%RWORK18(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK18=0D0
ALLOCATE(G%RWORK19(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK19=0D0
ALLOCATE(G%RWORK20(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK20=0D0
ALLOCATE(G%RWORK21(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK21=0D0
ALLOCATE(G%RWORK22(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK22=0D0
ALLOCATE(G%RWORK23(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK23=0D0
ALLOCATE(G%RWORK24(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK24=0D0
ALLOCATE(G%RWORK25(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK25=0D0
ALLOCATE(G%RWORK26(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK26=0D0
ALLOCATE(G%RWORK27(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK27=0D0
ALLOCATE(G%RWORK28(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK28=0D0
ALLOCATE(G%RWORK29(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK29=0D0
ALLOCATE(G%RWORK30(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK30=0D0
ALLOCATE(G%RWORK31(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK31=0D0
ALLOCATE(G%RWORK32(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK32=0D0
ALLOCATE(G%RWORK33(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK33=0D0
ALLOCATE(G%RWORK34(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK34=0D0
ALLOCATE(G%RWORK35(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK35=0D0
ALLOCATE(G%RWORK36(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK36=0D0
ALLOCATE(G%RWORK37(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK37=0D0
ALLOCATE(G%RWORK38(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK38=0D0
ALLOCATE(G%RWORK39(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK39=0D0
ALLOCATE(G%RWORK40(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK40=0D0
ALLOCATE(G%RWORK41(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK41=0D0
ALLOCATE(G%RWORK42(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK42=0D0
ALLOCATE(G%RWORK43(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK43=0D0
ALLOCATE(G%RWORK44(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK44=0D0
ALLOCATE(G%RWORK45(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK45=0D0
ALLOCATE(G%RWORK46(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK46=0D0
ALLOCATE(G%RWORK47(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK47=0D0
ALLOCATE(G%RWORK48(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK48=0D0
ALLOCATE(G%RWORK49(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK49=0D0
ALLOCATE(G%RWORK50(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK50=0D0
ALLOCATE(G%RWORK51(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK51=0D0
ALLOCATE(G%RWORK52(1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK52=0D0

ALLOCATE(G%RWORK100(1:NSSPEC,1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK100=0D0
ALLOCATE(G%RWORK101(1:NSSPEC,1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK101=0D0
ALLOCATE(G%RWORK102(1:NSSPEC,1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK102=0D0
ALLOCATE(G%RWORK103(1:NSSPEC,1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK103=0D0
ALLOCATE(G%RWORK104(1:NSSPEC,1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK104=0D0
ALLOCATE(G%RWORK105(1:NSSPEC,1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK105=0D0
ALLOCATE(G%RWORK106(1:NSSPEC,1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK106=0D0

ALLOCATE(G%RWORK110(1:NGSPEC,1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK110=0D0
ALLOCATE(G%RWORK111(1:NGSPEC,1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK111=0D0
ALLOCATE(G%RWORK112(1:NGSPEC,1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK112=0D0

ALLOCATE(G%RWORK120(1:SPROP%NRXN,1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK120=0D0
ALLOCATE(G%RWORK121(1:SPROP%NRXN,1:NCELLZ,1:NCELLX,1:NCELLY)); G%RWORK121=0D0

ALLOCATE(G%LWORK01(1:NCELLZ,1:NCELLX,1:NCELLY)); G%LWORK01=.FALSE.

ALLOCATE(G%WORK_DUMP(GPG%N_PROFILE_QUANTITIES))
DO I = 1, GPG%N_PROFILE_QUANTITIES
   STR_SIZE= MAX(G%NCELLX,G%NCELLY,G%NCELLZ)+1
   ALLOCATE(G%WORK_DUMP(I)%BUFFER(STR_SIZE))
   ALLOCATE(CHARACTER(STR_SIZE * 20) :: G%WORK_DUMP(I)%LINE)
ENDDO


ALLOCATE(G%CONV_INFO%CONVERGED_YIS(0:NSSPEC))
ALLOCATE(G%CONV_INFO%CONVERGED_YJG(0:NGSPEC))
ALLOCATE(G%CONV_INFO%ITER_YIS(0:NSSPEC))
ALLOCATE(G%CONV_INFO%ITER_YJG(0:NGSPEC))

!******************************************************************************	
END SUBROUTINE ALLOCATE_GPYRO
!******************************************************************************	

!******************************************************************************	
SUBROUTINE ALLOCATE_GPYRO_BOUNDARYS_INFO(NMESHES,IMESH,IDIM)
!******************************************************************************	

INTEGER, INTENT(IN) :: NMESHES, IMESH, IDIM
INTEGER :: ICOUNT !, LUSIZE

IF (.NOT. ALLOCATED(GP_BOUDARYS)) ALLOCATE (GP_BOUDARYS(0:NMESHES))

ALLOCATE (GP_BOUDARYS(IMESH)%IMESH_GPYRO(1:IDIM)); GP_BOUDARYS(IMESH)%IMESH_GPYRO(:) = 0
ALLOCATE (GP_BOUDARYS(IMESH)%IZ_GPYRO   (1:IDIM)); GP_BOUDARYS(IMESH)%IZ_GPYRO   (:) = 0
ALLOCATE (GP_BOUDARYS(IMESH)%IX_GPYRO   (1:IDIM)); GP_BOUDARYS(IMESH)%IX_GPYRO   (:) = 0
ALLOCATE (GP_BOUDARYS(IMESH)%IY_GPYRO   (1:IDIM)); GP_BOUDARYS(IMESH)%IY_GPYRO   (:) = 0
ALLOCATE (GP_BOUDARYS(IMESH)%IOR_GPYRO  (1:IDIM)); GP_BOUDARYS(IMESH)%IOR_GPYRO  (:) = 0
ALLOCATE (GP_BOUDARYS(IMESH)%IOR_FDS    (1:IDIM)); GP_BOUDARYS(IMESH)%IOR_FDS    (:) = 0
ALLOCATE (GP_BOUDARYS(IMESH)%IMESH_FDS  (1:IDIM)); GP_BOUDARYS(IMESH)%IMESH_FDS  (:) = 0
ALLOCATE (GP_BOUDARYS(IMESH)%I_FDS      (1:IDIM)); GP_BOUDARYS(IMESH)%I_FDS      (:) = 0
ALLOCATE (GP_BOUDARYS(IMESH)%J_FDS      (1:IDIM)); GP_BOUDARYS(IMESH)%J_FDS      (:) = 0
ALLOCATE (GP_BOUDARYS(IMESH)%K_FDS      (1:IDIM)); GP_BOUDARYS(IMESH)%K_FDS      (:) = 0
ALLOCATE (GP_BOUDARYS(IMESH)%IW_FDS     (1:IDIM)); GP_BOUDARYS(IMESH)%IW_FDS     (:) = 0
ALLOCATE (GP_BOUDARYS(IMESH)%RATIO      (1:IDIM)); GP_BOUDARYS(IMESH)%RATIO      (:) = 0
ALLOCATE (GP_BOUDARYS(IMESH)%DIFF       (1:IDIM)); GP_BOUDARYS(IMESH)%DIFF       (:) = 9D9
ALLOCATE (GP_BOUDARYS(IMESH)%XDIFF      (1:IDIM)); GP_BOUDARYS(IMESH)%XDIFF      (:) = 9D9
ALLOCATE (GP_BOUDARYS(IMESH)%YDIFF      (1:IDIM)); GP_BOUDARYS(IMESH)%YDIFF      (:) = 9D9
ALLOCATE (GP_BOUDARYS(IMESH)%ZDIFF      (1:IDIM)); GP_BOUDARYS(IMESH)%ZDIFF      (:) = 9D9
ALLOCATE (GP_BOUDARYS(IMESH)%Z_GPYRO    (1:IDIM)); GP_BOUDARYS(IMESH)%Z_GPYRO    (:) = 0D0
ALLOCATE (GP_BOUDARYS(IMESH)%X_GPYRO    (1:IDIM)); GP_BOUDARYS(IMESH)%X_GPYRO    (:) = 0D0
ALLOCATE (GP_BOUDARYS(IMESH)%Y_GPYRO    (1:IDIM)); GP_BOUDARYS(IMESH)%Y_GPYRO    (:) = 0D0
ALLOCATE (GP_BOUDARYS(IMESH)%X_FDS      (1:IDIM)); GP_BOUDARYS(IMESH)%X_FDS      (:) = 0D0
ALLOCATE (GP_BOUDARYS(IMESH)%Y_FDS      (1:IDIM)); GP_BOUDARYS(IMESH)%Y_FDS      (:) = 0D0
ALLOCATE (GP_BOUDARYS(IMESH)%Z_FDS      (1:IDIM)); GP_BOUDARYS(IMESH)%Z_FDS      (:) = 0D0
ALLOCATE (GP_BOUDARYS(IMESH)%COMPLETE_CELL_AT_BC(1:IDIM)); GP_BOUDARYS(IMESH)%COMPLETE_CELL_AT_BC(:) = .FALSE.
ALLOCATE (GP_BOUDARYS(IMESH)%GPBC       (1:IDIM))

DO ICOUNT = 1, IDIM
   GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%QE            = 0D0
   GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%QENET         = 0D0
   GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%HC0           = 0D0
   GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%TINF          = 0D0
   GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%TFIXED        = 0D0
   GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%HFIXED        = 0D0
   GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%PRES          = 0D0
   GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%MFLUX         = 0D0
   GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%HM0           = 0D0
   GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%QEG           = 0D0
   GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%HC0G          = 0D0
   GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%TINFG         = 0D0
   GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%TFIXEDG       = 0D0
   GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%HFIXEDG       = 0D0
   GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%RERAD         = .FALSE.
   GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%EMISSIVITY    = 0D0
   GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%T_SURFACE     = 0D0
   GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%T_SURFACE_OLD = 0D0
   GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%QRADOUT       = 0D0
   GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%QCONF         = 0D0
   ALLOCATE(GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%YJINF   (1:GPROP%NGSPEC)); GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%YJINF   (:) = 0D0
   ALLOCATE(GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%MASSFLUX(1:GPROP%NGSPEC)); GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%MASSFLUX(:) = 0D0
   ALLOCATE(GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%MLR     (1:GPROP%NGSPEC)); GP_BOUDARYS(IMESH)%GPBC(ICOUNT)%MLR     (:) = 0D0
ENDDO

!******************************************************************************	
END SUBROUTINE ALLOCATE_GPYRO_BOUNDARYS_INFO
!******************************************************************************	


! *****************************************************************************
SUBROUTINE ALLOCATE_MEMORY_FOR_FIELD(IMESH, NCELLZ, NCELLX, NCELLY, NSSPEC, NGSPEC)
! *****************************************************************************
INTEGER, INTENT(IN) :: IMESH, NCELLZ, NCELLX, NCELLY
INTEGER, INTENT(IN) :: NSSPEC, NGSPEC
TYPE(GPYRO_MESH_TYPE), POINTER :: G
TYPE(GPYRO_STORAGE_TYPE), POINTER :: S

G => GPM(IMESH)
S => G%STORAGE

! Arrays 3D (NCELLZ, NCELLX, NCELLY)
ALLOCATE(S%TP_W1(NCELLZ,NCELLX,NCELLY),     S%TP_W2(NCELLZ,NCELLX,NCELLY))
ALLOCATE(S%HP_W1(NCELLZ,NCELLX,NCELLY),     S%HP_W2(NCELLZ,NCELLX,NCELLY))
ALLOCATE(S%RP_W1(NCELLZ,NCELLX,NCELLY),     S%RP_W2(NCELLZ,NCELLX,NCELLY))
ALLOCATE(S%DLTZ_W1(NCELLZ,NCELLX,NCELLY),   S%DLTZ_W2(NCELLZ,NCELLX,NCELLY))
ALLOCATE(S%RDLTZ_W1(NCELLZ,NCELLX,NCELLY),  S%RDLTZ_W2(NCELLZ,NCELLX,NCELLY))

IF (GPG%SOLVE_GAS_YJ .OR. GPG%SOLVE_PRESSURE) THEN
   ALLOCATE(S%RG_W1(NCELLZ,NCELLX,NCELLY),     S%RG_W2(NCELLZ,NCELLX,NCELLY))
ENDIF
IF (GPG%SOLVE_PRESSURE) THEN
   ALLOCATE(S%P_W1(NCELLZ,NCELLX,NCELLY),      S%P_W2(NCELLZ,NCELLX,NCELLY))
   ALLOCATE(S%M_W1(NCELLZ,NCELLX,NCELLY),      S%M_W2(NCELLZ,NCELLX,NCELLY))
ENDIF

IF (GPG%SOLVE_GAS_ENERGY) THEN
   ALLOCATE(S%TG_W1(NCELLZ,NCELLX,NCELLY),     S%TG_W2(NCELLZ,NCELLX,NCELLY))
   ALLOCATE(S%HG_W1(NCELLZ,NCELLX,NCELLY),     S%HG_W2(NCELLZ,NCELLX,NCELLY))
ENDIF

IF (GPG%SOLVE_POROSITY) THEN
   ALLOCATE(S%RSP_W1(NCELLZ,NCELLX,NCELLY),    S%RSP_W2(NCELLZ,NCELLX,NCELLY))
   ALLOCATE(S%POROSS_W1(NCELLZ,NCELLX,NCELLY), S%POROSS_W2(NCELLZ,NCELLX,NCELLY))
ENDIF

! Arrays 4D (NSSPEC, NCELLZ, NCELLX, NCELLY)
ALLOCATE(S%YI_W1(NSSPEC,NCELLZ,NCELLX,NCELLY),       S%YI_W2(NSSPEC,NCELLZ,NCELLX,NCELLY))
ALLOCATE(S%XI_W1(NSSPEC,NCELLZ,NCELLX,NCELLY),       S%XI_W2(NSSPEC,NCELLZ,NCELLX,NCELLY))
ALLOCATE(S%RYIDZP_W1(NSSPEC,NCELLZ,NCELLX,NCELLY),   S%RYIDZP_W2(NSSPEC,NCELLZ,NCELLX,NCELLY))
ALLOCATE(S%RYIDZSIGMA_W1(NSSPEC,NCELLZ,NCELLX,NCELLY), S%RYIDZSIGMA_W2(NSSPEC,NCELLZ,NCELLX,NCELLY))

! Arrays 4D (NGSPEC, NCELLZ, NCELLX, NCELLY)
IF (GPG%SOLVE_GAS_YJ) THEN
   ALLOCATE(S%YJG_W1(NGSPEC,NCELLZ,NCELLX,NCELLY),      S%YJG_W2(NGSPEC,NCELLZ,NCELLX,NCELLY))
ENDIF
! *****************************************************************************
END SUBROUTINE ALLOCATE_MEMORY_FOR_FIELD
! *****************************************************************************


! *****************************************************************************
SUBROUTINE ASSIGN_FIELD_TO_MEMORY
! *****************************************************************************

! Assign field pointers to their working memory (W1/W2)
G%TP    => G%STORAGE%TP_W1
G%TPN   => G%STORAGE%TP_W2

G%HP    => G%STORAGE%HP_W1
G%HPN   => G%STORAGE%HP_W2

G%RP    => G%STORAGE%RP_W1
G%RPN   => G%STORAGE%RP_W2

G%DLTZ  => G%STORAGE%DLTZ_W1
G%DLTZN => G%STORAGE%DLTZ_W2

G%RDLTZ  => G%STORAGE%RDLTZ_W1
G%RDLTZN => G%STORAGE%RDLTZ_W2


IF (GPG%SOLVE_GAS_YJ .OR. GPG%SOLVE_PRESSURE) THEN
   G%RG    => G%STORAGE%RG_W1
   G%RGN   => G%STORAGE%RG_W2
ENDIF

IF (GPG%SOLVE_PRESSURE) THEN
   G%P     => G%STORAGE%P_W1
   G%PN    => G%STORAGE%P_W2
   G%M     => G%STORAGE%M_W1
   G%MN    => G%STORAGE%M_W2
ENDIF

IF (GPG%SOLVE_GAS_ENERGY) THEN
   G%TG    => G%STORAGE%TG_W1
   G%TGN   => G%STORAGE%TG_W2
   G%HG    => G%STORAGE%HG_W1
   G%HGN   => G%STORAGE%HG_W2
ENDIF

IF (GPG%SOLVE_POROSITY) THEN
   G%RSP      => G%STORAGE%RSP_W1
   G%RSPN     => G%STORAGE%RSP_W2
   G%POROSS   => G%STORAGE%POROSS_W1
   G%POROSSN  => G%STORAGE%POROSS_W2
ENDIF

G%YI    => G%STORAGE%YI_W1
G%YIN   => G%STORAGE%YI_W2

G%XI    => G%STORAGE%XI_W1
G%XIN   => G%STORAGE%XI_W2

G%RYIDZP  => G%STORAGE%RYIDZP_W1
G%RYIDZPN => G%STORAGE%RYIDZP_W2

G%RYIDZSIGMA  => G%STORAGE%RYIDZSIGMA_W1
G%RYIDZSIGMAN => G%STORAGE%RYIDZSIGMA_W2

IF (GPG%SOLVE_GAS_YJ) THEN
   G%YJG   => G%STORAGE%YJG_W1
   G%YJGN  => G%STORAGE%YJG_W2
ENDIF
! *****************************************************************************
END SUBROUTINE ASSIGN_FIELD_TO_MEMORY
! *****************************************************************************



!******************************************************************************	
SUBROUTINE INIT_GPYRO(ICASE,IMESH)
!******************************************************************************	

INTEGER, INTENT(IN) :: ICASE,IMESH


! This routine is called once for every mesh, so the following happens 
! multiple times but it's extremely fast so there is no reason to 
! add logic to only do it once:


CALL INIT_GPYRO_GENERAL(IMESH)

CALL INIT_REACTIONS


!========================================================!
!================== BUILD MESH AND IC ===================!
!========================================================!

! Initialise G%X, G%Y, G%Z, G%DLTZ, G%DLZN, G%DLTX, G%DLTY
! G%DZT, G%DZB, G%DXE, G%DXW, G%DYN, G%DYS
! G%DXDY, G%DXDZ,G%DYDZ, G%DV
CALL GENERATE_RAW_GRID(IMESH)
!Set initial mass fractions, temperature, and pressure:
! IMASK, YIN, YI, TPN, TP,
! P, PN, YJG, YJGN, TG, TGN
! GPBCP%T_SURFACE, GPBCP%T_SURFACE_OLD
CALL APPLY_GEOMETRY_TO_GRID(ICASE,IMESH)
! GP_BOUDARYS
CALL INIT_BOUNDARY_CELLS(IMESH)
!========================================================!
!========================================================!

CALL INIT_GPYRO_VARS(IMESH)


CALL INIT_DUMP_OUTPOUT

IF (IGPYRO_TYPE .EQ. 1 .AND. IMESH .EQ. 1) CALL WRITE_DOTOUT_FILE(ICASE,0D0,1) 

!******************************************************************************	
END SUBROUTINE INIT_GPYRO
!******************************************************************************	



!******************************************************************************	
SUBROUTINE INIT_GPYRO_GENERAL(IMESH)
!******************************************************************************
INTEGER, INTENT(IN) :: IMESH

G=>GPM(IMESH)

! Alpha coefficients for relaxation:
IF (GPG%ALPHA .GT. 0D0 .AND. GPG%ALPHA .LE. 1D0) THEN
   GPG%ALPHA_YIS = GPG%ALPHA
   GPG%ALPHA_YJG = GPG%ALPHA
   GPG%ALPHA_HG  = GPG%ALPHA
   GPG%ALPHA_H   = GPG%ALPHA
   GPG%ALPHA_P   = GPG%ALPHA
ENDIF

GPG%TDATUM = 200D0

GPG%POSINF =  1. / TINY(0D0)
GPG%NEGINF = -1. / TINY(0D0)

GPG%TUSED(:) = 0D0 !Zero out CPU timing array

GPG%TDUMPLAST_GA           = 0D0
GPG%TDUMPLAST_POINT    (:) = 0D0
GPG%TDUMPLAST_PROFILE  (:) = 0D0
GPG%TDUMPLAST_SMOKEVIEW(:) = 0D0
GPG%TDUMPLAST_CSVDUMP  (:) = 0D0
GPG%DTNEXT                 = GPG%DT0

! Initialise time
G%NTIMESTEPS= 0
G%TLAST_TMP = 0D0
G%TLAST_YIS = 0D0
G%TLAST_YJG = 0D0
G%TLAST_P   = 0D0
G%TLAST_HG  = 0D0

!******************************************************************************
END SUBROUTINE INIT_GPYRO_GENERAL
!******************************************************************************



!******************************************************************************	
SUBROUTINE GENERATE_RAW_GRID(IMESH)
!******************************************************************************
! Initialise G%X, G%Y, G%Z, G%DLTZ, G%DLZN, G%DLTX, G%DLTY
! G%DZT, G%DZB, G%DXE, G%DXW, G%DYN, G%DYS
! G%DXDY, G%DXDZ,G%DYDZ, G%DV
INTEGER, INTENT(IN) :: IMESH
INTEGER :: NCELLZ, NCELLX, NCELLY
REAL(EB) :: ZDIM, XDIM, YDIM
INTEGER :: IZ,IX,IY


G=>GPM(IMESH)

NCELLZ = G%NCELLZ !Number of cells in z direction
NCELLX = G%NCELLX !Number of cells in x direction
NCELLY = G%NCELLY !Number of cells in y direction

ZDIM   = G%ZDIM
XDIM   = G%XDIM
YDIM   = G%YDIM 

G%THICKNESS   = ZDIM

!Set z-direction spacing:
IF ((NCELLZ .GT. 1) .AND. G%HALF_CELLS_AT_BC) THEN
   !! Classical mode implemented by CL for building the geometry with half-boundary cells at the edges. !!

   ! Set DLTZN and DZN:
   G%DLTZN (:     ,:,:) = ZDIM/REAL(NCELLZ-1,EB)
   G%DLTZN (1     ,:,:) = 0.5D0 * G%DLTZN(1     ,:,:) ! Half Boundary cell
   G%DLTZN (NCELLZ,:,:) = 0.5D0 * G%DLTZN(NCELLZ,:,:) ! Half Boundary cell
   G%DLTZ  (:,:,:) = G%DLTZN(:,:,:)

   !Delta-z top and bottom:
   G%DZT(:,:,:) = ZDIM/REAL(NCELLZ-1,EB)
   G%DZB(:,:,:) = G%DZT(:,:,:)

   ! z-position of the center of the cells
   G%Z(1) = 0D0
   G%Z(2) = G%Z(1) + G%DLTZN(1,1,1) + 0.5D0*G%DLTZN(2,1,1)
   DO IZ = 3, NCELLZ-1
      G%Z(IZ) = G%Z(IZ-1) + 0.5D0*(G%DLTZN(IZ-1,1,1) + G%DLTZN(IZ,1,1))
   ENDDO
   G%Z(NCELLZ) = G%Z(NCELLZ-1) + 0.5D0 * G%DLTZN(NCELLZ-1,1,1) + G%DLTZN(NCELLZ,1,1)    

ELSEIF ((NCELLZ .GT. 1) .AND. (.NOT. G%HALF_CELLS_AT_BC)) THEN
   ! New mode for building mesh with complete cell at the edge and with ghost cell.
   ! NZ=1 and NZ=NCELLZ is the ghost cell.

   G%DLTZN (:,:,:) = ZDIM/REAL(NCELLZ-2,EB)
   G%DLTZ  (:,:,:) = G%DLTZN(:,:,:)
   
   !Delta-z top and bottom:
   G%DZT(:,:,:) = ZDIM/REAL(NCELLZ-2,EB)
   G%DZB(:,:,:) = G%DZT(:,:,:)

   ! z-position of the center of the cells
   G%Z(1) = -0.5D0 *G%DLTZN (1,1,1) ! Center of the gost cell (not used)
   DO IZ = 2, NCELLZ
      G%Z(IZ) = G%Z(IZ-1) + 0.5D0*(G%DLTZN(IZ-1,1,1) + G%DLTZN(IZ,1,1))
   ENDDO

ELSEIF (NCELLZ .EQ. 1) THEN
   G%DLTZN(:,:,:) = 1D0
   G%DLTZ (:,:,:) = 1D0
   G%DZT(:,:,:) = 1D0
   G%DZB(:,:,:) = 1D0
   G%Z(1) = 0D0
ENDIF

!Set x-direction spacing:
IF ((NCELLX .GT. 1) .AND. G%HALF_CELLS_AT_BC) THEN
   G%DLTX (:,:,:) = XDIM/REAL(NCELLX-1,EB)
   G%DLTX (:,1,:) = 0.5D0 * G%DLTX(:,1,:)
   G%DLTX (:,NCELLX,:) = 0.5D0 * G%DLTX(:,NCELLX,:) 
   !G%DLTXN  (:,:,:) = G%DLTX(:,:,:)

   !Delta-x east and west:
   G%DXE(:,:,:) = XDIM/REAL(NCELLX-1,EB)
   G%DXW(:,:,:) = G%DXE(:,:,:) 

   ! Set x of each cell:
   G%X(1) = 0D0
   G%X(2) = G%X(1) + G%DLTX (1,1,1) + 0.5D0*G%DLTX (1,2,1)
   DO IX = 3, NCELLX-1
      G%X(IX) = G%X(IX-1) + 0.5D0*(G%DLTX (1,IX-1,1) + G%DLTX (1,IX,1))
   ENDDO
   G%X(NCELLX) = G%X(NCELLX-1) + 0.5D0 * G%DLTX (1,NCELLX-1,1) + G%DLTX (1,NCELLX,1)    

ELSEIF ((NCELLX .GT. 1) .AND. (.NOT. G%HALF_CELLS_AT_BC)) THEN

   G%DLTX (:,:,:) = XDIM/REAL(NCELLX-2,EB) 
   !G%DLTXN(:,:,:) = G%DLTX(:,:,:)

   !Delta-x east and west:
   G%DXE(:,:,:) = XDIM/REAL(NCELLX-2,EB)
   G%DXW(:,:,:) = G%DXE(:,:,:) 

   ! x-position of the center of the cells
   G%X(1) = -0.5D0 *G%DLTX (1,1,1) ! Center of the gost cell (not used)
   DO IX = 2, NCELLX
      G%X(IX) = G%X(IX-1) + 0.5D0*(G%DLTX(IX-1,1,1) + G%DLTX(IX,1,1))
   ENDDO

ELSE ! NCELLX=1
   G%DLTX (:,:,:) = 1D0   
   !G%DLTXN(:,:,:) = 1D0
ENDIF

!Set y-direction spacing:
IF ((NCELLY .GT. 1) .AND. G%HALF_CELLS_AT_BC)  THEN
   G%DLTY (:,:,:) = YDIM/REAL(NCELLY-1,EB)
   G%DLTY (:,:,1) = 0.5D0 * G%DLTY(:,:,1)
   G%DLTY(:,:,NCELLY) = 0.5D0 * G%DLTY(:,:,NCELLY) 
   !G%DLTYN (:,:,:) = G%DLTY(:,:,:)

   !Delta-y north and south:
   G%DYN(:,:,:) = YDIM/REAL(NCELLY-1,EB)
   G%DYS(:,:,:) = G%DYN(:,:,:) 

   ! Set y of each cell:
   G%Y(1) = 0D0
   G%Y(2) = G%Y(1) + G%DLTY(1,1,1) + 0.5D0*G%DLTY(1,1,2)
   DO IY = 3, NCELLY-1
      G%Y(IY) = G%Y(IY-1) + 0.5D0*(G%DLTY(1,1,IY-1) + G%DLTY(1,1,IY))
   ENDDO
   G%Y(NCELLY) = G%Y(NCELLY-1) + 0.5D0 * G%DLTY(1,1,NCELLY-1) + G%DLTY(1,1,NCELLY)    

ELSEIF ((NCELLX .GT. 1) .AND. (.NOT. G%HALF_CELLS_AT_BC)) THEN
   G%DLTY  (:,:,:) = YDIM/REAL(NCELLY-2,EB)
   !G%DLTYN (:,:,:) = G%DLTY(:,:,:)
   
   !Delta-y north and south:
   G%DYN(:,:,:) = YDIM/REAL(NCELLY-2,EB)
   G%DYS(:,:,:) = G%DYN(:,:,:) 

   ! y-position of the center of the cells
   G%Y(1) = -0.5D0 *G%DLTY (1,1,1) ! Center of the gost cell (not used)
   DO IY = 2, NCELLY
      G%Y(IY) = G%Y(IY-1) + 0.5D0*(G%DLTY(1,1,IY-1) + G%DLTY(1,1,IY))
   ENDDO

ELSE !NCELLY=1
   G%DLTY (:,:,:) = 1D0   
   !G%DLTYN(:,:,:) = 1D0
ENDIF

! Set Volume and surfaces
IF (NCELLX .EQ. 1 .AND. NCELLY .EQ. 1 .AND. NCELLZ .EQ. 1) THEN
   G%DXDY(1,1,1) = 1D0
   G%DXDZ(1,1,1) = 1D0
   G%DYDZ(1,1,1) = 1D0
   G%DV  (1,1,1) = 1D0 
ELSE
   DO IX = 1, NCELLX
   DO IY = 1, NCELLY
   DO IZ = 1, NCELLZ
      G%DXDY(IZ,IX,IY) = G%DLTX(IZ,IX,IY) * G%DLTY(IZ,IX,IY)
      G%DXDZ(IZ,IX,IY) = G%DLTX(IZ,IX,IY) * G%DLTZ(IZ,IX,IY)
      G%DYDZ(IZ,IX,IY) = G%DLTY(IZ,IX,IY) * G%DLTZ(IZ,IX,IY)
      G%DV  (IZ,IX,IY) = G%DLTZ(IZ,IX,IY) * G%DLTX(IZ,IX,IY) * G%DLTY(IZ,IX,IY) 
   ENDDO
   ENDDO
   ENDDO
ENDIF

!******************************************************************************	
END SUBROUTINE GENERATE_RAW_GRID
!******************************************************************************	


!******************************************************************************	
SUBROUTINE APPLY_GEOMETRY_TO_GRID(ICASE,IMESH)
!******************************************************************************	
!Set initial mass fractions, temperature, and pressure:
! IMASK, YIN, YI, TPN, TP,
! P, PN, YJG, YJGN, TG, TGN
! GPBCP%T_SURFACE, GPBCP%T_SURFACE_OLD
INTEGER, INTENT(IN) :: ICASE,IMESH

REAL(EB) :: XB(6)
INTEGER :: ICNUM,IOBST,ISPEC,IGSPEC
INTEGER :: SURF_IDX(1:6)
INTEGER :: NCELLZ, NCELLX, NCELLY
LOGICAL :: VALID_OBSTS,GEOMETRY_FILE_EXISTS
CHARACTER(300) :: MESSAGE
CHARACTER(4) :: FOUR

!Variables for read of geometry:
INTEGER :: N_OBST

! Variables for orientation
LOGICAL :: GO
INTEGER :: IPOS, ID,I , J, K, L, ICOUNT, IOS
REAL(EB) :: A11, A12, A13, A21, A22, A23, A31, A32, A33

CHARACTER(60) :: COLOR,SURF_ID
NAMELIST /OBST/ XB, COLOR, SURF_ID, ICNUM, SURF_IDX


G=>GPM(IMESH)
GPBCP(1:,1:,1:,-3:)=>G%GPYRO_BOUNDARY_CONDITION(1:,1:,1:,-3:)

NCELLZ = G%NCELLZ !Number of cells in z direction
NCELLX = G%NCELLX !Number of cells in x direction
NCELLY = G%NCELLY !Number of cells in y direction


!========================================================!
!=========== SET DEFAULTS FOR ENTIRE DOMAIN =============!
!========================================================!

XB(1)         = 0D0
XB(2)         = G%XDIM
XB(3)         = 0D0
XB(4)         = G%YDIM
XB(5)         = 0D0
XB(6)         = G%ZDIM
ICNUM         = GPG%DEFAULT_IC(IMESH)
SURF_IDX(:)   = GPG%DEFAULT_SURF_IDX(IMESH,:) 
CALL SETUP_DETAILED_ICS(.FALSE.,NCELLX,NCELLY,NCELLZ)

!Now set initial conditions obstruction by obstruction
ICNUM = GPG%DEFAULT_IC(IMESH)

DO IOBST = 1, GPG%NOBST
   IF (GPG%GEOM(IOBST)%IMESH .NE. IMESH) CYCLE
   
   DO ISPEC = 1, SPROP%NSSPEC
      G%YIN(ISPEC,:,:,:) = GPG%INITIAL_CONDITIONS(GPG%DEFAULT_IC(IMESH))%YI0(ISPEC)
      G%YI (ISPEC,:,:,:) = GPG%INITIAL_CONDITIONS(GPG%DEFAULT_IC(IMESH))%YI0(ISPEC)
   ENDDO

   G%TPN(:,:,:) = GPG%INITIAL_CONDITIONS(GPG%DEFAULT_IC(IMESH))%TMP_INITIAL
   G%TP (:,:,:) = GPG%INITIAL_CONDITIONS(GPG%DEFAULT_IC(IMESH))%TMP_INITIAL

   GPBCP(:,:,:,:)%T_SURFACE_OLD = GPG%INITIAL_CONDITIONS(GPG%DEFAULT_IC(IMESH))%TMP_INITIAL
   GPBCP(:,:,:,:)%T_SURFACE     = GPG%INITIAL_CONDITIONS(GPG%DEFAULT_IC(IMESH))%TMP_INITIAL

   IF (GPG%SOLVE_GAS_YJ) THEN
      DO IGSPEC = 1, GPROP%NGSPEC
         G%YJGN(IGSPEC,:,:,:) = GPG%INITIAL_CONDITIONS(GPG%DEFAULT_IC(IMESH))%YJ0(IGSPEC)
         G%YJG (IGSPEC,:,:,:) = GPG%INITIAL_CONDITIONS(GPG%DEFAULT_IC(IMESH))%YJ0(IGSPEC)
      ENDDO
   ENDIF

   IF (GPG%SOLVE_PRESSURE) THEN
      G%P  (:,:,:) = GPG%INITIAL_CONDITIONS(GPG%DEFAULT_IC(IMESH))%P_INITIAL
      G%PN (:,:,:) = GPG%INITIAL_CONDITIONS(GPG%DEFAULT_IC(IMESH))%P_INITIAL
   ENDIF

   IF (GPG%SOLVE_GAS_ENERGY) THEN
      G%TG  (:,:,:) = GPG%INITIAL_CONDITIONS(GPG%DEFAULT_IC(IMESH))%TMPG_INITIAL
      G%TGN (:,:,:) = GPG%INITIAL_CONDITIONS(GPG%DEFAULT_IC(IMESH))%TMPG_INITIAL
   ENDIF

ENDDO

IF (GPG%ZEROD(ICASE)) RETURN 

!Begin code for masking:
G%IMASK(:,:,:) = .TRUE.

!========================================================!
!============ SET GEOMETRY FROM OBST ====================!
!========================================================!

! Read info from &GPYRO_GEOM. Note that the OBST's specified in the input deck are read in
! first, and then can be "overwritten" with geometry from the geometry file. 
VALID_OBSTS = .FALSE.
DO IOBST = 1, GPG%NOBST 
   IF (GPG%GEOM(IOBST)%IMESH .NE. IMESH) CYCLE
   VALID_OBSTS = .TRUE.
   XB(1)       = GPG%GEOM(IOBST)%X1
   XB(2)       = GPG%GEOM(IOBST)%X2
   XB(3)       = GPG%GEOM(IOBST)%Y1
   XB(4)       = GPG%GEOM(IOBST)%Y2
   XB(5)       = GPG%GEOM(IOBST)%Z1
   XB(6)       = GPG%GEOM(IOBST)%Z2
   ICNUM       = GPG%GEOM(IOBST)%ICNUM
   SURF_IDX(:) = GPG%GEOM(IOBST)%SURF_IDX(:)
   CALL SETUP_DETAILED_ICS(.TRUE.,NCELLX,NCELLY,NCELLZ)
ENDDO

!========================================================!
!=========== READ GEOMTRY FILE IF PRESENT ===============!
!========================================================!
! Check to see if geometry file exists. If it does, read it in.
GEOMETRY_FILE_EXISTS = .TRUE.
OPEN(LUINPUT,FILE=TRIM(GPG%GEOMETRY_FILE(IMESH)),FORM='FORMATTED',STATUS='OLD',IOSTAT=IOS)
IF (IOS .GT. 0) THEN
   IF (TRIM(GPG%GEOMETRY_FILE(IMESH)) .EQ. 'null') THEN
      CONTINUE
   !ELSE
   !   IF (IGPYRO_TYPE .NE. 3) WRITE(*,*) TRIM(GPG%GEOMETRY_FILE(IMESH)), ' geometry file not found, skipping.'
   ENDIF
   GEOMETRY_FILE_EXISTS = .FALSE.
ENDIF

! Read geometry file
IF (GEOMETRY_FILE_EXISTS) THEN 
   N_OBST = 0
   DO
      SURF_IDX(1:6) = GPG%DEFAULT_SURF_IDX(IMESH,1:6)
      ICNUM         = GPG%DEFAULT_IC(IMESH)

      READ(LUINPUT,NML=OBST,END=1,ERR=2,IOSTAT=IOS)
      N_OBST = N_OBST + 1

      CALL SETUP_DETAILED_ICS(.TRUE.,NCELLX,NCELLY,NCELLZ)

      2 IF (IOS .GT. 0) THEN
         WRITE(FOUR,'(I4.4)') N_OBST + 1
         MESSAGE='ERROR: Problem with OBSTruction number' // FOUR
         CALL SHUTDOWN_GPYRO(MESSAGE)
      ENDIF
   ENDDO 
   1 REWIND(LUINPUT)

   CLOSE(LUINPUT)

ENDIF !Geometry file exists

! Check to see if orientation file exists. If it does, read it in.
IF (GEOMETRY_FILE_EXISTS) THEN
   G%ORIENTATION_FILE_EXISTS = .TRUE.
   IPOS = SCAN(TRIM(GPG%GEOMETRY_FILE(IMESH)),".", BACK = .TRUE.)
   G%ORIENTATION_FILE = GPG%GEOMETRY_FILE(IMESH)(1:IPOS) // "ori"

   OPEN(LUINPUT,FILE=TRIM(G%ORIENTATION_FILE),FORM='FORMATTED',STATUS='OLD',IOSTAT=IOS)
   IF (IOS .GT. 0) THEN
      IF (IGPYRO_TYPE .NE. 3) WRITE(*,*) TRIM(GPG%GEOMETRY_FILE(IMESH)), ' orientation file not found, skipping.'
      G%ORIENTATION_FILE_EXISTS = .FALSE.
   ENDIF
ELSE
   G%ORIENTATION_FILE_EXISTS = .FALSE. 
ENDIF

IF (G%ORIENTATION_FILE_EXISTS) THEN

   ALLOCATE (G%ORI       (1:NCELLZ, 1:NCELLX, 1:NCELLY, 1:3, 1:3)); G%ORI     (:,:,:,:,:) = 0D0
   ALLOCATE (G%K_TENSOR  (1:NCELLZ, 1:NCELLX, 1:NCELLY, 1:3, 1:3)); G%K_TENSOR(:,:,:,:,:) = 0D0

   ! Count number of lines
   GO     = .TRUE.
   ICOUNT = 0 
   DO WHILE (GO)
      READ(LUINPUT, *, IOSTAT=IOS)
      IF (IOS .EQ. 0) THEN
         ICOUNT = ICOUNT + 1
      ELSE
         GO = .FALSE.
      ENDIF  
   ENDDO
   ICOUNT = ICOUNT - 1

   REWIND (LUINPUT)
   READ (LUINPUT,*)
   DO L = 1, ICOUNT

      READ(LUINPUT,*) ID, I, J, K, A11, A12, A13, A21, A22, A23, A31, A32, A33

      G%ORI(K,I,J,1,1) = A11
      G%ORI(K,I,J,1,2) = A12
      G%ORI(K,I,J,1,3) = A13

      G%ORI(K,I,J,2,1) = A21
      G%ORI(K,I,J,2,2) = A22
      G%ORI(K,I,J,2,3) = A23

      G%ORI(K,I,J,3,1) = A31
      G%ORI(K,I,J,3,2) = A32
      G%ORI(K,I,J,3,3) = A33

      IF (ID .NE. L) THEN
         WRITE(*,*) 'Problem with orientation file'
         STOP
      ENDIF

   ENDDO

   CLOSE(LUINPUT)

ENDIF

! If no geometry info present, assume the entire domain extents are unmasked:
IF ((.NOT. GEOMETRY_FILE_EXISTS) .AND. (.NOT. VALID_OBSTS)) THEN
   XB(1)         = 0D0
   XB(2)         = G%XDIM
   XB(3)         = 0D0
   XB(4)         = G%YDIM
   XB(5)         = 0D0
   XB(6)         = G%ZDIM
   ICNUM         = GPG%DEFAULT_IC(IMESH)
   SURF_IDX(:)   = GPG%DEFAULT_SURF_IDX(IMESH,:) 
   CALL SETUP_DETAILED_ICS(.TRUE.,NCELLX,NCELLY,NCELLZ)
ENDIF


CONTAINS

! *****************************************************************************
SUBROUTINE SETUP_DETAILED_ICS(SET_IMASK,NCELLX,NCELLY,NCELLZ)
! *****************************************************************************

INTEGER, INTENT(IN) :: NCELLX,NCELLY,NCELLZ
LOGICAL, INTENT(IN) :: SET_IMASK

REAL(EB), DIMENSION(1:NCELLX) :: COORDX1, COORDX2
REAL(EB), DIMENSION(1:NCELLY) :: COORDY1, COORDY2
REAL(EB), DIMENSION(1:NCELLZ) :: COORDZ1, COORDZ2

INTEGER :: IX1,IX2,IY1,IY2,IZ1,IZ2
REAL(EB) :: XB1, XB2, XB3, XB4, XB5, XB6, TARGX1, TARGX2, TARGY1, TARGY2, TARGZ1, TARGZ2

XB1=XB(1)
XB2=XB(2)
XB3=XB(3)
XB4=XB(4)
XB5=XB(5)
XB6=XB(6)

IZ1 = 1
IZ2 = 1
IX1 = 1
IX2 = 1
IY1 = 1
IY2 = 1

IF (GPG%GEOMETRY_IS_UPSIDE_DOWN) THEN
   COORDZ1(:) = G%Z(:) - 0.5*G%DLTZN(:,1,1)
   TARGZ1     = G%ZDIM - XB6 
   IZ1        = IJK_FROM_XYZ(COORDZ1, NCELLZ, TARGZ1, 1)

   COORDZ2(:) = G%Z(:) + 0.5*G%DLTZN(:,1,1)
   TARGZ2     = G%ZDIM - XB5 
   IZ2        = IJK_FROM_XYZ(COORDZ2, NCELLZ, TARGZ2, 2)
ELSE
   COORDZ1(:) = G%Z(:) - 0.5*G%DLTZN(:,1,1)
   TARGZ1     = XB5 
   IZ1        = IJK_FROM_XYZ(COORDZ1, NCELLZ, TARGZ1, 1)

   COORDZ2(:) = G%Z(:) + 0.5*G%DLTZN(:,1,1)
   TARGZ2     = XB6
   IZ2        = IJK_FROM_XYZ(COORDZ2, NCELLZ, TARGZ2, 2)
ENDIF

IF (NCELLX .GT. 1) THEN 
   COORDX1(:) = G%X(:) - 0.5*G%DLTX (1,:,1)
   TARGX1     = XB1 
   IX1        = IJK_FROM_XYZ(COORDX1, NCELLX, TARGX1, 1)

   COORDX2(:) = G%X(:) + 0.5*G%DLTX (1,:,1)
   TARGX2     = XB2 
   IX2        = IJK_FROM_XYZ(COORDX2, NCELLX, TARGX2, 2)
ENDIF

IF (NCELLX .GT. 1 .AND. NCELLY .GT. 1) THEN
   COORDY1(:) = G%Y(:) - 0.5*G%DLTY(1,1,:)
   TARGY1     = XB3 
   IY1        = IJK_FROM_XYZ(COORDY1, NCELLY, TARGY1, 1)

   COORDY2(:) = G%Y(:) + 0.5*G%DLTY(1,1,:)
   TARGY2     = XB4 
   IY2        = IJK_FROM_XYZ(COORDY2, NCELLY, TARGY2, 2)
ENDIF

IF (SET_IMASK) THEN 
   G%IMASK(IZ1:IZ2,IX1:IX2,IY1:IY2) = .FALSE. ! Cells not masked
ENDIF


IF (SURF_IDX(5) .GT. 0) G%SURF_IDX_BCT(IZ1:IZ2,IX1:IX2,IY1:IY2) = SURF_IDX(5)
IF (SURF_IDX(6) .GT. 0) G%SURF_IDX_BCB(IZ1:IZ2,IX1:IX2,IY1:IY2) = SURF_IDX(6)
IF (SURF_IDX(1) .GT. 0) G%SURF_IDX_BCW(IZ1:IZ2,IX1:IX2,IY1:IY2) = SURF_IDX(1)
IF (SURF_IDX(2) .GT. 0) G%SURF_IDX_BCE(IZ1:IZ2,IX1:IX2,IY1:IY2) = SURF_IDX(2)
IF (SURF_IDX(3) .GT. 0) G%SURF_IDX_BCS(IZ1:IZ2,IX1:IX2,IY1:IY2) = SURF_IDX(3)
IF (SURF_IDX(4) .GT. 0) G%SURF_IDX_BCN(IZ1:IZ2,IX1:IX2,IY1:IY2) = SURF_IDX(4)

IF (ICNUM .GT. 0) THEN

   G%TPN(IZ1:IZ2,IX1:IX2,IY1:IY2) = GPG%INITIAL_CONDITIONS(ICNUM)%TMP_INITIAL
   G%TP (IZ1:IZ2,IX1:IX2,IY1:IY2) = GPG%INITIAL_CONDITIONS(ICNUM)%TMP_INITIAL

   GPBCP(IZ1:IZ2,IX1:IX2,IY1:IY2,:)%T_SURFACE_OLD = GPG%INITIAL_CONDITIONS(ICNUM)%TMP_INITIAL
   GPBCP(IZ1:IZ2,IX1:IX2,IY1:IY2,:)%T_SURFACE     = GPG%INITIAL_CONDITIONS(ICNUM)%TMP_INITIAL

   DO ISPEC = 1, SPROP%NSSPEC
      G%YIN(ISPEC,IZ1:IZ2,IX1:IX2,IY1:IY2) = GPG%INITIAL_CONDITIONS(ICNUM)%YI0(ISPEC)
      G%YI (ISPEC,IZ1:IZ2,IX1:IX2,IY1:IY2) = GPG%INITIAL_CONDITIONS(ICNUM)%YI0(ISPEC)
   ENDDO

   IF (GPG%NEED_GAS_YJ) THEN
      DO IGSPEC = 1, GPROP%NGSPEC
         G%YJGN(IGSPEC,IZ1:IZ2,IX1:IX2,IY1:IY2) = GPG%INITIAL_CONDITIONS(ICNUM)%YJ0(IGSPEC)
         G%YJG (IGSPEC,IZ1:IZ2,IX1:IX2,IY1:IY2) = GPG%INITIAL_CONDITIONS(ICNUM)%YJ0(IGSPEC)
      ENDDO
   ENDIF

   IF (GPG%SOLVE_PRESSURE) THEN
      G%PN (IZ1:IZ2,IX1:IX2,IY1:IY2) = GPG%INITIAL_CONDITIONS(ICNUM)%P_INITIAL
      G%P  (IZ1:IZ2,IX1:IX2,IY1:IY2) = GPG%INITIAL_CONDITIONS(ICNUM)%P_INITIAL
   ENDIF

   IF (GPG%SOLVE_GAS_ENERGY) THEN
      G%HGN(IZ1:IZ2,IX1:IX2,IY1:IY2) = GPG%INITIAL_CONDITIONS(ICNUM)%TMPG_INITIAL
      G%HG (IZ1:IZ2,IX1:IX2,IY1:IY2) = GPG%INITIAL_CONDITIONS(ICNUM)%TMPG_INITIAL
   ENDIF
ENDIF

! *****************************************************************************
END SUBROUTINE SETUP_DETAILED_ICS
! *****************************************************************************

!******************************************************************************	
END SUBROUTINE APPLY_GEOMETRY_TO_GRID
!******************************************************************************	


!******************************************************************************	
SUBROUTINE INIT_BOUNDARY_CELLS(IMESH)
!******************************************************************************	

INTEGER, INTENT(IN) :: IMESH
INTEGER :: IZ,IX,IY,NCELLZ,NCELLX,NCELLY
INTEGER :: ICOUNT

NCELLZ = G%NCELLZ
NCELLX = G%NCELLX
NCELLY = G%NCELLY

IF (.NOT. G%HALF_CELLS_AT_BC) THEN
   ! Mask the Gost cells
   IF (NCELLZ .GT. 1) G%IMASK(1     ,:,:) = .TRUE.
   IF (NCELLZ .GT. 1) G%IMASK(NCELLZ,:,:) = .TRUE.
   IF (NCELLX .GT. 1) G%IMASK(:,1     ,:) = .TRUE.
   IF (NCELLX .GT. 1) G%IMASK(:,NCELLX,:) = .TRUE.
   IF (NCELLY .GT. 1) G%IMASK(:,:,1     ) = .TRUE.
   IF (NCELLY .GT. 1) G%IMASK(:,:,NCELLY) = .TRUE.
ENDIF

! Determine cells that need boundary conditions.

! First Count the number of Boundary Cells
ICOUNT = 0

!z-direction
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 2, NCELLZ-1
   IF ((.NOT. G%IMASK(IZ,IX,IY)) .AND. G%IMASK(IZ-1,IX,IY))  ICOUNT = ICOUNT +1
   IF ((.NOT. G%IMASK(IZ,IX,IY)) .AND. G%IMASK(IZ+1,IX,IY))  ICOUNT = ICOUNT +1
ENDDO
ENDDO
ENDDO
!z-edges
DO IY = 1, NCELLY
DO IX = 1, NCELLX
   IZ=1     ; IF(.NOT. G%IMASK(IZ,IX,IY) ) ICOUNT = ICOUNT +1
   IZ=NCELLZ; IF(.NOT. G%IMASK(IZ,IX,IY) ) ICOUNT = ICOUNT +1
ENDDO
ENDDO

!x-direction
IF (NCELLX .GT. 1) THEN
   DO IY = 1, NCELLY
   DO IX = 2, NCELLX-1
   DO IZ = 1, NCELLZ
      IF ((.NOT. G%IMASK(IZ,IX,IY)) .AND. G%IMASK(IZ,IX-1,IY)) ICOUNT = ICOUNT +1
      IF ((.NOT. G%IMASK(IZ,IX,IY)) .AND. G%IMASK(IZ,IX+1,IY)) ICOUNT = ICOUNT +1
   ENDDO
   ENDDO
   ENDDO
   !x-edges
   DO IY = 1, NCELLY
   DO IZ = 1, NCELLZ
      IX=1     ; IF(.NOT. G%IMASK(IZ,IX,IY) ) ICOUNT = ICOUNT +1
      IX=NCELLX; IF(.NOT. G%IMASK(IZ,IX,IY) ) ICOUNT = ICOUNT +1
   ENDDO
   ENDDO
ENDIF

!y-direction
IF (NCELLY .GT. 1) THEN
   DO IY = 2, NCELLY-1
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
      IF ((.NOT. G%IMASK(IZ,IX,IY)) .AND. G%IMASK(IZ,IX,IY-1)) ICOUNT = ICOUNT +1
      IF ((.NOT. G%IMASK(IZ,IX,IY)) .AND. G%IMASK(IZ,IX,IY+1)) ICOUNT = ICOUNT +1
   ENDDO
   ENDDO
   ENDDO
   !y-edges
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
      IY=1     ; IF(.NOT. G%IMASK(IZ,IX,IY) ) ICOUNT = ICOUNT +1
      IY=NCELLY; IF(.NOT. G%IMASK(IZ,IX,IY) ) ICOUNT = ICOUNT +1
   ENDDO
   ENDDO
ENDIF

NGPYRO_FACES_NEEDING_BCS(IMESH) = ICOUNT

CALL ALLOCATE_GPYRO_BOUNDARYS_INFO(GPG%NUM_GPYRO_MESHES, IMESH, NGPYRO_FACES_NEEDING_BCS(IMESH))
ICOUNT=0

! Start with interior nodes:
!z-direction
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 2, NCELLZ-1
   IF ((.NOT. G%IMASK(IZ,IX,IY)) .AND. G%IMASK(IZ-1,IX,IY))  CALL INIT_BOUNDARY(ICOUNT,IMESH,IZ,IX,IY, 3,G%NEEDSBCT(IZ,IX,IY))
   IF ((.NOT. G%IMASK(IZ,IX,IY)) .AND. G%IMASK(IZ+1,IX,IY))  CALL INIT_BOUNDARY(ICOUNT,IMESH,IZ,IX,IY,-3,G%NEEDSBCB(IZ,IX,IY))
ENDDO
ENDDO
ENDDO

!x-direction
IF (NCELLX .GT. 1) THEN
   DO IY = 1, NCELLY
   DO IX = 2, NCELLX-1
   DO IZ = 1, NCELLZ
      IF ((.NOT. G%IMASK(IZ,IX,IY)) .AND. G%IMASK(IZ,IX-1,IY)) CALL INIT_BOUNDARY(ICOUNT,IMESH,IZ,IX,IY,-1,G%NEEDSBCW(IZ,IX,IY))
      IF ((.NOT. G%IMASK(IZ,IX,IY)) .AND. G%IMASK(IZ,IX+1,IY)) CALL INIT_BOUNDARY(ICOUNT,IMESH,IZ,IX,IY, 1,G%NEEDSBCE(IZ,IX,IY))
   ENDDO
   ENDDO
   ENDDO
ENDIF

!y-direction
IF (NCELLY .GT. 1) THEN
DO IY = 2, NCELLY-1
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ
   IF ((.NOT. G%IMASK(IZ,IX,IY)) .AND. G%IMASK(IZ,IX,IY-1)) CALL INIT_BOUNDARY(ICOUNT,IMESH,IZ,IX,IY,-2,G%NEEDSBCS(IZ,IX,IY))
   IF ((.NOT. G%IMASK(IZ,IX,IY)) .AND. G%IMASK(IZ,IX,IY+1)) CALL INIT_BOUNDARY(ICOUNT,IMESH,IZ,IX,IY, 2,G%NEEDSBCN(IZ,IX,IY))
ENDDO
ENDDO
ENDDO
ENDIF

!Now do edges:
DO IY = 1, NCELLY
DO IX = 1, NCELLX
   IZ=1     ; IF(.NOT. G%IMASK(IZ,IX,IY) ) CALL INIT_BOUNDARY(ICOUNT,IMESH,IZ,IX,IY, 3,G%NEEDSBCT(IZ,IX,IY))
   IZ=NCELLZ; IF(.NOT. G%IMASK(IZ,IX,IY) ) CALL INIT_BOUNDARY(ICOUNT,IMESH,IZ,IX,IY,-3,G%NEEDSBCB(IZ,IX,IY))
ENDDO
ENDDO

IF (NCELLX .GT. 1) THEN
   DO IY = 1, NCELLY
   DO IZ = 1, NCELLZ
      IX=1     ; IF(.NOT. G%IMASK(IZ,IX,IY) ) CALL INIT_BOUNDARY(ICOUNT,IMESH,IZ,IX,IY,-1,G%NEEDSBCW(IZ,IX,IY))
      IX=NCELLX; IF(.NOT. G%IMASK(IZ,IX,IY) ) CALL INIT_BOUNDARY(ICOUNT,IMESH,IZ,IX,IY, 1,G%NEEDSBCE(IZ,IX,IY))
   ENDDO
   ENDDO
ENDIF

IF (NCELLY .GT. 1) THEN
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
      IY=1     ; IF(.NOT. G%IMASK(IZ,IX,IY) ) CALL INIT_BOUNDARY(ICOUNT,IMESH,IZ,IX,IY,-2,G%NEEDSBCS(IZ,IX,IY))
      IY=NCELLY; IF(.NOT. G%IMASK(IZ,IX,IY) ) CALL INIT_BOUNDARY(ICOUNT,IMESH,IZ,IX,IY, 2,G%NEEDSBCN(IZ,IX,IY))
   ENDDO
   ENDDO
ENDIF

! For each cell get the ID of the bounday cell were the gas will escape
CALL GET_EXIT_GAS_CELL(IMESH)
CONTAINS
! *****************************************************************************
SUBROUTINE INIT_BOUNDARY(ICOUNTP,IMESHP,IZP,IXP,IYP,IORP,NEEDSBC)
! *****************************************************************************


INTEGER, INTENT(INOUT) :: ICOUNTP
INTEGER, INTENT(IN) :: IMESHP, IZP, IXP, IYP, IORP
LOGICAL, INTENT(OUT) :: NEEDSBC
REAL(EB) :: ZGPYRO, XGPYRO, YGPYRO
LOGICAL :: COMPLETE_CELL_AT_BC
G => GPM(IMESHP)

IF (ABS(IORP) .EQ. 1) THEN
   XGPYRO = G%X0 + G%X(IXP)
   ZGPYRO = G%Z0 + G%ZDIM - G%Z(IZP)
   YGPYRO = G%Y0 + G%Y(IYP)
   COMPLETE_CELL_AT_BC = .NOT. (IXP == 1 .OR. IXP == G%NCELLX)
ENDIF

IF (ABS(IORP) .EQ. 2) THEN
   YGPYRO = G%Y0 + G%Y(IYP)
   ZGPYRO = G%Z0 + G%ZDIM - G%Z(IZP)
   XGPYRO = G%X0 + G%X(IXP)
   COMPLETE_CELL_AT_BC = .NOT. (IYP == 1 .OR. IYP == G%NCELLY)

ENDIF

IF (ABS(IORP) .EQ. 3) THEN
   ZGPYRO = G%Z0 - G%Z(IZP) + G%ZDIM 
   YGPYRO = G%Y0 + G%Y(IYP)
   XGPYRO = G%X0 + G%X(IXP)
   COMPLETE_CELL_AT_BC = .NOT. (IZP == 1 .OR. IZP == NCELLZ)
ENDIF

ICOUNTP = ICOUNTP + 1
GP_BOUDARYS(IMESHP)%IMESH_GPYRO (ICOUNTP) = IMESHP
GP_BOUDARYS(IMESHP)%IZ_GPYRO    (ICOUNTP) = IZP
GP_BOUDARYS(IMESHP)%IX_GPYRO    (ICOUNTP) = IXP
GP_BOUDARYS(IMESHP)%IY_GPYRO    (ICOUNTP) = IYP
GP_BOUDARYS(IMESHP)%IOR_GPYRO   (ICOUNTP) = IORP
GP_BOUDARYS(IMESHP)%Z_GPYRO     (ICOUNTP) = ZGPYRO
GP_BOUDARYS(IMESHP)%X_GPYRO     (ICOUNTP) = XGPYRO
GP_BOUDARYS(IMESHP)%Y_GPYRO     (ICOUNTP) = YGPYRO
GP_BOUDARYS(IMESHP)%COMPLETE_CELL_AT_BC(ICOUNTP)= COMPLETE_CELL_AT_BC
NEEDSBC                      = .TRUE. 

! *****************************************************************************
END SUBROUTINE INIT_BOUNDARY
! *****************************************************************************

!******************************************************************************	
END SUBROUTINE INIT_BOUNDARY_CELLS
!******************************************************************************	



! *****************************************************************************
 SUBROUTINE INIT_GPYRO_VARS(IMESH)
! *****************************************************************************
INTEGER, INTENT(IN) :: IMESH

INTEGER :: NCELLZ, NCELLX, NCELLY, NSSPEC
INTEGER :: IZ, IX, IY, ISPEC
REAL(EB) :: HI, MWGN,RSUMYI

G=>GPM(IMESH)

NCELLZ = G%NCELLZ
NCELLX = G%NCELLX
NCELLY = G%NCELLY
NSSPEC = SPROP%NSSPEC


! Calculate ambient density (only used if gravity is present)
IF (IGPYRO_TYPE .NE. 2 .AND. GPG%SOLVE_GAS_YJ) THEN !In FDS sim we set GPG%RHOINF = RHOA in main.f90
   MWGN = 0D0 
   DO ISPEC = 1, GPROP%NGSPEC
      MWGN = MWGN + GPG%INITIAL_CONDITIONS(GPG%DEFAULT_IC(IMESH))%YJ0(ISPEC) / GPROP%M(ISPEC)
   ENDDO
   MWGN = 1D0 / MWGN !Mean molecular weight
   GPG%RHOINF = GPG%P0 * MWGN / (8314D0 * GPG%TAMB)
ENDIF

DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ

   G%HP(IZ,IX,IY) = 0D0
   DO ISPEC = 1, NSSPEC
      HI = HOFT(ISPEC,G%TP(IZ,IX,IY) )
      G%HP(IZ,IX,IY) = G%HP(IZ,IX,IY) + G%YI(ISPEC,IZ,IX,IY)*HI
   ENDDO

   ! Calculate the weighted condensed-phase density and volume fractions
   RSUMYI  = 0D0 !For bulk density
   DO ISPEC = 1, NSSPEC
      IF (G%YIN(ISPEC,IZ,IX,IY) .EQ. 0) CYCLE
      RSUMYI = RSUMYI  + G%YIN(ISPEC,IZ,IX,IY)/RHOOFT(ISPEC,G%TP(IZ,IX,IY))
   ENDDO
   G%RP (IZ,IX,IY) = 1D0 / RSUMYI

   DO ISPEC = 1, NSSPEC
      G%XI(ISPEC,IZ,IX,IY) = G%RP(IZ,IX,IY)* G%YI(ISPEC,IZ,IX,IY) / RHOOFT(ISPEC,G%TP(IZ,IX,IY))
   ENDDO

   G%RPN (IZ,IX,IY)   = G%RP (IZ,IX,IY)
   G%XIN (:,IZ,IX,IY) = G%XI (:,IZ,IX,IY)

   ! Calculate rho*D
   G%RDLTZ (IZ,IX,IY) = G%RP   (IZ,IX,IY) * G%DLTZ(IZ,IX,IY)
   G%RDLTZN(IZ,IX,IY) = G%RDLTZ(IZ,IX,IY)

   ! Setup for TGA simulation
   IF (NCELLZ .EQ. 1) THEN
      G%DLTZ(IZ,IX,IY)   = 1D0/G%RP(IZ,IX,IY)
      G%DLTZN(IZ,IX,IY)  = G%DLTZ(IZ,IX,IY)
      G%RDLTZ(IZ,IX,IY)  = G%RP(IZ,IX,IY) * G%DLTZ(IZ,IX,IY)
      G%RDLTZN(IZ,IX,IY) = G%RDLTZ(IZ,IX,IY)
      G%IMASK(IZ,IX,IY) = .FALSE.
   ENDIF

   ! Calculate rho*Yi*Dz
   DO ISPEC = 1, NSSPEC
      G%RYIDZP(ISPEC,IZ,IX,IY) = G%RDLTZ(IZ,IX,IY) * G%YI(ISPEC,IZ,IX,IY)
   ENDDO
   G%RYIDZPN(:,IZ,IX,IY) = G%RYIDZP(:,IZ,IX,IY)

   ! Calculate initial mass of each species in each cell
   DO ISPEC = 1, NSSPEC
      G%RYIDZ0(ISPEC,IZ,IX,IY) = G%RP(IZ,IX,IY)*G%YI(ISPEC,IZ,IX,IY)*G%DLTZ(IZ,IX,IY)
   ENDDO
   G%RYIDZ0(0,IZ,IX,IY) = G%RP(IZ,IX,IY)*G%DLTZ(IZ,IX,IY)
   ! Initialize rho*Yi*Dz sigma (star)	     
   G%RYIDZSIGMA (:,IZ,IX,IY) = G%RYIDZ0(1:NSSPEC,IZ,IX,IY)
   G%RYIDZSIGMAN(:,IZ,IX,IY) = G%RYIDZ0(1:NSSPEC,IZ,IX,IY)

   IF (GPG%CONVENTIONAL_RXN_ORDER) THEN
      DO ISPEC = 1, NSSPEC
         G%RYIDZSIGMA (ISPEC,IZ,IX,IY) = G%RDLTZ(IZ,IX,IY)
         G%RYIDZSIGMAN(ISPEC,IZ,IX,IY) = G%RDLTZ(IZ,IX,IY)
      ENDDO
   ENDIF
         
   ! Set weighted thermal conductivity and specific heat
   G%CPS(IZ,IX,IY) = 0D0;  G%KZ (IZ,IX,IY) = 0D0
   DO ISPEC = 1, NSSPEC
      G%KZ (IZ,IX,IY) =  G%KZ (IZ,IX,IY)  + G%XI(ISPEC,IZ,IX,IY) * KZOFT(ISPEC,G%TP(IZ,IX,IY))
      G%CPS(IZ,IX,IY) =  G%CPS(IZ,IX,IY)  + G%YI(ISPEC,IZ,IX,IY) * CPOFT(ISPEC,G%TP(IZ,IX,IY))
      G%KAPPA(IZ,IX,IY)= G%KAPPA(IZ,IX,IY) +G%XI(ISPEC,IZ,IX,IY) * SPROP%KAPPA(ISPEC)
   ENDDO
 
   IF (NCELLX .GT. 1) THEN
      G%KX(IZ,IX,IY) = 0D0
      DO ISPEC = 1, NSSPEC
         G%KX(IZ,IX,IY) = G%KX(IZ,IX,IY) + G%XI(ISPEC,IZ,IX,IY) * KXOFT(ISPEC,G%TP(IZ,IX,IY))
      ENDDO
   ENDIF

   IF (NCELLY .GT. 1) THEN
      G%KY(IZ,IX,IY) = 0D0
      DO ISPEC = 1, NSSPEC
         G%KY(IZ,IX,IY) = G%KY(IZ,IX,IY) + G%XI(ISPEC,IZ,IX,IY) * KYOFT(ISPEC,G%TP(IZ,IX,IY))
      ENDDO
   ENDIF


   ! Initialize gas-phase quantities
   IF (GPG%SOLVE_GAS_YJ) THEN   
      MWGN = 0D0 
      DO ISPEC = 1, GPROP%NGSPEC
         MWGN = MWGN + G%YJG(ISPEC,IZ,IX,IY) / GPROP%M(ISPEC)
      ENDDO
      MWGN = 1D0 / MWGN !Mean molecular weight
      ! Set hydrostatic pressure if gravity is present
      IF (GPG%SOLVE_PRESSURE) THEN
         G%RGN(IZ,IX,IY) = RHOGOFT(G%PN(IZ,IX,IY), MWGN, G%TPN(IZ,IX,IY))
         G%RG (IZ,IX,IY) = G%RGN(IZ,IX,IY)

         IF (G%GZ .NE. 0D0) THEN
            G%P (IZ,IX,IY) = G%P(IZ,IX,IY) + GPG%RHOINF*G%GZ*G%Z(IZ)
            G%PN(IZ,IX,IY) = G%P(IZ,IX,IY)
         ENDIF
         IF (G%GX .NE. 0D0) THEN
            G%P (IZ,IX,IY) = G%P(IZ,IX,IY) + GPG%RHOINF*G%GX*G%X(IX)
            G%PN(IZ,IX,IY) = G%P(IZ,IX,IY)
         ENDIF
         IF (G%GY .NE. 0D0) THEN
            G%P (IZ,IX,IY) = G%P(IZ,IX,IY) + GPG%RHOINF*G%GY*G%Y(IY)
            G%PN(IZ,IX,IY) = G%P(IZ,IX,IY)
         ENDIF
      ELSE
         G%RG (IZ,IX,IY) = RHOGOFT(GPG%P0, MWGN, G%TPN(IZ,IX,IY))
         G%RGN(IZ,IX,IY) = G%RG(IZ,IX,IY)
      ENDIF
                           
      IF (GPG%SOLVE_GAS_ENERGY) THEN

         G%TG (IZ,IX,IY) = G%TP(IZ,IX,IY) !Gas temperature
         G%TGN(IZ,IX,IY) = G%TG(IZ,IX,IY) 
                  
         G%HG (IZ,IX,IY) = GPROP%CPG * (G%TG(IZ,IX,IY)-GPG%TDATUM) !Gas enthalpy
         G%HGN(IZ,IX,IY) = G%HG(IZ,IX,IY)
      ENDIF

   ENDIF !SOLVE_GAS_YJ 

   IF (GPG%SOLVE_PRESSURE) THEN
      DO ISPEC = 1, NSSPEC
         G%PERMZ(IZ,IX,IY) = G%PERMZ(IZ,IX,IY) + G%XI(ISPEC,IZ,IX,IY)*SPROP%PERMZ(ISPEC)
      ENDDO
   ENDIF
 

   IF (GPG%SOLVE_POROSITY) THEN
      RSUMYI = 0D0
      DO ISPEC = 1, NSSPEC
         RSUMYI = RSUMYI + G%YI(ISPEC,IZ,IX,IY) / SPROP%RS0(ISPEC)
      ENDDO
      G%RSP(IZ,IX,IY) = 1D0 / RSUMYI
      G%POROSS(IZ,IX,IY) = 1D0 - G%RP(IZ,IX,IY) / G%RSP(IZ,IX,IY)

      G%RSPN(IZ,IX,IY)    = G%RSP(IZ,IX,IY)
      G%POROSSN(IZ,IX,IY) = G%POROSS(IZ,IX,IY)
   ENDIF

ENDDO !IZ
ENDDO !IX
ENDDO !IY



!Record initial mass in mesh:
G%INITIAL_MASS = TOTAL_MASS(0)

GPG%MAXTMP = MAXVAL(G%TPN(:,:,:))
IF (GPG%SOLVE_GAS_ENERGY .AND. (.NOT. GPG%THERMAL_EQUILIBRIUM)) THEN
   GPG%MAXTMP = MAX(GPG%MAXTMP,MAXVAL(G%TGN(:,:,:)))
ENDIF

! Calculate variables for newton iteration
DO ISPEC = 1, NSSPEC
   GPG%NEWTON_B(ISPEC) = 1D0 + SPROP%NC(ISPEC)
   GPG%NEWTON_A(ISPEC) = SPROP%C0(ISPEC) / (GPG%NEWTON_B(ISPEC) * GPG%TREF**SPROP%NC(ISPEC))
   GPG%NEWTON_C(ISPEC) = GPG%NEWTON_A(ISPEC) * GPG%TDATUM**(SPROP%NC(ISPEC) + 1D0)
   GPG%NEWTON_D(ISPEC) = 0.5D0 * SPROP%DHMELT(ISPEC)
   GPG%NEWTON_E(ISPEC) = SPROP%TMELT(ISPEC)
   GPG%NEWTON_F(ISPEC) = MAX(SQRT(2D0*SPROP%SIGMA2MELT(ISPEC)),1D-10)
   GPG%NEWTON_G(ISPEC) = GPG%NEWTON_D(ISPEC) * DERF((GPG%TDATUM - SPROP%TMELT(ISPEC)) / GPG%NEWTON_F(ISPEC))
ENDDO

! *****************************************************************************
END SUBROUTINE INIT_GPYRO_VARS
! *****************************************************************************


! *****************************************************************************
SUBROUTINE INIT_REACTIONS
! *****************************************************************************
INTEGER :: IRXN, NRXN, NHGRXN
INTEGER :: ISPEC, NSSPEC, NGSPEC
LOGICAL :: FOUND_FROM, FOUND_TO, FOUND_REACTANT1, FOUND_REACTANT2
CHARACTER(300) :: MESSAGE
CHARACTER(3) :: THREE

NSSPEC =  SPROP%NSSPEC !Number of solid species
NGSPEC =  GPROP%NGSPEC !Number of gaseous species
NRXN   =  SPROP%NRXN   !Number of condensed-phase reactions
NHGRXN =  GPROP%NHGRXN !Number of homogeneous gaseous reactions

! Figure out condensed-phase reaction IFROM and ITO
DO IRXN = 1, NRXN
   FOUND_FROM = .FALSE.
   FOUND_TO   = .FALSE. 
   DO ISPEC = 1, NSSPEC
      IF (TRIM(RXN(IRXN)%CFROM).EQ.TRIM(SPROP%NAME(ISPEC)))THEN
         RXN(IRXN)%IFROM = ISPEC
         FOUND_FROM = .TRUE.
      ENDIF

      IF (TRIM(RXN(IRXN)%CTO  ).EQ.TRIM(SPROP%NAME(ISPEC)))THEN
         FOUND_TO = .TRUE.
         RXN(IRXN)%ITOS = ISPEC
      ENDIF

   ENDDO !ISPEC

   IF (.NOT. FOUND_FROM) THEN
      WRITE(THREE,'(I3.3)') IRXN
      MESSAGE='Error: for IRXN= ' // THREE // ', FROM species not found.'
      CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF

   IF (.NOT. FOUND_TO   ) THEN
      IF (TRIM(RXN(IRXN)%CTO) .EQ. 'gases' .OR. TRIM(RXN(IRXN)%CTO) .EQ. 'GASES') THEN
         IF (IGPYRO_TYPE .NE. 3) WRITE(0,"(1X,A,I2,A)") "Setting up reaction ", IRXN, " as noncharring because TO='gases'"
         RXN(IRXN)%ITOS = 0
      ELSE
         WRITE(THREE,'(I3.3)') IRXN
         MESSAGE='Error: for IRXN= ' // THREE // "TO species not found. Specify either a condensed phase product or 'gases' for a noncharring reaction."
         CALL SHUTDOWN_GPYRO(MESSAGE)
      ENDIF
   ENDIF

ENDDO !IRXN

! Figure out homogeneous gas reactions reactant indices:
DO IRXN = 1, NHGRXN
   FOUND_REACTANT1 = .FALSE.
   FOUND_REACTANT2 = .FALSE. 

   DO ISPEC = 1, NGSPEC

      IF (TRIM(HGRXN(IRXN)%CREACTANT1) .EQ. TRIM (GPROP%NAME(ISPEC))) THEN
         HGRXN(IRXN)%IREACTANT1 = ISPEC
         FOUND_REACTANT1 = .TRUE.
      ENDIF

      IF (TRIM(HGRXN(IRXN)%CREACTANT2) .EQ. TRIM (GPROP%NAME(ISPEC))) THEN
         FOUND_REACTANT2 = .TRUE.
         HGRXN(IRXN)%IREACTANT2 = ISPEC
      ENDIF

   ENDDO !ISPEC

   IF (.NOT. FOUND_REACTANT1) THEN
      WRITE(THREE,'(I3.3)') IRXN
         MESSAGE='Error: for IHGRXN = '  // THREE // 'REACTANT1 not found.'
      CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF

   IF (.NOT. FOUND_REACTANT2) THEN
      WRITE(THREE,'(I3.3)') IRXN
      MESSAGE='Error: for IHGRXN = '  // THREE // 'REACTANT2 not found.'
      CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF
               
   ! Check the yield of reactant 1:
   IF (GPROP%HGYIELDS(HGRXN(IRXN)%IREACTANT1,IRXN) .NE. -1D0) THEN
      WRITE(THREE,'(I3.3)') IRXN
      MESSAGE='Error: for IHGRXN = '  // THREE // ' hgyield of reactant 1 must be -1'
      CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF

ENDDO !IRXN

! *****************************************************************************
END SUBROUTINE INIT_REACTIONS
! *****************************************************************************



! *****************************************************************************
 SUBROUTINE INIT_DUMP_OUTPOUT
! *****************************************************************************

INTEGER :: I, N, IMESH,J


! Set FIRST_SF_IN_MESH and FIRST_PROF_IN_MESH
GPG%FIRST_SF_IN_MESH(:) = -1
DO J = 1, 100 !J is mesh #
   DO I = 1, GPG%N_SMOKEVIEW_QUANTITIES
      IF ((GPG%SMOKEVIEW_IMESH(I) .EQ. J .OR. GPG%SMOKEVIEW_IMESH(I) .EQ. 0) .AND. GPG%FIRST_SF_IN_MESH  (J) .LT. 1 ) GPG%FIRST_SF_IN_MESH  (J) = I
   ENDDO
ENDDO

GPG%FIRST_PROF_IN_MESH(:) = -1
DO J = 1, 100 !J is mesh #
   DO I = 1, GPG%N_PROFILE_QUANTITIES
      IF ((GPG%PROFILE_IMESH  (I) .EQ. J .OR. GPG%PROFILE_IMESH  (I)  .EQ. 0) .AND. GPG%FIRST_PROF_IN_MESH(J) .LT. 1) GPG%FIRST_PROF_IN_MESH(J) = I
   ENDDO
ENDDO


! Initialize point dumps:
DO I = 1, GPG%N_POINT_QUANTITIES
   IMESH  = GPG%POINT_IMESH(I)
   IF (IMESH .EQ. 0) IMESH=1
   G => GPM(IMESH)

   ! Start with z:
   CALL GET_CLOSEST_CELL(G%NCELLZ,G%Z,GPG%POINT_Z(I),GPG%POINT_IZ(I))
   ! Now on to x:
   IF (G%NCELLX .EQ. 1) THEN
      GPG%POINT_IX(I) = 1
   ELSE
      CALL GET_CLOSEST_CELL(G%NCELLX,G%X,GPG%POINT_X(I),GPG%POINT_IX(I))
   ENDIF
   
   ! And on to y:
   IF (G%NCELLY .EQ. 1) THEN
      GPG%POINT_IY(I) = 1
   ELSE
      CALL GET_CLOSEST_CELL(G%NCELLY,G%Y,GPG%POINT_Y(I),GPG%POINT_IY(I))
   ENDIF
ENDDO

! Initialize profile dumps:
DO I = 1, GPG%N_PROFILE_QUANTITIES
   IMESH  = GPG%PROFILE_IMESH(I)
   IF (IMESH .EQ. 0) IMESH=1
   G => GPM(IMESH)
   IF (GPG%PROFILE_DIRECTION(I) .EQ. 'z') THEN
      CALL GET_CLOSEST_CELL(G%NCELLX,G%X,GPG%PROFILE_COORD1(I),GPG%PROFILE_IX(I))
      CALL GET_CLOSEST_CELL(G%NCELLY,G%Y,GPG%PROFILE_COORD2(I),GPG%PROFILE_IY(I))   

   ELSEIF (GPG%PROFILE_DIRECTION(I) .EQ. 'x') THEN
      CALL GET_CLOSEST_CELL(G%NCELLY,G%Y,GPG%PROFILE_COORD1(I),GPG%PROFILE_IY(I))   
      CALL GET_CLOSEST_CELL(G%NCELLZ,G%Z,GPG%PROFILE_COORD2(I),GPG%PROFILE_IZ(I))

   ELSEIF (GPG%PROFILE_DIRECTION(I) .EQ. 'y') THEN
      CALL GET_CLOSEST_CELL(G%NCELLX,G%X,GPG%PROFILE_COORD1(I),GPG%PROFILE_IX(I))   
      CALL GET_CLOSEST_CELL(G%NCELLZ,G%Z,GPG%PROFILE_COORD2(I),GPG%PROFILE_IZ(I))
   ENDIF
ENDDO


!Initialize Smokeview dumps:
DO N = 1, GPG%N_SMOKEVIEW_QUANTITIES
   IMESH  = GPG%SMOKEVIEW_IMESH(N)
   IF (IMESH .EQ. 0) IMESH=1
   G => GPM(IMESH)
   IF (GPG%SMOKEVIEW_PLANE(N) .EQ. 'xz') CALL GET_CLOSEST_CELL(G%NCELLY,G%Y,GPG%SMOKEVIEW_LOCATION(N),GPG%SMOKEVIEW_ICELL(N)) !y=const plane
   IF (GPG%SMOKEVIEW_PLANE(N) .EQ. 'yz') CALL GET_CLOSEST_CELL(G%NCELLX,G%X,GPG%SMOKEVIEW_LOCATION(N),GPG%SMOKEVIEW_ICELL(N)) !x=const plane
   IF (GPG%SMOKEVIEW_PLANE(N) .EQ. 'xy') CALL GET_CLOSEST_CELL(G%NCELLZ,G%Z,GPG%SMOKEVIEW_LOCATION(N),GPG%SMOKEVIEW_ICELL(N)) !z=const plane
ENDDO


! *****************************************************************************
END SUBROUTINE INIT_DUMP_OUTPOUT
! *****************************************************************************

!******************************************************************************	
END MODULE GPYRO_INIT
!******************************************************************************