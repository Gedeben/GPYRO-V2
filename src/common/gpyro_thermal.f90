! *****************************************************************************
MODULE GPYRO_SOLID_THERMAL
! *****************************************************************************

USE PREC
USE GPYRO_VARS
USE GPYRO_FUNCS
USE GPYRO_BC, ONLY : SET_BC_FIXED_VALUE
IMPLICIT NONE

CONTAINS


! *****************************************************************************
SUBROUTINE SOLID_ENTHALPY_SOLVER(IMESH,DTIME)
! *****************************************************************************
INTEGER, INTENT(IN) :: IMESH
REAL(EB), INTENT(IN) :: DTIME
REAL(EB) :: TSTART,TMID1,TMID2
REAL(EB), POINTER, DIMENSION(:,:,:) :: AP,AB,AT,AE,AW,AN,AS,B
REAL(EB), POINTER, DIMENSION(:,:,:) :: PTR
INTEGER :: NCELLX,NCELLY,NCELLZ,ISSPEC

CALL GET_CPU_TIME(TSTART)
G => GPM(IMESH)

AP      => G%RWORK02
AB      => G%RWORK03
AT      => G%RWORK04
AE      => G%RWORK05
AW      => G%RWORK06
AN      => G%RWORK07
AS      => G%RWORK08
B       => G%RWORK09

PTR     => G%HPN

NCELLX = G%NCELLX
NCELLY = G%NCELLY
NCELLZ = G%NCELLZ

!================= THERMAL_PROPERTIES =======================!
! Get conductivity, heat capacity, absorption, specific enthalpy
CALL GET_CPU_TIME(TMID1)
CALL CALCULATE_WEIGHTED_POINT_THERMAL_PROPERTIES(IMESH)
CALL GET_CPU_TIME(TMID2);  GPG%TUSED(9) = GPG%TUSED(9) + TMID2 - TMID1

!==================== SOURCE TERMS =========================!
CALL SOLID_ENERGY_SOURCE_TERMS(IMESH,NCELLZ,NCELLX,NCELLY,G%DTIME_TMP)
CALL GET_CPU_TIME(TMID1);  GPG%TUSED(2) = GPG%TUSED(2) + TMID1 - TMID2


!============= Specify coefficients for TDMA ===============!
IF (GPG%USE_ANISOTROPIC_SOLID_ENTHALPY_SOLVER) THEN
   CALL GET_COEFF_OF_SOLID_ENTHALPY_SOLVER_ANISOTROPIC(DTIME,B,AP,AB,AT,AE,AW,AN,AS)
ELSE
   CALL GET_COEFF_OF_SOLID_ENTHALPY_SOLVER      (IMESH,DTIME,B,AP,AB,AT,AE,AW,AN,AS)
ENDIF
CALL GET_CPU_TIME(TMID2)
GPG%TUSED(3) = GPG%TUSED(3) + TMID2 - TMID1 !CPU time to build coefficents for TDMA


!=================== SHYI CORRECTION =======================!

IF (GPG%SHYI_CORRECTION ) THEN
   DO ISSPEC = 1, SPROP%NSSPEC
      CALL CALCULATE_THERMAL_INTERFACE_QUANTITIES(IMESH,'KOCHI    ',ISSPEC)
   ENDDO      
   CALL HANDLE_SHYI_CORRECTION(IMESH,B)
ENDIF

CALL GET_CPU_TIME(TMID1)
GPG%TUSED(4) = GPG%TUSED(4) + TMID1 - TMID2  ! time for SHYI CORRECTION Part 2



!=============== BOUNDARY CONDITIONS =======================!
CALL BC_FOR_SOLID_ENTHALPY_SOLVER(IMESH,G%SHNET,AP,AB,AT,AE,AW,AN,AS,B)
CALL GET_CPU_TIME(TMID2)
GPG%TUSED(5) = GPG%TUSED(5) + TMID2 - TMID1 ! Time to handle BC

!= SEPARATE POSITVE AND NEGATIVE SOURCE TERME FOR THE SOLVER =!
CALL HANDLE_POSITIVE_AND_NEGATIVE_SOURCE_TERM(G%SHNET,B)!,AP)
CALL GET_CPU_TIME(TMID1)
GPG%TUSED(3) = GPG%TUSED(3) + TMID1 - TMID2  

!========= SOLVE TDMA ACCORDING SWEEP DIRECTION ==========!
CALL SOLVE_TDMA_ACCORDING_SWEEP_DIRECTION(AP,AE,AW,AS,AN,AT,AB,B,PTR)
CALL GET_CPU_TIME(TMID2)
GPG%TUSED(6) = GPG%TUSED(6) + TMID2 - TMID1 ! Time for TDMA solver

!============ GET TEMPERATURE FROM ENTHALPY ==============!
CALL GET_T_FROM_H(IMESH)
CALL GET_CPU_TIME(TMID1)
GPG%TUSED(7) = GPG%TUSED(7) + TMID1 - TMID2  ! Time to Get T from H

!================= CHECK CONVERGENCE =====================!
CALL CHECK_T_AND_H_CONVERGENCE(IMESH)

! *****************************************************************************
END SUBROUTINE SOLID_ENTHALPY_SOLVER
! *****************************************************************************


! *****************************************************************************
SUBROUTINE IN_DEPTH_RADIATION_ABSORPTION(QRADNET, SHP_RAD, IZ, IX, IY, IOR, IMESH)
! *****************************************************************************

REAL(EB), INTENT(IN   ) :: QRADNET
INTEGER , INTENT(IN   ) :: IZ, IX, IY, IOR, IMESH
REAL(EB), INTENT(INOUT), POINTER, DIMENSION (:,:,:) :: SHP_RAD

REAL(EB), POINTER    , DIMENSION(:)  ::  DELTA, DQRDN, KAPPA
LOGICAL , POINTER    , DIMENSION(:)  ::  IMASK
INTEGER  :: NCELL, I, I_IN, I_OUT, ISIGN
REAL(EB) :: Q1, Q2, DQ

G => GPM(IMESH)

SELECT CASE (ABS(IOR)) ! Configure specific parameters depending on the direction
   CASE (1) ! x direction
      I_IN      =  IX  
      NCELL    =  G%NCELLX
      DELTA    => G%DLTX (IZ,:,IY)
      IMASK    => G%IMASK(IZ,:,IY)
      DQRDN    => SHP_RAD(IZ,:,IY)
      KAPPA    => G%KAPPA(IZ,:,IY)
    
   CASE (2) ! y direction
      I_IN     =  IY
      NCELL    =  G%NCELLY
      DELTA    => G%DLTY (IZ,IX,:)
      IMASK    => G%IMASK(IZ,IX,:)
      DQRDN    => SHP_RAD(IZ,IX,:)
      KAPPA    => G%KAPPA(IZ,IX,:)

   CASE (3) ! z direction
      I_IN     =  IZ
      NCELL    =  G%NCELLZ
      DELTA    => G%DLTZN(:,IX,IY)
      IMASK    => G%IMASK(:,IX,IY)
      DQRDN    => SHP_RAD(:,IX,IY)
      KAPPA    => G%KAPPA(:,IX,IY)
END SELECT

! Configure the direction of the absorption
IF ((IOR== 3) .OR. (IOR==-2) .OR. (IOR==-1)) THEN ; I_OUT = NCELL; ISIGN = 1 ; ENDIF
IF ((IOR==-3) .OR. (IOR== 2) .OR. (IOR== 1)) THEN ; I_OUT =  1   ; ISIGN =-1 ; ENDIF
Q1 = QRADNET

DO I = I_IN, I_OUT, ISIGN
   IF (IMASK(I))  EXIT !Reached end of the sample
   IF (Q1 .LT. EPSILON_FB) EXIT !No more Energy
   Q2 = Q1 * EXP(-KAPPA(I)*DELTA(I))
   DQ = Q1 - Q2
   DQRDN(I) = DQRDN(I) + DQ/DELTA(I)
   Q1 = Q2
ENDDO

! *****************************************************************************
END SUBROUTINE IN_DEPTH_RADIATION_ABSORPTION
! *****************************************************************************


! *****************************************************************************
SUBROUTINE CALCULATE_WEIGHTED_POINT_THERMAL_PROPERTIES(IMESH)
! *****************************************************************************

INTEGER, INTENT(IN) :: IMESH
INTEGER :: IZ,IX,IY,ISPEC,I,J
INTEGER :: NCELLX, NCELLY, NCELLZ,NSSPEC


G => GPM(IMESH)

NCELLX = G%NCELLX
NCELLY = G%NCELLY
NCELLZ = G%NCELLZ
NSSPEC = SPROP%NSSPEC



! Calculate conductivity
IF (G%ORIENTATION_FILE_EXISTS) THEN
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ,ISPEC) SHARED(G,SPROP,NCELLX,NCELLY,NCELLZ,NSSPEC) COLLAPSE(3)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ

      G%KZ (IZ,IX,IY) = 0D0
      G%KX (IZ,IX,IY) = 0D0 
      G%KY (IZ,IX,IY) = 0D0 
      G%K_TENSOR(IZ,IX,IY,:,:) = 0D0
      IF (G%IMASK(IZ,IX,IY)) CYCLE

      DO ISPEC = 1, NSSPEC

         ! Thermal conductivity
         IF (GPG%ANISOTROPIC_SPECIES(ISPEC)) THEN
            DO I = 1, 3 
            DO J = 1, 3
               G%K_TENSOR(IZ,IX,IY,I,J) = G%XIN(ISPEC,IZ,IX,IY) * ( G%ORI(IZ,IX,IY,1,I) * G%ORI(IZ,IX,IY,1,J) * G%XIN(ISPEC,IZ,IX,IY) * KXOFT(ISPEC,G%TP(IZ,IX,IY)) + & 
                                                                     G%ORI(IZ,IX,IY,2,I) * G%ORI(IZ,IX,IY,2,J) * G%XIN(ISPEC,IZ,IX,IY) * KYOFT(ISPEC,G%TP(IZ,IX,IY)) + &
                                                                     G%ORI(IZ,IX,IY,3,I) * G%ORI(IZ,IX,IY,3,J) * G%XIN(ISPEC,IZ,IX,IY) * KZOFT(ISPEC,G%TP(IZ,IX,IY))     )
            ENDDO
            ENDDO
         ELSE 
            G%KZ (IZ,IX,IY) = G%KZ (IZ,IX,IY) + G%XIN(ISPEC,IZ,IX,IY)*KZOFT(ISPEC,G%TP(IZ,IX,IY))
            G%KX (IZ,IX,IY) = G%KX (IZ,IX,IY) + G%XIN(ISPEC,IZ,IX,IY)*KXOFT(ISPEC,G%TP(IZ,IX,IY))
            G%KY (IZ,IX,IY) = G%KY (IZ,IX,IY) + G%XIN(ISPEC,IZ,IX,IY)*KYOFT(ISPEC,G%TP(IZ,IX,IY))
         ENDIF

      ENDDO !ISPEC

      G%K_TENSOR(IZ,IX,IY,1,1) = G%K_TENSOR(IZ,IX,IY,1,1) + G%KX(IZ,IX,IY)
      G%K_TENSOR(IZ,IX,IY,2,2) = G%K_TENSOR(IZ,IX,IY,2,2) + G%KY(IZ,IX,IY)
      G%K_TENSOR(IZ,IX,IY,3,3) = G%K_TENSOR(IZ,IX,IY,3,3) + G%KZ(IZ,IX,IY)

      ! This is a temporary hack:
      G%KX(IZ,IX,IY) = G%K_TENSOR(IZ,IX,IY,1,1)
      G%KY(IZ,IX,IY) = G%K_TENSOR(IZ,IX,IY,2,2)
      G%KZ(IZ,IX,IY) = G%K_TENSOR(IZ,IX,IY,3,3)


   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

ELSE !No orientation file

   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ,ISPEC) SHARED(G,SPROP,NCELLX,NCELLY,NCELLZ,NSSPEC) COLLAPSE(3)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
      G%KZ (IZ,IX,IY)  = 0D0
      G%KX (IZ,IX,IY)  = 0D0 
      G%KY (IZ,IX,IY)  = 0D0 
      IF (G%IMASK(IZ,IX,IY)) CYCLE
      DO ISPEC = 1, NSSPEC
         IF (G%XIN(ISPEC,IZ,IX,IY) .LT. EPSILON_FB) CYCLE
         G%KZ (IZ,IX,IY) = G%KZ (IZ,IX,IY) + G%XIN(ISPEC,IZ,IX,IY)*KZOFT(ISPEC,G%TP(IZ,IX,IY))
         G%KX (IZ,IX,IY) = G%KX (IZ,IX,IY) + G%XIN(ISPEC,IZ,IX,IY)*KXOFT(ISPEC,G%TP(IZ,IX,IY))
         G%KY (IZ,IX,IY) = G%KY (IZ,IX,IY) + G%XIN(ISPEC,IZ,IX,IY)*KYOFT(ISPEC,G%TP(IZ,IX,IY))
         IF (SPROP%GAMMA(ISPEC) .GT. 0D0) THEN 
            G%KZ(IZ,IX,IY) = G%KZ(IZ,IX,IY) + G%XIN(ISPEC,IZ,IX,IY)*SPROP%GAMMA(ISPEC)*5.67D-8*G%TP(IZ,IX,IY)*G%TP(IZ,IX,IY)*G%TP(IZ,IX,IY)
            G%KX(IZ,IX,IY) = G%KX(IZ,IX,IY) + G%XIN(ISPEC,IZ,IX,IY)*SPROP%GAMMA(ISPEC)*5.67D-8*G%TP(IZ,IX,IY)*G%TP(IZ,IX,IY)*G%TP(IZ,IX,IY)
            G%KY(IZ,IX,IY) = G%KY(IZ,IX,IY) + G%XIN(ISPEC,IZ,IX,IY)*SPROP%GAMMA(ISPEC)*5.67D-8*G%TP(IZ,IX,IY)*G%TP(IZ,IX,IY)*G%TP(IZ,IX,IY)
         ENDIF
      ENDDO
   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDIF


!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ,ISPEC) SHARED(G,SPROP,NCELLX,NCELLY,NCELLZ,NSSPEC) COLLAPSE(3)
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ
   G%KAPPA(IZ,IX,IY) = 0D0
   G%CPS  (IZ,IX,IY) = 0D0
   G%HI (:,IZ,IX,IY) = 0D0
   IF (G%IMASK(IZ,IX,IY)) CYCLE
   DO ISPEC = 1, NSSPEC
      ! Heat capacity
      G%CPS(IZ,IX,IY) = G%CPS(IZ,IX,IY) + G%YIN(ISPEC,IZ,IX,IY)*CPOFT(ISPEC,G%TP(IZ,IX,IY))
      ! Absorption
      G%KAPPA(IZ,IX,IY)=G%KAPPA(IZ,IX,IY)+G%XIN(ISPEC,IZ,IX,IY)*SPROP%KAPPA(ISPEC)
      ! Specific enthalpy
      G%HI(ISPEC,IZ,IX,IY) = HOFT(ISPEC,G%TP (IZ,IX,IY))
   ENDDO
ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO

! *****************************************************************************
END SUBROUTINE CALCULATE_WEIGHTED_POINT_THERMAL_PROPERTIES
! *****************************************************************************


! *****************************************************************************
SUBROUTINE CALCULATE_THERMAL_INTERFACE_QUANTITIES(IMESH,CQUANTITY,ISPEC)
! *****************************************************************************

INTEGER, INTENT(IN) :: IMESH,ISPEC
CHARACTER(9), INTENT(IN) :: CQUANTITY

REAL(EB), POINTER, DIMENSION(:,:,:) :: GAMPZ, GAMPX, GAMPY !"Gamma", from Patankar nomenclature 
REAL(EB), POINTER, DIMENSION(:,:,:) :: GAMT,GAMB,GAME,GAMW,GAMN,GAMS

INTEGER :: IZ,IX,IY,NCELLZ,NCELLX,NCELLY

REAL(EB), PARAMETER :: EPS = EPSILON_EB

G => GPM(IMESH)

NCELLZ = G%NCELLZ
NCELLX = G%NCELLX
NCELLY = G%NCELLY

GAMPZ => G%RWORK11
GAMPX => G%RWORK12
GAMPY => G%RWORK13

SELECT CASE(TRIM(CQUANTITY))
   CASE('KOC')

      IF (NCELLZ .GT. 1) THEN
         !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) SHARED(G,GAMPZ,NCELLX,NCELLY,NCELLZ) COLLAPSE(3)
         DO IY = 1, NCELLY
         DO IX = 1, NCELLX
         DO IZ = 1, NCELLZ
            IF (G%IMASK(IZ,IX,IY)) THEN; GAMPZ(IZ,IX,IY)=0D0; CYCLE; ENDIF
            GAMPZ(IZ,IX,IY) = G%KZ(IZ,IX,IY) / G%CPS(IZ,IX,IY)  !Gamma-z (point)
         ENDDO
         ENDDO
         ENDDO
         !$OMP END PARALLEL DO
         GAMT => G%KOCT
         GAMB => G%KOCB
      ENDIF

      IF (NCELLX .GT. 1) THEN
         !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) SHARED(G,GAMPX,NCELLX,NCELLY,NCELLZ) COLLAPSE(3)
         DO IY = 1, NCELLY
         DO IX = 1, NCELLX
         DO IZ = 1, NCELLZ
            IF (G%IMASK(IZ,IX,IY)) THEN; GAMPX(IZ,IX,IY)=0D0; CYCLE; ENDIF
            GAMPX(IZ,IX,IY) = G%KX(IZ,IX,IY) / G%CPS(IZ,IX,IY)  !Gamma-x (point)
         ENDDO
         ENDDO
         ENDDO
         !$OMP END PARALLEL DO
         GAME => G%KOCE
         GAMW => G%KOCW
      ENDIF

      IF (NCELLY .GT. 1) THEN
         !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) SHARED(G,GAMPY,NCELLX,NCELLY,NCELLZ) COLLAPSE(3)
         DO IY = 1, NCELLY
         DO IX = 1, NCELLX
         DO IZ = 1, NCELLZ
            IF (G%IMASK(IZ,IX,IY)) THEN; GAMPY(IZ,IX,IY)=0D0; CYCLE; ENDIF
            GAMPY(IZ,IX,IY) = G%KY(IZ,IX,IY) / G%CPS(IZ,IX,IY)  !Gamma-y (point)
         ENDDO
         ENDDO
         ENDDO
         !$OMP END PARALLEL DO
         GAMN => G%KOCN
         GAMS => G%KOCS
      ENDIF
  
   CASE('KOCHI')

      IF (NCELLZ .GT. 1) THEN
         !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) SHARED(G,GAMPZ,NCELLX,NCELLY,NCELLZ) COLLAPSE(3)
         DO IY = 1, NCELLY
         DO IX = 1, NCELLX
         DO IZ = 1, NCELLZ
            IF (G%IMASK(IZ,IX,IY)) THEN; GAMPZ(IZ,IX,IY)=0D0; CYCLE; ENDIF
            GAMPZ(IZ,IX,IY) = ((G%KZ(IZ,IX,IY)) / G%CPS(IZ,IX,IY)) * G%HI(ISPEC,IZ,IX,IY) !Gamma-z (point)
         ENDDO
         ENDDO
         ENDDO
         !$OMP END PARALLEL DO
         GAMT => G%KOCHIT(ISPEC,:,:,:)
         GAMB => G%KOCHIB(ISPEC,:,:,:)
      ENDIF

      IF (NCELLX .GT. 1) THEN
         !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) SHARED(G,GAMPX,NCELLX,NCELLY,NCELLZ) COLLAPSE(3)
         DO IY = 1, NCELLY
         DO IX = 1, NCELLX
         DO IZ = 1, NCELLZ
            IF (G%IMASK(IZ,IX,IY)) THEN; GAMPX(IZ,IX,IY)=0D0; CYCLE; ENDIF
            GAMPX(IZ,IX,IY) = ((G%KX(IZ,IX,IY)) / G%CPS(IZ,IX,IY)) * G%HI(ISPEC,IZ,IX,IY) !Gamma-x (point)
         ENDDO
         ENDDO
         ENDDO
         !$OMP END PARALLEL DO
         GAME => G%KOCHIE(ISPEC,:,:,:)
         GAMW => G%KOCHIW(ISPEC,:,:,:)
      ENDIF

      IF (NCELLY .GT. 1) THEN
         !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) SHARED(G,GAMPY,NCELLX,NCELLY,NCELLZ) COLLAPSE(3)
         DO IY = 1, NCELLY
         DO IX = 1, NCELLX
         DO IZ = 1, NCELLZ
            IF (G%IMASK(IZ,IX,IY)) THEN; GAMPY(IZ,IX,IY)=0D0; CYCLE; ENDIF
            GAMPY(IZ,IX,IY) = ((G%KY(IZ,IX,IY)) / G%CPS(IZ,IX,IY)) * G%HI(ISPEC,IZ,IX,IY) !Gamma-y (point)
         ENDDO
         ENDDO
         ENDDO
         !$OMP END PARALLEL DO
         GAMN => G%KOCHIN(ISPEC,:,:,:)
         GAMS => G%KOCHIS(ISPEC,:,:,:)
      ENDIF

   CASE DEFAULT

      CONTINUE
          
END SELECT

CALL CALCULATE_INTERFACE_QUANTITIES(IMESH,GAMPZ,GAMPX,GAMPY, &
                                    GAMT,GAMB,GAME,GAMW,GAMN,GAMS)

! *****************************************************************************
END SUBROUTINE CALCULATE_THERMAL_INTERFACE_QUANTITIES
! *****************************************************************************


! *****************************************************************************
SUBROUTINE HANDLE_SHYI_CORRECTION(IMESH,B)
! *****************************************************************************
INTEGER, INTENT(IN) :: IMESH
REAL(EB),INTENT(INOUT), DIMENSION(:,:,:) :: B
INTEGER :: NCELLX, NCELLY, NCELLZ,NSPEC
INTEGER :: IX,IY,IZ,ISPEC
REAL(EB) :: DYIB, DYIT, DYIE, DYIW, DYIN, DYIS
REAL(EB) :: ABHI, ATHI, AEHI, AWHI, ANHI, ASHI


G=>GPM(IMESH)

NCELLX = G%NCELLX
NCELLY = G%NCELLY
NCELLZ = G%NCELLZ
NSPEC  = SPROP%NSSPEC

!$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(3)
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ
   G%SHYI(IZ,IX,IY) = 0D0
ENDDO
ENDDO 
ENDDO 
!$OMP END PARALLEL DO


!Calculate coefficients for z-direction correction:
IF (NCELLZ .GT. 1) THEN
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ,ISPEC,DYIT,DYIB,ATHI,ABHI) &
   !$OMP SHARED(G,NSPEC,NCELLX,NCELLY,NCELLZ) COLLAPSE(3)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
   IF (G%IMASK(IZ,IX,IY)) CYCLE
   IF (.NOT. G%NEEDSBCT(IZ,IX,IY)) THEN
      DO ISPEC = 1, NSPEC
         DYIT = G%YIN(ISPEC,IZ,IX,IY)-G%YIN(ISPEC,IZ-1,IX,IY)
         IF (ABS(DYIT).LT. EPSILON_FB) CYCLE
         ATHI = G%KOCHIT(ISPEC,IZ,IX,IY)/G%DZT(IZ,IX,IY)
         G%SHYI(IZ,IX,IY)=G%SHYI(IZ,IX,IY)+ ATHI * DYIT * G%DXDY(IZ,IX,IY)
      ENDDO
   ENDIF
   IF (.NOT. G%NEEDSBCB(IZ,IX,IY)) THEN
      DO ISPEC = 1, NSPEC
         DYIB = G%YIN(ISPEC,IZ,IX,IY)-G%YIN(ISPEC,IZ+1,IX,IY)
         IF (ABS(DYIB).LT. EPSILON_FB) CYCLE
         ABHI = G%KOCHIB(ISPEC,IZ,IX,IY)/G%DZB(IZ,IX,IY)  
         G%SHYI(IZ,IX,IY)=G%SHYI(IZ,IX,IY) + ABHI * DYIB * G%DXDY(IZ,IX,IY)
      ENDDO
   ENDIF

   ENDDO
   ENDDO 
   ENDDO 
   !$OMP END PARALLEL DO
ENDIF 

IF (NCELLX .GT. 1) THEN
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ,ISPEC,DYIW,DYIE,AWHI,AEHI) &
   !$OMP SHARED(G,NSPEC,NCELLX,NCELLY,NCELLZ) COLLAPSE(3)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
      IF (G%IMASK(IZ,IX,IY)) CYCLE
      IF (.NOT. G%NEEDSBCW(IZ,IX,IY)) THEN
         DO ISPEC = 1, NSPEC
            DYIW=G%YIN(ISPEC,IZ,IX,IY)-G%YIN(ISPEC,IZ,IX-1,IY)
            IF (ABS(DYIW).LT. EPSILON_FB) CYCLE
            AWHI= G%KOCHIW(ISPEC,IZ,IX,IY)/G%DXW(IZ,IX,IY)
            G%SHYI(IZ,IX,IY)=G%SHYI(IZ,IX,IY)+ AWHI * DYIW * G%DYDZ(IZ,IX,IY)
         ENDDO
      ENDIF

      IF (.NOT. G%NEEDSBCE(IZ,IX,IY)) THEN
         DO ISPEC = 1, NSPEC
            DYIE= G%YIN(ISPEC,IZ,IX,IY)-G%YIN(ISPEC,IZ,IX+1,IY)
            IF (ABS(DYIE).LT. EPSILON_FB) CYCLE
            AEHI = G%KOCHIE(ISPEC,IZ,IX,IY)/G%DXE(IZ,IX,IY)
            G%SHYI(IZ,IX,IY)=G%SHYI(IZ,IX,IY) + AEHI * DYIE * G%DYDZ(IZ,IX,IY)
         ENDDO
      ENDIF
   ENDDO
   ENDDO 
   ENDDO 
   !$OMP END PARALLEL DO
ENDIF

IF (NCELLY .GT. 1 ) THEN
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ,ISPEC,DYIS,DYIN,ASHI,ANHI) &
   !$OMP SHARED(G,NSPEC,NCELLX,NCELLY,NCELLZ) COLLAPSE(3)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
      IF (G%IMASK(IZ,IX,IY)) CYCLE
      IF (.NOT. G%NEEDSBCS(IZ,IX,IY)) THEN
         DO ISPEC = 1, NSPEC
            DYIS=G%YIN(ISPEC,IZ,IX,IY)-G%YIN(ISPEC,IZ,IX,IY-1)
            IF (ABS(DYIS).LT. EPSILON_FB) CYCLE
            ASHI = G%KOCHIS(ISPEC,IZ,IX,IY)/G%DYS(IZ,IX,IY)
            G%SHYI(IZ,IX,IY)=G%SHYI(IZ,IX,IY) + ASHI * DYIS * G%DXDZ(IZ,IX,IY)
         ENDDO
      ENDIF

      IF (.NOT. G%NEEDSBCN(IZ,IX,IY)) THEN
         DO ISPEC = 1, NSPEC
            DYIN= G%YIN(ISPEC,IZ,IX,IY)-G%YIN(ISPEC,IZ,IX,IY+1)
            IF (ABS(DYIN).LT. EPSILON_FB) CYCLE
            ANHI = G%KOCHIN(ISPEC,IZ,IX,IY)/G%DYN(IZ,IX,IY)

            G%SHYI(IZ,IX,IY)=G%SHYI(IZ,IX,IY) + ANHI * DYIN * G%DXDZ(IZ,IX,IY)
         ENDDO
      ENDIF
   ENDDO
   ENDDO 
   ENDDO 
   !$OMP END PARALLEL DO
ENDIF


! Add SHYI_CORRECTION coefficients to the TDMA parrameter
!$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(3) DEFAULT(SHARED) PRIVATE(IX,IY,IZ) 
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ
   IF (G%IMASK(IZ,IX,IY)) CYCLE
   B (IZ,IX,IY) = B (IZ,IX,IY) + G%SHYI(IZ,IX,IY)
ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO


! *****************************************************************************
END SUBROUTINE HANDLE_SHYI_CORRECTION
! *****************************************************************************



! *****************************************************************************
SUBROUTINE SOLID_ENERGY_SOURCE_TERMS(IMESH,NCELLZ,NCELLX,NCELLY,DTIME)
! *****************************************************************************

INTEGER, INTENT(IN) :: IMESH,NCELLZ,NCELLX,NCELLY
REAL(EB), INTENT(IN) :: DTIME
INTEGER :: ISPEC,IZ,IX,IY,IRXN,IOR,NRXN
REAL(EB):: MDOTBACK(1:NCELLX,1:NCELLY)
REAL(EB), POINTER, DIMENSION(:,:,:) :: MDOTPPT, FLUX
REAL(EB):: DH

NRXN=SPROP%NRXN
! Calculate positive and negative components of source term
! for solid enthalpy equation:

G => GPM(IMESH)
GPBCP(1:,1:,1:,-3:)=>G%GPYRO_BOUNDARY_CONDITION(1:,1:,1:,-3:)

MDOTPPT=>G%RWORK01
FLUX   =>G%RWORK02


G%SHNET= 0D0

! Gas/solid heat transfer:
IF (GPG%SOLVE_GAS_ENERGY .AND. (.NOT. GPG%THERMAL_EQUILIBRIUM)) THEN
      !$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(3)
      DO IY = 1, NCELLY
      DO IX = 1, NCELLX
      DO IZ = 1, NCELLZ
         G%SHNET(IZ,IX,IY)= G%HCV(IZ,IX,IY)*(G%TGN(IZ,IX,IY)-G%TP(IZ,IX,IY))
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO

ELSEIF (GPG%SOLVE_GAS_ENERGY .AND. GPG%THERMAL_EQUILIBRIUM ) THEN
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
      IF (GPG%SOLVE_PRESSURE) THEN !Mass flux calcd at cell 'top':
         !$OMP PARALLEL DO SCHEDULE(STATIC)
         DO IZ = 1, NCELLZ
            MDOTPPT(IZ,IX,IY) =  G%MDOTPPDARCYT(IZ,IX,IY)
         ENDDO
         !$OMP END PARALLEL DO

      ELSE
         !$OMP PARALLEL DO SCHEDULE(STATIC)
         DO IZ = 1, NCELLZ
            MDOTPPT(IZ,IX,IY) = -G%MDOTPPZ(0,IZ,IX,IY)
         ENDDO
         !$OMP END PARALLEL DO
      ENDIF

      IF (GPG%FULL_QSG) THEN !Detailed way of calculating QSG
         !$OMP PARALLEL DO SCHEDULE(STATIC)
         DO IZ = 1, NCELLZ
            G%QSG(IZ,IX,IY) = G%RG(IZ,IX,IY) * G%POROSS(IZ,IX,IY) * GPROP%CPG * (G%TPN(IZ,IX,IY) - G%TP(IZ,IX,IY)) / DTIME
         ENDDO
         !$OMP END PARALLEL DO

         !$OMP PARALLEL DO SCHEDULE(STATIC)
         DO IZ = NCELLZ-1, 1, -1
            G%QSG(IZ,IX,IY)= G%QSG(IZ,IX,IY) + MDOTPPT(IZ+1,IX,IY) * GPROP%CPG * (G%TP(IZ+1,IX,IY)-G%TP(IZ,IX,IY)) / G%DZT(IZ,IX,IY)
         ENDDO
         !$OMP END PARALLEL DO

   
         FLUX(NCELLZ,IX,IY) = 0D0 !Defined at bottom of cell. This should be done more precisely, FIX!!!
         !$OMP PARALLEL DO SCHEDULE(STATIC)
         DO IZ = 1, NCELLZ-1
            FLUX(IZ,IX,IY) = G%PSIRGDB(IZ,IX,IY) * GPROP%CPG * (G%TP(IZ+1,IX,IY)-G%TP(IZ,IX,IY)) / G%DZT(IZ,IX,IY) 
         ENDDO
         !$OMP END PARALLEL DO

         !$OMP PARALLEL DO SCHEDULE(STATIC)
         DO IZ = 1, NCELLZ-1
            G%QSG(IZ,IX,IY) = G%QSG(IZ,IX,IY) - (FLUX(IZ+1,IX,IY) - FLUX (IZ,IX,IY)) / G%DZT(IZ,IX,IY)
         ENDDO
         !$OMP END PARALLEL DO

               
      ELSE !Default way of calculating QSG
         !$OMP PARALLEL DO SCHEDULE(STATIC)
         DO IZ = 1, NCELLZ
            G%QSG(IZ,IX,IY)= 0D0
         ENDDO
         !$OMP END PARALLEL DO

         !$OMP PARALLEL DO SCHEDULE(STATIC)
         DO IZ = NCELLZ-1, 1, -1
               G%QSG(IZ,IX,IY)= MDOTPPT(IZ+1,IX,IY) * GPROP%CPG * (G%TP(IZ+1,IX,IY)-G%TP(IZ,IX,IY)) / G%DZT(IZ,IX,IY)
         ENDDO
         !$OMP END PARALLEL DO

      ENDIF
   ENDDO      
   ENDDO 

   ! If thermal equilibrium, add heat release from homogeneous
   ! gas-phase reactions to condensed-phase
   IF (GPG%SOLVE_GAS_ENERGY .AND. GPROP%NHGRXN .GT. 0) THEN 
      !$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(4)
      DO IRXN = 1, GPROP%NHGRXN
      DO IY = 1, NCELLY
      DO IX = 1, NCELLX
      DO IZ = 1, NCELLZ
         G%QSG(IZ,IX,IY) = G%QSG(IZ,IX,IY) + G%HGRR(IRXN,IZ,IX,IY) * HGRXN(IRXN)%DH 
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   ENDIF

   ! If gases coming in back face are not in thermal equilibrium with 
   ! solid-phase, then QSG must be accounted for:
   !$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(2)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
      IOR = -3
      IF (GPBCP(NCELLZ,IX,IY,IOR)%PRES .LT. 0D0)  THEN
         MDOTBACK(IX,IY) = -1D-3 * GPBCP(NCELLZ,IX,IY,IOR)%MFLUX
         G%QSG(NCELLZ,IX,IY) = MDOTBACK(IX,IY) * GPROP%CPG &
                              * (GPBCP(NCELLZ,IX,IY,IOR)%TINF - G%TP(NCELLZ,IX,IY)) / G%DLTZN(NCELLZ,IX,IY)
      ENDIF
   ENDDO      
   ENDDO 
   !$OMP END PARALLEL DO

   !$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(2)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX   
   DO IZ = 1, NCELLZ
      G%SHNET(IZ,IX,IY)= G%SHNET(IZ,IX,IY) - G%QSG(IZ,IX,IY)
   ENDDO
   ENDDO
   ENDDO 
   !$OMP END PARALLEL DO

   
ELSE ! NOT GPG%SOLVE_GAS_ENERGY
   !$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(3)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
      G%QSG(IZ,IX,IY) = 0D0
   ENDDO 
   ENDDO      
   ENDDO      
   !$OMP END PARALLEL DO  
ENDIF


! This term shows up when the gas-phase mass conservation equation
! is multiplied by HPN and then subtracted from the energy
! conservation equation 

!G%SHNET(:,:,:) = G%SHNET(:,:,:) + G%OMEGASFG(:,:,:) * G%HP(:,:,:)

! Get enthalpy of each species and then do 
! formation and destruction of solid-phase species due to reactions
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IZ,IX,IY,ISPEC,IRXN) SHARED(NCELLY,NCELLX,NCELLZ,G) COLLAPSE(3)
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ
   IF (G%IMASK(IZ,IX,IY)) CYCLE
   G%SHNET(IZ,IX,IY) = G%SHNET(IZ,IX,IY) + G%OMEGASFG(IZ,IX,IY) * G%HP(IZ,IX,IY)
   DO ISPEC = 1, SPROP%NSSPEC 
      IF (.NOT. G%IS_REACTING(0,IZ,IX,IY)) CYCLE
      G%SHNET(IZ,IX,IY) = G%SHNET(IZ,IX,IY) +G%SOMEGA(3,ISPEC,IZ,IX,IY)* G%HI(ISPEC,IZ,IX,IY)
   ENDDO
ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO

DO IRXN = 1, NRXN
   IF (RXN(IRXN)%DHS .GT. -1D-6 .AND. RXN(IRXN)%DHS .LT. 1D-6) CYCLE
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IZ,IX,IY) SHARED(NCELLY,NCELLX,NCELLZ,G,IRXN,RXN) COLLAPSE(3)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
      IF (G%IMASK(IZ,IX,IY)) CYCLE
      IF (.NOT. G%IS_REACTING(IRXN,IZ,IX,IY)) CYCLE
      G%SHNET(IZ,IX,IY) = G%SHNET(IZ,IX,IY) - G%OMEGASFBK(IRXN,IZ,IX,IY) * RXN(IRXN)%DHS
   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDDO

DO IRXN = 1, NRXN
   IF (RXN(IRXN)%DHV .GT. -1D-6 .AND. RXN(IRXN)%DHV .LT. 1D-6) CYCLE
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IZ,IX,IY) SHARED(NCELLY,NCELLX,NCELLZ,G,IRXN,RXN) COLLAPSE(3)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
      IF ( G%IMASK(IZ,IX,IY)                ) CYCLE
      IF (.NOT. G%IS_REACTING(IRXN,IZ,IX,IY)) CYCLE
      G%SHNET(IZ,IX,IY) = G%SHNET(IZ,IX,IY) - G%OMEGASFGK(IRXN,IZ,IX,IY) * RXN(IRXN)%DHV
   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDDO

! Gas-phase heat of reaction:
IF (.NOT. GPG%CONSTANT_DHVOL) THEN
   IF (GPG%THERMAL_EQUILIBRIUM .OR. (.NOT. GPG%SOLVE_GAS_ENERGY)) THEN
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IRXN,IX,IY,IZ,DH) DEFAULT(SHARED) COLLAPSE(3)
      DO IY = 1, NCELLY
      DO IX = 1, NCELLX
      DO IZ = 1, NCELLZ
      DO IRXN = 1, NRXN
         IF (.NOT. G%IS_REACTING(IRXN,IZ,IX,IY)) CYCLE
         DH = GPROP%CPG*G%TP(IZ,IX,IY) - G%HI(RXN(IRXN)%IFROM,IZ,IX,IY)
         G%SHNET(IZ,IX,IY) = G%SHNET(IZ,IX,IY) - G%OMEGASFGK(IRXN,IZ,IX,IY) * DH
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO

   ELSE
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IRXN,IX,IY,IZ,DH) DEFAULT(SHARED) COLLAPSE(3)
      DO IY = 1, NCELLY
      DO IX = 1, NCELLX
      DO IZ = 1, NCELLZ
      DO IRXN = 1, NRXN
         IF (.NOT. G%IS_REACTING(IRXN,IZ,IX,IY)) CYCLE
         DH = GPROP%CPG*G%TGN(IZ,IX,IY) - G%HI(RXN(IRXN)%IFROM,IZ,IX,IY)
         G%SHNET(IZ,IX,IY) = G%SHNET(IZ,IX,IY) - G%OMEGASFGK(IRXN,IZ,IX,IY) * DH
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   ENDIF
ENDIF

! Volumetric heat losses:	
IF (GPG%VHLC .GT. 0D0) THEN
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) DEFAULT(SHARED) COLLAPSE(3)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
      IF (G%IMASK(IZ,IX,IY)) CYCLE
      G%SHNET(IZ,IX,IY) = G%SHNET(IZ,IX,IY) - GPG%VHLC*(G%TPN(IZ,IX,IY)-GPG%TAMB)
   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDIF
! *****************************************************************************
END SUBROUTINE SOLID_ENERGY_SOURCE_TERMS
! *****************************************************************************




! *****************************************************************************
SUBROUTINE GET_COEFF_OF_SOLID_ENTHALPY_SOLVER(IMESH,DTIME,B,AP,AB,AT,AE,AW,AN,AS)
! *****************************************************************************
INTEGER, INTENT(IN) :: IMESH
REAL(EB),INTENT(IN) :: DTIME
REAL(EB),INTENT(OUT), DIMENSION(:,:,:) :: AP,AB,AT,AE,AW,AN,AS,B
REAL(EB), POINTER, DIMENSION(:,:,:) :: AP0
INTEGER :: IZ, IX, IY, NCELLX, NCELLY, NCELLZ

NCELLX = G%NCELLX
NCELLY = G%NCELLY
NCELLZ = G%NCELLZ

AP0     => G%RWORK01

CALL CALCULATE_THERMAL_INTERFACE_QUANTITIES(IMESH,'KOC      ',0)

! Specify coefficients for TDMA:
!$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(3) DEFAULT(SHARED) PRIVATE(IX,IY,IZ) 
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ
   IF (G%IMASK(IZ,IX,IY)) THEN
      B  (IZ,IX,IY) = G%HP(IZ,IX,IY)
      AP (IZ,IX,IY) = 1D0
      IF (NCELLZ .GT. 1) THEN ; AT (IZ,IX,IY) = 0D0 ; AB (IZ,IX,IY) = 0D0 ; ENDIF
      IF (NCELLX .GT. 1) THEN ; AE (IZ,IX,IY) = 0D0 ; AW (IZ,IX,IY) = 0D0 ; ENDIF
      IF (NCELLZ .GT. 1) THEN ; AS (IZ,IX,IY) = 0D0 ; AN (IZ,IX,IY) = 0D0 ; ENDIF
      CYCLE
   ENDIF

   AP0(IZ,IX,IY) = G%RDLTZN(IZ,IX,IY) * G%DXDY(IZ,IX,IY) / DTIME
   AP (IZ,IX,IY) = AP0(IZ,IX,IY)
   B  (IZ,IX,IY) = AP0(IZ,IX,IY) * G%HP(IZ,IX,IY)
   IF (NCELLZ .GT. 1) THEN
      AB (IZ,IX,IY) = G%KOCB(IZ,IX,IY) * G%DXDY(IZ,IX,IY) / G%DZB(IZ,IX,IY)
      AT (IZ,IX,IY) = G%KOCT(IZ,IX,IY) * G%DXDY(IZ,IX,IY) / G%DZT(IZ,IX,IY)
   ENDIF
   IF (NCELLX .GT. 1) THEN
      AE (IZ,IX,IY) = G%KOCE(IZ,IX,IY) * G%DYDZ(IZ,IX,IY) / G%DXE(IZ,IX,IY)
      AW (IZ,IX,IY) = G%KOCW(IZ,IX,IY) * G%DYDZ(IZ,IX,IY) / G%DXW(IZ,IX,IY)
   ENDIF
   IF (NCELLY .GT. 1) THEN
      AN (IZ,IX,IY) = G%KOCN(IZ,IX,IY) * G%DXDZ(IZ,IX,IY) / G%DYN(IZ,IX,IY)
      AS (IZ,IX,IY) = G%KOCS(IZ,IX,IY) * G%DXDZ(IZ,IX,IY) / G%DYS(IZ,IX,IY)
   ENDIF
ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO
! *****************************************************************************
END SUBROUTINE GET_COEFF_OF_SOLID_ENTHALPY_SOLVER
! *****************************************************************************


! *****************************************************************************
 SUBROUTINE GET_COEFF_OF_SOLID_ENTHALPY_SOLVER_ANISOTROPIC(DTIME,B,AP,AB,AT,AE,AW,AN,AS)
! *****************************************************************************
REAL(EB),INTENT(IN) :: DTIME
REAL(EB),INTENT(OUT), DIMENSION(:,:,:) :: AP,AB,AT,AE,AW,AN,AS,B
REAL(EB), POINTER, DIMENSION(:,:,:) :: AP0
REAL(EB), POINTER, DIMENSION(:,:,:) :: A11W, A12W, A13W, A11E, A12E, A13E, &
                                       A21S, A22S, A23S, A21N, A22N, A23N, &
                                       A31T, A32T, A33T, A31B, A32B, A33B


INTEGER :: IZ, IX, IY, NCELLX, NCELLY, NCELLZ

NCELLX = G%NCELLX
NCELLY = G%NCELLY
NCELLZ = G%NCELLZ

AP0     => G%RWORK01

A11W    => G%RWORK19
A12W    => G%RWORK20
A13W    => G%RWORK21
A11E    => G%RWORK22
A12E    => G%RWORK23
A13E    => G%RWORK24

A21S    => G%RWORK25
A22S    => G%RWORK26
A23S    => G%RWORK27
A21N    => G%RWORK28
A22N    => G%RWORK29
A23N    => G%RWORK30

A31T    => G%RWORK31
A32T    => G%RWORK32
A33T    => G%RWORK33
A31B    => G%RWORK34
A32B    => G%RWORK35
A33B    => G%RWORK36



! Specify coefficients for TDMA:
!$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(3) DEFAULT(SHARED) PRIVATE(IX,IY,IZ) 
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ
   IF (G%IMASK(IZ,IX,IY)) THEN
      B  (IZ,IX,IY) = G%HP(IZ,IX,IY)
      AP (IZ,IX,IY) = 1D0
      AT (IZ,IX,IY) = 0D0 ; AB (IZ,IX,IY) = 0D0
      AE (IZ,IX,IY) = 0D0 ; AW (IZ,IX,IY) = 0D0
      AS (IZ,IX,IY) = 0D0 ; AN (IZ,IX,IY) = 0D0
      CYCLE
   ENDIF
   A11W(IZ,IX,IY) = G%K_TENSOR(IZ,IX,IY,1,1) * G%DYDZ(IZ,IX,IY) / ( G%DXW(IZ,IX,IY) * G%CPS(IZ,IX,IY) )
   A12W(IZ,IX,IY) = G%K_TENSOR(IZ,IX,IY,1,2) * G%DYDZ(IZ,IX,IY) / ( G%DYS(IZ,IX,IY) * G%CPS(IZ,IX,IY) )
   A13W(IZ,IX,IY) = G%K_TENSOR(IZ,IX,IY,1,3) * G%DYDZ(IZ,IX,IY) / ( G%DZT(IZ,IX,IY) * G%CPS(IZ,IX,IY) )

   A11E(IZ,IX,IY) = A11W(IZ,IX,IY)
   A12E(IZ,IX,IY) = A12W(IZ,IX,IY)
   A13E(IZ,IX,IY) = A13W(IZ,IX,IY)

   A21S(IZ,IX,IY) = G%K_TENSOR(IZ,IX,IY,2,1) * G%DXDZ(IZ,IX,IY) / ( G%DXW(IZ,IX,IY) * G%CPS(IZ,IX,IY) )
   A22S(IZ,IX,IY) = G%K_TENSOR(IZ,IX,IY,2,2) * G%DXDZ(IZ,IX,IY) / ( G%DYS(IZ,IX,IY) * G%CPS(IZ,IX,IY) )
   A23S(IZ,IX,IY) = G%K_TENSOR(IZ,IX,IY,2,3) * G%DXDZ(IZ,IX,IY) / ( G%DZT(IZ,IX,IY) * G%CPS(IZ,IX,IY) )

   A21N(IZ,IX,IY) = A21S(IZ,IX,IY)
   A22N(IZ,IX,IY) = A22S(IZ,IX,IY)
   A23N(IZ,IX,IY) = A23S(IZ,IX,IY)

   A31T(IZ,IX,IY) = G%K_TENSOR(IZ,IX,IY,3,1) * G%DXDY(IZ,IX,IY) / ( G%DXW(IZ,IX,IY) * G%CPS(IZ,IX,IY) )
   A32T(IZ,IX,IY) = G%K_TENSOR(IZ,IX,IY,3,2) * G%DXDY(IZ,IX,IY) / ( G%DYS(IZ,IX,IY) * G%CPS(IZ,IX,IY) )
   A33T(IZ,IX,IY) = G%K_TENSOR(IZ,IX,IY,3,3) * G%DXDY(IZ,IX,IY) / ( G%DZT(IZ,IX,IY) * G%CPS(IZ,IX,IY) )

   A31B(IZ,IX,IY) = A31T(IZ,IX,IY)
   A32B(IZ,IX,IY) = A32T(IZ,IX,IY)
   A33B(IZ,IX,IY) = A33T(IZ,IX,IY)

   AW(IZ,IX,IY) = A11W(IZ,IX,IY) + A21S(IZ,IX,IY) + A31T(IZ,IX,IY)
   AS(IZ,IX,IY) = A12W(IZ,IX,IY) + A22S(IZ,IX,IY) + A32T(IZ,IX,IY)
   AT(IZ,IX,IY) = A13W(IZ,IX,IY) + A23S(IZ,IX,IY) + A33T(IZ,IX,IY)

   AE(IZ,IX,IY) = A11E(IZ,IX,IY) + A21N(IZ,IX,IY) + A31B(IZ,IX,IY)
   AN(IZ,IX,IY) = A12E(IZ,IX,IY) + A22N(IZ,IX,IY) + A32B(IZ,IX,IY)
   AB(IZ,IX,IY) = A13E(IZ,IX,IY) + A23N(IZ,IX,IY) + A33B(IZ,IX,IY)

   IF (G%NEEDSBCE(IZ,IX,IY)) AE(IZ,IX,IY) = 0D0
   IF (G%NEEDSBCW(IZ,IX,IY)) AW(IZ,IX,IY) = 0D0
   IF (G%NEEDSBCS(IZ,IX,IY)) AS(IZ,IX,IY) = 0D0
   IF (G%NEEDSBCN(IZ,IX,IY)) AN(IZ,IX,IY) = 0D0
   IF (G%NEEDSBCT(IZ,IX,IY)) AT(IZ,IX,IY) = 0D0
   IF (G%NEEDSBCB(IZ,IX,IY)) AB(IZ,IX,IY) = 0D0

   AP0(IZ,IX,IY) = G%RDLTZN(IZ,IX,IY) * G%DXDY(IZ,IX,IY) / DTIME
   AP (IZ,IX,IY) = AP0(IZ,IX,IY)
   B  (IZ,IX,IY) = AP0(IZ,IX,IY) * G%HP(IZ,IX,IY)

ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO
! *****************************************************************************
END SUBROUTINE GET_COEFF_OF_SOLID_ENTHALPY_SOLVER_ANISOTROPIC
! *****************************************************************************


! *****************************************************************************
SUBROUTINE BC_FOR_SOLID_ENTHALPY_SOLVER(IMESH,SHNET,AP,AB,AT,AE,AW,AN,AS,B)
! *****************************************************************************
INTEGER, INTENT(IN   ) :: IMESH
REAL(EB),INTENT(INOUT), DIMENSION(:,:,:) :: SHNET
REAL(EB),INTENT(INOUT), DIMENSION(:,:,:) :: AP,AB,AT,AE,AW,AN,AS,B
INTEGER  :: IX, IY, IZ, IOR, ICOUNT
INTEGER  :: NCELLX, NCELLY, NCELLZ
REAL(EB), DIMENSION(:,:,:), POINTER :: SHP_THREAD

NCELLY = G%NCELLY
NCELLX = G%NCELLX
NCELLZ = G%NCELLZ

GPBCP(1:,1:,1:,-3:)=>G%GPYRO_BOUNDARY_CONDITION(1:,1:,1:,-3:)

! Boundary conditions are considered as sources term, and is added to SHNET,
! For the coefficient of TDMA solver, A_IOR = 0 in the direction of the boundary.

!$OMP PARALLEL & 
!$OMP PRIVATE(ICOUNT,IZ, IX, IY, IOR, SHP_THREAD) & 
!$OMP SHARED(NGPYRO_FACES_NEEDING_BCS,GP_BOUDARYS,G, GPBCP, AP, AT, AB, AE, AW, AN, AS, B,IMESH) &
!$OMP REDUCTION(+:SHNET)
ALLOCATE(SHP_THREAD(1:NCELLZ,1:NCELLX,1:NCELLY)); SHP_THREAD = 0D0
!$OMP DO SCHEDULE(STATIC)
DO ICOUNT = 1, NGPYRO_FACES_NEEDING_BCS(IMESH)
   IZ  = GP_BOUDARYS(IMESH)%IZ_GPYRO (ICOUNT)
   IX  = GP_BOUDARYS(IMESH)%IX_GPYRO (ICOUNT)
   IY  = GP_BOUDARYS(IMESH)%IY_GPYRO (ICOUNT)
   IOR = GP_BOUDARYS(IMESH)%IOR_GPYRO (ICOUNT)

   IF(GP_BOUDARYS(IMESH)%COMPLETE_CELL_AT_BC(ICOUNT)) THEN
      !WRITE(0,*) "IT IS A FULL CELL !!!!", IZ
      CALL BC_FOR_FULL_CELLS(IMESH,IZ,IX,IY,IOR,SHP_THREAD,AP,AT,AB,AE,AW,AN,AS,B)
   ELSE
      CALL BC_FOR_HALF_CELLS(IMESH,IZ,IX,IY,IOR,SHP_THREAD,SHNET,AP,AT,AB,AE,AW,AN,AS,B)
   ENDIF
END DO
!$OMP END DO

!$OMP CRITICAL
SHNET = SHNET + SHP_THREAD
DEALLOCATE(SHP_THREAD)
!$OMP END CRITICAL
!$OMP END PARALLEL

! *****************************************************************************
END SUBROUTINE BC_FOR_SOLID_ENTHALPY_SOLVER
! *****************************************************************************


! *****************************************************************************
SUBROUTINE BC_FOR_FULL_CELLS(IMESH,IZ,IX,IY,IOR,SHP_THREAD,AP,AT,AB,AE,AW,AN,AS,B)
! *****************************************************************************
REAL(EB), INTENT(INOUT), DIMENSION(:,:,:), POINTER :: SHP_THREAD
REAL(EB), INTENT(INOUT), DIMENSION(:,:,:),TARGET :: AP, AB, AT, AE, AW, AN, AS, B
INTEGER, INTENT(IN) :: IMESH,IZ,IX,IY,IOR

INTEGER :: IX_G, IY_G, IZ_G
INTEGER :: ISPEC, NSSPEC
REAL(EB) :: QCONV1, EMIS, KAPPA, SUMYIC0I,D_T4
REAL(EB), POINTER :: A_IOR_GOHST
REAL(EB)          :: A_IOR_BC
REAL(EB)          :: DS_IOR, DELTA_IOR
REAL(EB)          :: HTC_BC


SELECT CASE(IOR)
   CASE( 3)
      IZ_G        =  IZ-1
      IX_G        =  IX 
      IY_G        =  IY
      DELTA_IOR   =  G%DLTZN(IZ  ,IX  ,IY  )
      DS_IOR      =  G%DXDY (IZ  ,IX  ,IY  )
      A_IOR_BC    =  AT     (IZ  ,IX  ,IY  )
      A_IOR_GOHST => AB     (IZ_G,IX_G,IY_G)
   CASE(-3)
      IZ_G        =  IZ+1
      IX_G        =  IX
      IY_G        =  IY
      DELTA_IOR   =  G%DLTZN(IZ  ,IX  ,IY  )
      DS_IOR      =  G%DXDY (IZ  ,IX  ,IY  )
      A_IOR_BC    =  AB     (IZ  ,IX  ,IY  )
      A_IOR_GOHST => AT     (IZ_G,IX_G,IY_G)
   CASE( 1)
      IZ_G        =  IZ
      IX_G        =  IX+1
      IY_G        =  IY
      DELTA_IOR   =  G%DLTX (IZ  ,IX  ,IY  )
      DS_IOR      =  G%DYDZ (IZ  ,IX  ,IY  )
      A_IOR_BC    =  AE     (IZ  ,IX  ,IY  )
      A_IOR_GOHST => AW     (IZ_G,IX_G,IY_G)
   CASE(-1)
      IZ_G        =  IZ
      IX_G        =  IX-1
      IY_G        =  IY
      DELTA_IOR   =  G%DLTX (IZ  ,IX  ,IY  )
      DS_IOR      =  G%DYDZ (IZ  ,IX  ,IY  )
      A_IOR_BC    =  AW     (IZ  ,IX  ,IY  )
      A_IOR_GOHST => AE     (IZ_G,IX_G,IY_G)
   CASE( 2)
      IZ_G        =  IZ
      IX_G        =  IX
      IY_G        =  IY+1
      DELTA_IOR   =  G%DLTY (IZ  ,IX  ,IY  )
      DS_IOR      =  G%DXDZ (IZ  ,IX  ,IY  )
      A_IOR_BC    =  AN     (IZ  ,IX  ,IY  )
      A_IOR_GOHST => AS     (IZ_G,IX_G,IY_G)
   CASE(-2)
      IZ_G        =  IZ
      IX_G        =  IX
      IY_G        =  IY-1
      DELTA_IOR   =  G%DLTY (IZ  ,IX  ,IY  )
      DS_IOR      =  G%DXDZ (IZ  ,IX  ,IY  )
      A_IOR_BC    =  AS     (IZ  ,IX  ,IY  )
      A_IOR_GOHST => AN     (IZ_G,IX_G,IY_G)
END SELECT
NSSPEC = SPROP%NSSPEC


!=========================================!
!====== Case with fixed Temperature ======!
!=========================================!

IF (GPBCP(IZ,IX,IY,IOR)%TFIXED .GE. 0D0) THEN
   GPBCP(IZ,IX,IY,IOR)%HFIXED = 0D0
   DO ISPEC = 1, NSSPEC
      GPBCP(IZ,IX,IY,IOR)%HFIXED = GPBCP(IZ,IX,IY,IOR)%HFIXED + &
                                 G%YIN(ISPEC,IZ,IX,IY) * HOFT(ISPEC,GPBCP(IZ,IX,IY,IOR)%TFIXED)
   ENDDO
   A_IOR_GOHST        = -1D0
   AP(IZ_G,IX_G,IY_G) =  1D0
   B (IZ_G,IX_G,IY_G) = 2 * GPBCP(IZ,IX,IY,IOR)%HFIXED
   RETURN
ENDIF


!=========================================!
!============ BC WITH FLUX ===============!
!=========================================!

! HTC_BC (Heat Transfer Coefficient introduce by Boundary Condition).
HTC_BC = 0D0
! Clean the default value of B of IMASK
B (IZ_G,IX_G,IY_G) = 0D0
GPBCP(IZ,IX,IY,IOR)%QENET = 0D0
!=========== Convection Part =============!

QCONV1 = GPBCP(IZ,IX,IY,IOR)%HC0 *(GPBCP(IZ,IX,IY,IOR)%TINF - GPG%TDATUM)
!QCONV1 = GPBCP(IZ,IX,IY,IOR)%HC0 *(GPBCP(IZ,IX,IY,IOR)%TINF - GPBCP(IZ,IX,IY,IOR)%T_SURFACE_OLD)
B (IZ_G,IX_G,IY_G) =B (IZ_G,IX_G,IY_G) + QCONV1 * DS_IOR

HTC_BC = HTC_BC +GPBCP(IZ,IX,IY,IOR)%HC0


!=========== Radiation Part =============!
EMIS  = CALC_EMIS(IZ,IX,IY)

! Radiation emitted by the surface
IF (GPBCP(IZ,IX,IY,IOR)%RERAD) THEN 
   D_T4=(GPBCP(IZ,IX,IY,IOR)%T_SURFACE_OLD**4-GPBCP(IZ,IX,IY,IOR)%TINF**4)
   !D_T4 = G%TP(IZ,IX,IY)**4-GPBCP(IZ,IX,IY,IOR)%TINF**4
   GPBCP(IZ,IX,IY,IOR)%QENET = GPBCP(IZ,IX,IY,IOR)%QENET - EMIS * SIGMA * D_T4
ENDIF


! Incoming Radiation Flux
GPBCP(IZ,IX,IY,IOR)%QENET = GPBCP(IZ,IX,IY,IOR)%QENET + EMIS*GPBCP(IZ,IX,IY,IOR)%QE

KAPPA = G%KAPPA(IZ,IX,IY) !kappa at surface cell

IF ((GPBCP(IZ,IX,IY,IOR)%QENET .GT. 0D0) .AND. (KAPPA .LT. 1D6)) THEN !In depth Radiation
   CALL IN_DEPTH_RADIATION_ABSORPTION(GPBCP(IZ,IX,IY,IOR)%QENET,SHP_THREAD,IZ,IX,IY,IOR,IMESH)

ELSE
   B(IZ_G,IX_G,IY_G) = B(IZ_G,IX_G,IY_G) + GPBCP(IZ,IX,IY,IOR)%QENET * DS_IOR
ENDIF


!======== Change TDMA coeff ==========!

SUMYIC0I = 0D0
DO ISPEC = 1, NSSPEC
   SUMYIC0I = SUMYIC0I + G%YIN(ISPEC,IZ,IX,IY)*SPROP%C0(ISPEC)
ENDDO

HTC_BC = HTC_BC * DS_IOR/ (2*SUMYIC0I)

A_IOR_GOHST        = A_IOR_BC - HTC_BC
AP(IZ_G,IX_G,IY_G) = A_IOR_BC + HTC_BC

! *****************************************************************************
END SUBROUTINE BC_FOR_FULL_CELLS
! *****************************************************************************


! *****************************************************************************
SUBROUTINE BC_FOR_HALF_CELLS(IMESH,IZ,IX,IY,IOR,SHP_THREAD,SHNET,AP,AT,AB,AE,AW,AN,AS,B)
! *****************************************************************************
REAL(EB), INTENT(INOUT), DIMENSION(:,:,:), POINTER :: SHP_THREAD
REAL(EB), INTENT(INOUT), DIMENSION(:,:,:) :: SHNET
REAL(EB), INTENT(INOUT), DIMENSION(:,:,:) :: AP, AB, AT, AE, AW, AN, AS, B
INTEGER, INTENT(IN) :: IMESH,IZ,IX,IY,IOR

INTEGER :: ISPEC, NSSPEC
REAL(EB) :: QCONV, EMIS, KAPPA
REAL(EB) :: DELTA_IOR

NSSPEC = SPROP%NSSPEC

IF (ABS(IOR) == 3) DELTA_IOR = G%DLTZN(IZ,IX,IY)
IF (ABS(IOR) == 1) DELTA_IOR = G%DLTX (IZ,IX,IY)
IF (ABS(IOR) == 2) DELTA_IOR = G%DLTY (IZ,IX,IY)


!====== Case with fixed Temperature ======!
IF (GPBCP(IZ,IX,IY,IOR)%TFIXED .GE. 0D0) THEN
   GPBCP(IZ,IX,IY,IOR)%HFIXED = 0D0
   DO ISPEC = 1, NSSPEC
      GPBCP(IZ,IX,IY,IOR)%HFIXED = GPBCP(IZ,IX,IY,IOR)%HFIXED + &
                                 G%YIN(ISPEC,IZ,IX,IY) * HOFT(ISPEC,GPBCP(IZ,IX,IY,IOR)%TFIXED)
   ENDDO
   CALL SET_BC_FIXED_VALUE(AP(IZ,IX,IY),AT(IZ,IX,IY),AB(IZ,IX,IY),AE(IZ,IX,IY),AW(IZ,IX,IY),&
                           AN(IZ,IX,IY),AS(IZ,IX,IY),B(IZ,IX,IY),GPBCP(IZ,IX,IY,IOR)%HFIXED)
   RETURN
ENDIF


!========== BC WITH FLUX =================!

   ! Convection Part
   QCONV = GPBCP(IZ,IX,IY,IOR)%HC0 * (GPBCP(IZ,IX,IY,IOR)%TINF - G%TP(IZ,IX,IY))
   SHNET(IZ,IX,IY) = SHNET(IZ,IX,IY) + QCONV / DELTA_IOR

   ! Radiation Part
   EMIS  = CALC_EMIS(IZ,IX,IY)
   GPBCP(IZ,IX,IY,IOR)%QENET = EMIS * GPBCP(IZ,IX,IY,IOR)%QE

IF (GPBCP(IZ,IX,IY,IOR)%RERAD) THEN 
   GPBCP(IZ,IX,IY,IOR)%QENET = GPBCP(IZ,IX,IY,IOR)%QENET &
                           - EMIS*5.67D-8 *(G%TP(IZ,IX,IY)**4D0-GPBCP(IZ,IX,IY,IOR)%TINF**4D0)
ENDIF

   KAPPA = G%KAPPA(IZ,IX,IY) !kappa at surface cell
   IF ((GPBCP(IZ,IX,IY,IOR)%QENET .GT. 0.) .AND. (KAPPA .LT. 1D6)) THEN !In depth Radiation
      CALL IN_DEPTH_RADIATION_ABSORPTION(GPBCP(IZ,IX,IY,IOR)%QENET,SHP_THREAD,IZ,IX,IY,IOR,IMESH)
   ELSE! No In_depth_radiation
      SHNET(IZ,IX,IY) = SHNET(IZ,IX,IY) + GPBCP(IZ,IX,IY,IOR)%QENET / DELTA_IOR
   ENDIF

! Set A_IOR to zero
SELECT CASE (IOR) 

CASE ( 3) ; AT(IZ,IX,IY)= 0D0
CASE (-3) ; AB(IZ,IX,IY)= 0D0

CASE ( 2) ; AN(IZ,IX,IY)= 0D0
CASE (-2) ; AS(IZ,IX,IY)= 0D0

CASE ( 1) ; AE(IZ,IX,IY)= 0D0
CASE (-1) ; AW(IZ,IX,IY)= 0D0
END SELECT


! *****************************************************************************
END SUBROUTINE BC_FOR_HALF_CELLS
! *****************************************************************************

! *****************************************************************************
SUBROUTINE  HANDLE_POSITIVE_AND_NEGATIVE_SOURCE_TERM(SHNET,B)!,AP)
! *****************************************************************************
REAL(EB),INTENT(IN), DIMENSION(:,:,:) :: SHNET
REAL(EB),INTENT(INOUT), DIMENSION(:,:,:) :: B !, AP
INTEGER :: IX, IY, IZ
INTEGER :: NCELLX, NCELLY, NCELLZ

NCELLY = G%NCELLY
NCELLX = G%NCELLX
NCELLZ = G%NCELLZ

!To prevent divergence, add positive Source Terme SHP to B,
! and negative Source SHM Terme to Ap

!$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(3) DEFAULT(SHARED) PRIVATE(IX,IY,IZ)
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ
   IF (G%IMASK(IZ,IX,IY)) CYCLE
   B (IZ,IX,IY) = B (IZ,IX,IY) + SHNET(IZ,IX,IY)*G%DV(IZ,IX,IY)
   IF (SHNET(IZ,IX,IY) .GE. 0D0) THEN
      !B (IZ,IX,IY) = B (IZ,IX,IY) + SHNET(IZ,IX,IY)*G%DV(IZ,IX,IY)
      G%SHP(IZ,IX,IY) = SHNET(IZ,IX,IY)
   ELSE
      !write(0,*) "I am here", IZ
      !AP(IZ,IX,IY) = AP(IZ,IX,IY) - SHNET(IZ,IX,IY)*G%DV(IZ,IX,IY) / G%HPN(IZ,IX,IY)
      G%SHM(IZ,IX,IY) = -SHNET(IZ,IX,IY)
   ENDIF
ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO
! *****************************************************************************
END SUBROUTINE  HANDLE_POSITIVE_AND_NEGATIVE_SOURCE_TERM
! *****************************************************************************
      

! *****************************************************************************
SUBROUTINE SOLVE_TDMA_ACCORDING_SWEEP_DIRECTION(AP,AE,AW,AS,AN,AT,AB,B,PTR)
! *****************************************************************************
REAL(EB), INTENT(OUT), DIMENSION(:,:,:) :: PTR
REAL(EB), INTENT(IN ), DIMENSION(:,:,:) :: AP, AE, AW, AS, AN, AT, AB, B
REAL(EB), POINTER, DIMENSION(:,:,:) :: APSTART,BSTART, BX, BY, BZ

INTEGER :: IX, IY, IZ
INTEGER :: NCELLX, NCELLY, NCELLZ
REAL(EB):: HP_IJK
LOGICAL :: SWEEPZ,SWEEPX,SWEEPY
!REAL(EB) :: FLUX, DH


BZ      => G%RWORK13
BX      => G%RWORK14
BY      => G%RWORK15
BSTART  => G%RWORK16
APSTART => G%RWORK17

NCELLZ = G%NCELLZ
NCELLX = G%NCELLX
NCELLY = G%NCELLY


! Vary sweep direction across iterations 
CALL SET_SWEEP_DIRECTION(SWEEPZ,SWEEPX,SWEEPY)

! SOLVE TDMA according sweep direction
IF (SWEEPX) THEN
   !$OMP PARALLEL DO PRIVATE(IX,IY,IZ) DEFAULT(SHARED) SCHEDULE(STATIC) COLLAPSE(2)
   DO IY = 1, NCELLY
   DO IZ = 1, NCELLZ
      DO IX = 1, NCELLX
         BZ(IZ,IX,IY) = 0D0; BY(IZ,IX,IY) = 0D0; BSTART(IZ,IX,IY) = 0D0
         APSTART (IZ,IX,IY) = AP(IZ,IX,IY) + AE(IZ,IX,IY) + AW(IZ,IX,IY)

         HP_IJK = G%HP(IZ,IX,IY)
         
         IF (IZ .NE. 1     ) BZ(IZ,IX,IY) = BZ(IZ,IX,IY) + AT(IZ,IX,IY)*(G%HP(IZ-1,IX,IY)-HP_IJK)
         IF (IZ .NE. NCELLZ) BZ(IZ,IX,IY) = BZ(IZ,IX,IY) + AB(IZ,IX,IY)*(G%HP(IZ+1,IX,IY)-HP_IJK)

         IF (IY .NE. 1     ) BY(IZ,IX,IY) = BY(IZ,IX,IY) + AS(IZ,IX,IY)*(G%HP(IZ,IX,IY-1)-HP_IJK)
         IF (IY .NE. NCELLY) BY(IZ,IX,IY) = BY(IZ,IX,IY) + AN(IZ,IX,IY)*(G%HP(IZ,IX,IY+1)-HP_IJK)
         
         APSTART(IZ,IX,IY) = APSTART(IZ,IX,IY)
         BSTART (IZ,IX,IY) = B (IZ,IX,IY) + BY (IZ,IX,IY) + BZ (IZ,IX,IY)
      ENDDO
      CALL TDMA_SOLVER_GENERAL(NCELLX,PTR(IZ,:,IY),APSTART(IZ,:,IY),AE(IZ,:,IY),AW(IZ,:,IY),BSTART(IZ,:,IY))
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDIF

IF (SWEEPY) THEN
   !$OMP PARALLEL DO PRIVATE(IX,IY,IZ) DEFAULT(SHARED) SCHEDULE(STATIC) COLLAPSE(2)
   DO IZ = 1, NCELLZ
   DO IX = 1, NCELLX
      DO IY = 1, NCELLY
         BZ(IZ,IX,IY) = 0D0; BX(IZ,IX,IY) = 0D0; BSTART(IZ,IX,IY) = 0D0
         APSTART (IZ,IX,IY)= AP(IZ,IX,IY) + AN(IZ,IX,IY) + AS(IZ,IX,IY)

         HP_IJK = G%HP(IZ,IX,IY)

         IF (IZ .NE. 1     ) BZ(IZ,IX,IY) = BZ(IZ,IX,IY) + AT(IZ,IX,IY)*(G%HP(IZ-1,IX,IY)-HP_IJK)
         IF (IZ .NE. NCELLZ) BZ(IZ,IX,IY) = BZ(IZ,IX,IY) + AB(IZ,IX,IY)*(G%HP(IZ+1,IX,IY)-HP_IJK)

         IF (IX .NE. 1     ) BX(IZ,IX,IY) = BX(IZ,IX,IY) + AW(IZ,IX,IY)*(G%HP(IZ,IX-1,IY)-HP_IJK)
         IF (IX .NE. NCELLX) BX(IZ,IX,IY) = BX(IZ,IX,IY) + AE(IZ,IX,IY)*(G%HP(IZ,IX+1,IY)-HP_IJK)

         BSTART (IZ,IX,IY) = B (IZ,IX,IY) + BX (IZ,IX,IY) + BZ (IZ,IX,IY)
      ENDDO
      CALL TDMA_SOLVER_GENERAL(NCELLY,PTR(IZ,IX,:),APSTART(IZ,IX,:),AN(IZ,IX,:),AS(IZ,IX,:),BSTART(IZ,IX,:))
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

ENDIF

IF (SWEEPZ) THEN
   !$OMP PARALLEL DO PRIVATE(IX,IY,IZ) DEFAULT(SHARED) SCHEDULE(STATIC) COLLAPSE(2)
   DO IX = 1, NCELLX
   DO IY = 1, NCELLY
      DO IZ = 1, NCELLZ
         BY(IZ,IX,IY) = 0D0; BX(IZ,IX,IY) = 0D0; BSTART(IZ,IX,IY) = 0D0
         IF (G%IMASK(IZ,IX,IY)) THEN
            ! This is for gohst cell
            APSTART (IZ,IX,IY)= AP(IZ,IX,IY)
         ELSE
            APSTART (IZ,IX,IY)= AP(IZ,IX,IY) + AB(IZ,IX,IY) + AT(IZ,IX,IY)
         ENDIF

         HP_IJK = G%HP(IZ,IX,IY)

         IF (IX .NE. 1     ) BX(IZ,IX,IY) = BX(IZ,IX,IY) + AW(IZ,IX,IY)*(G%HP(IZ,IX-1,IY)-HP_IJK)
         IF (IX .NE. NCELLX) BX(IZ,IX,IY) = BX(IZ,IX,IY) + AE(IZ,IX,IY)*(G%HP(IZ,IX+1,IY)-HP_IJK)

         IF (IY .NE. 1     ) BY(IZ,IX,IY) = BY(IZ,IX,IY) + AS(IZ,IX,IY)*(G%HP(IZ,IX,IY-1)-HP_IJK)
         IF (IY .NE. NCELLY) BY(IZ,IX,IY) = BY(IZ,IX,IY) + AN(IZ,IX,IY)*(G%HP(IZ,IX,IY+1)-HP_IJK)

         BSTART (IZ,IX,IY) = B (IZ,IX,IY) + BX (IZ,IX,IY) + BY (IZ,IX,IY)

         IF (G%NEEDSBCT(IZ,IX,IY)) THEN
            ! WRITE(0,*) ""
            ! WRITE(0,*) ""
            ! WRITE(0,*) "FOR IZ=", IZ,", AP=", APSTART(IZ,IX,IY), "B=",B(IZ,IX,IY)
            ! WRITE(0,*) "             AT=", AT(IZ,IX,IY) , "AB=", AB(IZ,IX,IY)
            ! WRITE(0,*) ""
            ! WRITE(0,*) "FOR IZ=", IZ-1,", AP=", APSTART(IZ-1,IX,IY), "B=",B(IZ-1,IX,IY)
            ! WRITE(0,*) "             AT=", AT(IZ-1,IX,IY) , "AB=", AB(IZ-1,IX,IY)
            
         ENDIF
      ENDDO
      CALL TDMA_SOLVER_GENERAL(NCELLZ,PTR(:,IX,IY),APSTART(:,IX,IY),AB(:,IX,IY),AT(:,IX,IY),BSTART(:,IX,IY))

      DO IZ = 1, NCELLZ
         IF (G%NEEDSBCT(IZ,IX,IY)) THEN
            !DH = PTR(IZ-1,IX,IY) -PTR(IZ,IX,IY)
            !FLUX= DH * G%KOCT(IZ,IX,IY) * G%DXDY(IZ,IX,IY) / G%DZT(IZ,IX,IY)
            !WRITE(0,*) "IZ=", IZ, "FLUX= ", FLUX
         ENDIF
      ENDDO

   ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDIF

! *****************************************************************************
END SUBROUTINE SOLVE_TDMA_ACCORDING_SWEEP_DIRECTION
! *****************************************************************************

! *****************************************************************************
END MODULE GPYRO_SOLID_THERMAL
! *****************************************************************************

