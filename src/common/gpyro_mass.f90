! *****************************************************************************
MODULE GPYRO_MASS
! *****************************************************************************

USE PREC
USE GPYRO_VARS
USE GPYRO_FUNCS, ONLY :  GET_CPU_TIME
USE GPYRO_IO, ONLY: SHUTDOWN_GPYRO

IMPLICIT NONE

CONTAINS

! *****************************************************************************
SUBROUTINE UPDATE_SOLID_SPECIES_EVOLUTION(IMESH)
! This subroutine updates the solid-phase species concentrations by evaluating 
! reaction rates, computing the solid source terms, solving the transport equations
! for species, and updating the mesh.
! *****************************************************************************
INTEGER,  INTENT(IN)  :: IMESH
REAL(EB) :: TMID1, TMID2, TMID3
TYPE(SOLVER_CONVERGENCE_INFO), POINTER :: CI
REAL(EB), POINTER, DIMENSION(:,:,:) :: TS,TG

G=>GPM(IMESH)
CI =>G%CONV_INFO

!The temperature used is that at the last time step, as the new one is not calculated
TS => G%TP
TG => G%TG


IF ((.NOT. CI%CONVERGED_YIS(0)) .AND. (SPROP%NRXN .GT. 0)) THEN

   CALL GET_CPU_TIME(TMID3)

   CALL REACTION_RATE_T(TS)
   CALL REACTION_RATE_Y(G%DTIME_YIS)
   CALL GET_CPU_TIME(TMID1) ; GPG%TUSED(11) = GPG%TUSED(11) + TMID1 - TMID3

   CALL SPECIES_SOURCE_TERMS(TS,TG)
   CALL GET_CPU_TIME(TMID2) ; GPG%TUSED(12) = GPG%TUSED(12) + TMID2 - TMID1

   CALL SOLID_SPECIES_SOLVER(TS, G%DTIME_YIS)
   CALL GET_CPU_TIME(TMID1) ; GPG%TUSED(13) = GPG%TUSED(13) + TMID1 - TMID2

   IF (GPG%SOLVE_POROSITY) CALL SOLVE_POROSITY
   CALL UPDATE_MESH(IMESH)
   IF (G%DEFORMATION .OR. G%NTIMESTEPS .EQ. 0) THEN
      CALL CALCULATE_INTERFACE_MESH_FACTORS
   ENDIF
   CALL GET_CPU_TIME(TMID2) ; GPG%TUSED(14) = GPG%TUSED(14) + TMID2 - TMID1

   GPG%TUSED(10) = GPG%TUSED(10) + TMID2 - TMID3

ENDIF
 ! *****************************************************************************
 END SUBROUTINE UPDATE_SOLID_SPECIES_EVOLUTION
 ! *****************************************************************************

 
! *****************************************************************************
SUBROUTINE REACTION_RATE_T(TS)
! *****************************************************************************

REAL(EB), POINTER, DIMENSION(:,:,:), INTENT(IN) :: TS

REAL(EB), PARAMETER  :: R0 = 8.314E-3
INTEGER :: IZ,IX,IY,IRXN
REAL(EB) :: MULT
INTEGER :: NCELLX,NCELLY,NCELLZ,NRXN


NCELLX = G%NCELLX
NCELLY = G%NCELLY
NCELLZ = G%NCELLZ
NRXN   = SPROP%NRXN


! Calculate temperature part of forward reaction rate (destruction)
! RRT=Z*exp(-TA/T) * RYIDZSIGMA / DLTZ plus O2 dependency, if any
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ,IRXN,MULT) &
!$OMP SHARED(G,SPROP,TS,RXN,GPG,NCELLX,NCELLY,NCELLZ,NRXN) COLLAPSE(3)
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ
   G%IS_REACTING(:,IZ,IX,IY)=.FALSE.
   DO IRXN = 1, NRXN

      G%RRT(IRXN,IZ,IX,IY) = 0D0 
      IF (G%IMASK   (IZ,IX,IY)) CYCLE
      IF (G%CONSUMED(IZ,IX,IY)) CYCLE
      IF (TS(IZ,IX,IY) .GE. RXN(IRXN)%TCRIT) THEN
         MULT = MIN( TS(IZ,IX,IY)-RXN(IRXN)%TCRIT, 1D0)
         G%RRT(IRXN,IZ,IX,IY) = MULT*RXN(IRXN)%Z*EXP(-RXN(IRXN)%E/(R0*TS(IZ,IX,IY))) &
                               * G%RYIDZSIGMA(RXN(IRXN)%IFROM,IZ,IX,IY)/G%DLTZ(IZ,IX,IY)
      ENDIF
      IF (RXN(IRXN)%ORDERO2 .NE. 0D0) THEN
         IF (RXN(IRXN)%IO2TYPE .EQ. 0) THEN ! [(1+YO2)^nO2-1]
            MULT=(1D0+G%YJG(GPROP%IO2,IZ,IX,IY))**RXN(IRXN)%ORDERO2-1D0
         ELSE !YO2^nO2
            MULT=G%YJGN(GPROP%IO2,IZ,IX,IY)**RXN(IRXN)%ORDERO2
         ENDIF
         G%RRT(IRXN,IZ,IX,IY) = G%RRT(IRXN,IZ,IX,IY) * MULT
      ENDIF

      IF (ABS(G%RRT(IRXN,IZ,IX,IY)).GT. EPSILON_FB) THEN
         G%IS_REACTING(IRXN,IZ,IX,IY)=.TRUE.
         G%IS_REACTING(0   ,IZ,IX,IY)=.TRUE.
      ENDIF
   ENDDO
ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO

! *****************************************************************************      
END SUBROUTINE REACTION_RATE_T
! *****************************************************************************

! *****************************************************************************
SUBROUTINE REACTION_RATE_Y(DTIME)
! *****************************************************************************
REAL(EB), INTENT(IN) :: DTIME
REAL(EB)             :: M,N,OMN
REAL(EB)             :: ALPHA, OMALPHA
INTEGER              :: IRXN,IKM,IZ,IX,IY,ISPEC
INTEGER              :: NCELLX,NCELLY,NCELLZ,NSSPEC,NRXN


NCELLX = G%NCELLX
NCELLY = G%NCELLY
NCELLZ = G%NCELLZ
NSSPEC = SPROP%NSSPEC
NRXN   = SPROP%NRXN

! Update RYIDZSIGMA
IF (.NOT. GPG%CONVENTIONAL_RXN_ORDER) THEN
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ,ISPEC) SHARED(G,DTIME,NCELLX,NCELLY,NCELLZ) COLLAPSE(3)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
   IF (.NOT. G%IS_REACTING(0,IZ,IX,IY)) THEN
      G%RYIDZSIGMAN(:,IZ,IX,IY)=G%RYIDZSIGMA(:,IZ,IX,IY) 
      CYCLE
   ENDIF
   DO ISPEC = 1, NSSPEC
      G%RYIDZSIGMAN(ISPEC,IZ,IX,IY)=G%RYIDZSIGMA(ISPEC,IZ,IX,IY) + G%SOMEGA(1,ISPEC,IZ,IX,IY)*G%DLTZ (IZ,IX,IY)*DTIME
   ENDDO
   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDIF

! Determine unreactedness, i.e. 1 - conversion

!$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(3) PRIVATE(IX,IY,IZ,ISPEC) &
!$OMP SHARED(G,SPROP,NCELLX,NCELLY,NCELLZ,NSSPEC)
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ
DO ISPEC = 1, NSSPEC
   IF(G%CONSUMED(IZ,IX,IY) .OR. G%IMASK(IZ,IX,IY) .OR. (.NOT. G%IS_REACTING(0,IZ,IX,IY))) CYCLE
   IF (G%RYIDZSIGMA  (ISPEC,IZ,IX,IY) .GT. 0D0) THEN
      G%UNREACTEDNESS(ISPEC,IZ,IX,IY) = G%RYIDZP(ISPEC,IZ,IX,IY) / G%RYIDZSIGMA(ISPEC,IZ,IX,IY)
   ELSEIF (G%UNREACTEDNESS(ISPEC,IZ,IX,IY) .LT. 0D0) THEN
      G%UNREACTEDNESS(ISPEC,IZ,IX,IY) = 0D0
   ELSE
      G%UNREACTEDNESS(ISPEC,IZ,IX,IY) = 1D0
   ENDIF
ENDDO
ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO

! Calculate Y part of reaction rate

!$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(4) &
!$OMP PRIVATE(IX, IY, IZ, IRXN, M, N, OMN, IKM, OMALPHA, ALPHA) &
!$OMP SHARED(G, RXN, NCELLX, NCELLY, NCELLZ, NRXN)
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ
DO IRXN = 1, NRXN
  ! In some specific cases do nothing
   IF(G%CONSUMED(IZ,IX,IY) .OR. G%IMASK(IZ,IX,IY) .OR. (.NOT. G%IS_REACTING(IRXN,IZ,IX,IY))) THEN
      G%RRY(IRXN,IZ,IX,IY) = 0D0 
      CYCLE
   ENDIF

   M       = RXN(IRXN)%M       ! m 
   N       = RXN(IRXN)%ORDER   ! n 
   OMN     = 1D0 - N           ! 1 - n
   IKM = RXN(IRXN)%IKINETICMODEL
     
   OMALPHA =       G%UNREACTEDNESS(RXN(IRXN)%IFROM,IZ,IX,IY) ! 1 - alpha
   ALPHA   = 1D0 - G%UNREACTEDNESS(RXN(IRXN)%IFROM,IZ,IX,IY) ! alpha
   IF (OMALPHA .LT. 0D0) THEN ; OMALPHA = 0D0 ;  ALPHA   = 1D0 ; ENDIF
   IF (  ALPHA .LT. 0D0) THEN ; ALPHA   = 0D0 ;  OMALPHA = 1D0 ; ENDIF
   
   SELECT CASE (RXN(IRXN)%IKINETICMODEL)
      CASE(0) !Default reaction order treatment
         IF (N .GT. 0.99999 .AND. N .LT. 1.00001) THEN
            G%RRY(IRXN,IZ,IX,IY)= OMALPHA 
         ELSE
            G%RRY(IRXN,IZ,IX,IY)= OMALPHA ** N 
         ENDIF
      CASE(1) !Nucleation and nucleus growing
         G%RRY(IRXN,IZ,IX,IY)= (1D0/N) * OMALPHA * (-LOG(MIN(OMALPHA,1D0-1D-12)))**OMN
      CASE(2) !Phase boundary - same as CASE(0)
         IF (N .GT. 0.99999 .AND. N .LT. 1.00001) THEN
            G%RRY(IRXN,IZ,IX,IY)= OMALPHA
         ELSE
            G%RRY(IRXN,IZ,IX,IY)= OMALPHA ** N 
         ENDIF
      CASE(3) !Diffusion - plane symmetry
         G%RRY(IRXN,IZ,IX,IY) = 0.5D0 / MAX(ALPHA, 1D-12)
      CASE(4) !Diffusion - cylindrical symmetry
         G%RRY(IRXN,IZ,IX,IY) = -1D0 / LOG(MIN(OMALPHA,1D0-1D-12))
      CASE(5) !Diffusion - spherical symmetry
         G%RRY(IRXN,IZ,IX,IY) = 1.5D0/MAX(OMALPHA**(-1D0/3D0)-1D0,1D-12)
      CASE(6) !Diffusion - Jander
         G%RRY(IRXN,IZ,IX,IY) = 1.5D0 * MIN(OMALPHA,1D0-1D-12)**(-1D0/3D0) / MAX(ALPHA,1D-12)
      CASE(7) !Potential law
         G%RRY(IRXN,IZ,IX,IY) = (1D0/N) * MAX(ALPHA,1D-12) ** OMN 
      CASE(8) !Reaction order
         G%RRY(IRXN,IZ,IX,IY) = (1D0/N) * OMALPHA ** OMN
      CASE(9) !Catalytic
         G%RRY(IRXN,IZ,IX,IY) = (OMALPHA ** N) *  (1D0 + RXN(IRXN)%KCAT * (1D0 - G%UNREACTEDNESS(RXN(IRXN)%ICAT,IZ,IX,IY))) 
      CASE DEFAULT 
         CONTINUE 
   END SELECT

   IF (.NOT. (M .GT. -0.00001 .AND. M .LT. 0.00001)) THEN
      ALPHA = MAX(1D-10,ALPHA)
      G%RRY(IRXN,IZ,IX,IY) = G%RRY(IRXN,IZ,IX,IY) * ALPHA ** M 
   ENDIF

ENDDO
ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO

! *****************************************************************************
END SUBROUTINE REACTION_RATE_Y
! *****************************************************************************

! *****************************************************************************      
SUBROUTINE SPECIES_SOURCE_TERMS(TS,TG)
! *****************************************************************************
REAL(EB), POINTER, DIMENSION(:,:,:),INTENT(IN):: TS,TG
REAL(EB) :: RHOFROM, RHOTO
INTEGER, DIMENSION(1:20), SAVE :: IFROM,ITOS
REAL(EB) :: P,Q,B,Z,TA
REAL(EB) :: CA,CB
INTEGER :: NSSPEC, NGSPEC, NCELLZ,NCELLX,NCELLY
INTEGER :: IGSPEC,IZ,IX,IY,IRXN,ISPEC,IA,IB, NRXN


NCELLX = G%NCELLX
NCELLY = G%NCELLY
NCELLZ = G%NCELLZ
NSSPEC = SPROP%NSSPEC
NGSPEC = GPROP%NGSPEC
NRXN   = SPROP%NRXN


IF (.NOT. G%CONV_INFO%CONVERGED_YIS(0)) THEN

   !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) &
   !$OMP PRIVATE(IX,IY,IZ,IRXN,IFROM,ITOS,RHOFROM,RHOTO,IGSPEC,ISPEC) COLLAPSE(3)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
      G%SOMEGA  (:,:,IZ,IX,IY) = 0D0
      G%GOMEGA  (:,:,IZ,IX,IY) = 0D0
      G%OMEGASFG(IZ,IX,IY)     = 0D0

      IF (G%IMASK(IZ,IX,IY)) CYCLE

      IF (.NOT. G%IS_REACTING(0,IZ,IX,IY)) THEN
         G%OMEGASDAK  (:,IZ,IX,IY) = 0D0
         G%OMEGASFBK  (:,IZ,IX,IY) = 0D0
         G%OMEGASFGK  (:,IZ,IX,IY) = 0D0
         G%OMEGASFJK(:,:,IZ,IX,IY) = 0D0
         CYCLE
      ENDIF

      ! First, determine the quantities SOLIDFRAC (solid fraction) and 
      ! OMSOLIDFRAC (one minus solid fraction). Note that to save time,
      ! this could be pre-calculated and stored as a function of temperature
      ! for each reaction. 

      DO IRXN = 1, NRXN
         IF (.NOT. G%IS_REACTING(IRXN,IZ,IX,IY)) CYCLE
         !IF (ITER .EQ. 1) THEN
         IFROM(IRXN) = RXN(IRXN)%IFROM
         ITOS (IRXN) = RXN(IRXN)%ITOS
         RHOFROM = RHOOFT(IFROM(IRXN),TS(IZ,IX,IY))
         IF (ITOS(IRXN) .GT. 0) RHOTO = RHOOFT(ITOS (IRXN),TS(IZ,IX,IY))
         IF (ITOS(IRXN) .LE. 0) RHOTO = 0. 
         G%OMSOLIDFRAC(IRXN,IZ,IX,IY) = -RXN(IRXN)%CHI * (RHOTO/RHOFROM - 1D0)
         !ENDIF 
         G%OMEGASDAK (IRXN,IZ,IX,IY) = G%RRY(IRXN,IZ,IX,IY) * G%RRT(IRXN,IZ,IX,IY)
         G%OMEGASFBK (IRXN,IZ,IX,IY) = (1D0-G%OMSOLIDFRAC(IRXN,IZ,IX,IY))*G%OMEGASDAK(IRXN,IZ,IX,IY)
         G%OMEGASFGK (IRXN,IZ,IX,IY) = G%OMSOLIDFRAC(IRXN,IZ,IX,IY)*G%OMEGASDAK(IRXN,IZ,IX,IY)
         G%OMEGASFG  (     IZ,IX,IY) = G%OMEGASFG(IZ,IX,IY) + G%OMEGASFGK(IRXN,IZ,IX,IY)

         DO IGSPEC = 1, NGSPEC
            IF (ABS(GPROP%YIELDS(IGSPEC,IRXN)) .LT. EPSILON_FB) CYCLE
            IF (GPROP%YIELDS(IGSPEC,IRXN) .GT. 0D0) THEN !Formation of gaseous species
                  G%OMEGASFJK(IGSPEC,IRXN,IZ,IX,IY) = G%OMEGASFGK(IRXN ,IZ,IX,IY) * GPROP%YIELDS(IGSPEC,IRXN)
                  G%GOMEGA   (1   ,IGSPEC,IZ,IX,IY) = G%GOMEGA(1,IGSPEC,IZ,IX,IY) + G%OMEGASFJK (IGSPEC,IRXN,IZ,IX,IY)
            ELSE !Destruction of gaseous species
                  G%OMEGASDJK(IGSPEC,IRXN,IZ,IX,IY) =-G%OMEGASFGK(IRXN,IZ,IX,IY) * GPROP%YIELDS(IGSPEC,IRXN)
                  G%GOMEGA(2,IGSPEC,IZ,IX,IY) = G%GOMEGA(2,IGSPEC,IZ,IX,IY) + G%OMEGASDJK(IGSPEC,IRXN,IZ,IX,IY)
            ENDIF
         ENDDO !IGSPEC
   
         IF (ITOS(IRXN) .NE. 0) THEN
            G%SOMEGA(1,ITOS(IRXN),IZ,IX,IY) = G%SOMEGA(1,ITOS(IRXN),IZ,IX,IY) + G%OMEGASFBK(IRXN,IZ,IX,IY)
         ENDIF
         G%SOMEGA(2,IFROM(IRXN),IZ,IX,IY) = G%SOMEGA(2,IFROM(IRXN),IZ,IX,IY) + G%OMEGASDAK(IRXN,IZ,IX,IY)

      ENDDO !IRXN
   
      DO IGSPEC = 1, NGSPEC
         G%GOMEGA(3,IGSPEC,IZ,IX,IY) = G%GOMEGA(1,IGSPEC,IZ,IX,IY) - G%GOMEGA(2,IGSPEC,IZ,IX,IY)
         G%GOMEGA(3,0     ,IZ,IX,IY) = G%GOMEGA(3,0     ,IZ,IX,IY) + G%GOMEGA(3,IGSPEC,IZ,IX,IY)
      ENDDO
      DO ISPEC = 1, NSSPEC
         G%SOMEGA(3,ISPEC,IZ,IX,IY) = G%SOMEGA(1,ISPEC,IZ,IX,IY) - G%SOMEGA(2,ISPEC,IZ,IX,IY)
      ENDDO
   ENDDO
   ENDDO
   ENDDO   
   !$OMP END PARALLEL DO
ENDIF 

IF (.NOT. G%CONV_INFO%CONVERGED_YJG(0) .AND. GPROP%NHGRXN .GT. 0) THEN

   !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(PRIVATE) SHARED(G,GPROP,HGRXN,GPG,NCELLY,NCELLX,NCELLZ) COLLAPSE(3)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
      ! Calculate homogeneous gas reaction rates
      ! 1 = FORMATION
      ! 2 = DESTRUCTION
      ! 3 = NET
      IF (G%IMASK(IZ,IX,IY)) CYCLE
      DO IRXN = 1, GPROP%NHGRXN
         G%HGOMEGA(:,IRXN,IZ,IX,IY) = 0D0
         IA = HGRXN(IRXN)%IREACTANT1
         IB = HGRXN(IRXN)%IREACTANT2
         P  = HGRXN(IRXN)%P
         Q  = HGRXN(IRXN)%Q
         B  = HGRXN(IRXN)%B
         Z  = HGRXN(IRXN)%Z
         TA = 1D3*HGRXN(IRXN)%E / 8.314D0

         CA = 1D3 * G%RGN(IZ,IX,IY) * G%YJGN(IA,IZ,IX,IY) / GPROP%M(IA)
         CB = 1D3 * G%RGN(IZ,IX,IY) * G%YJGN(IB,IZ,IX,IY) / GPROP%M(IB)
            
         ! Try to limit exponentiation where possible (cpu time):
         IF (P .EQ. 1D0) G%HGRR(IRXN,IZ,IX,IY) = CA
         IF (P .NE. 1D0) G%HGRR(IRXN,IZ,IX,IY) = CA**P
            
         IF (Q .EQ. 1D0) G%HGRR(IRXN,IZ,IX,IY) = G%HGRR(IRXN,IZ,IX,IY) * CB
         IF (Q .NE. 1D0) G%HGRR(IRXN,IZ,IX,IY) = G%HGRR(IRXN,IZ,IX,IY) * CB**Q
            
         IF (B .NE. 0D0) THEN
            IF (GPG%THERMAL_EQUILIBRIUM) THEN
               G%HGRR(IRXN,IZ,IX,IY) = G%HGRR(IRXN,IZ,IX,IY) * TS(IZ,IX,IY)**B
            ELSE
               G%HGRR(IRXN,IZ,IX,IY) = G%HGRR(IRXN,IZ,IX,IY) * TG(IZ,IX,IY)**B
            ENDIF
         ENDIF
            
         IF (GPG%THERMAL_EQUILIBRIUM) THEN
            G%HGRR(IRXN,IZ,IX,IY) = G%HGRR(IRXN,IZ,IX,IY) * Z * EXP(-TA/TS(IZ,IX,IY))
         ELSE
            G%HGRR(IRXN,IZ,IX,IY) = G%HGRR(IRXN,IZ,IX,IY) * Z * EXP(-TA/TG(IZ,IX,IY))
         ENDIF
            
         G%HGRR(IRXN,IZ,IX,IY) = G%HGRR(IRXN,IZ,IX,IY) * G%POROSS(IZ,IX,IY)
                          
         DO IGSPEC = 1, GPROP%NGSPEC
            IF (GPROP%HGYIELDS(IGSPEC,IRXN) .EQ. 0D0) CYCLE
            IF (GPROP%HGYIELDS(IGSPEC,IRXN) .GT. 0D0) THEN !Formation
               G%HGOMEGA(1,IGSPEC,IZ,IX,IY) = G%HGOMEGA(1,IGSPEC,IZ,IX,IY) &
                                            + G%HGRR(IRXN,IZ,IX,IY)*GPROP%HGYIELDS(IGSPEC,IRXN)
            ELSE !Destruction
               G%OMEGAGDJL(IGSPEC,IRXN,IZ,IX,IY) = -G%HGRR(IRXN,IZ,IX,IY) * GPROP%HGYIELDS(IGSPEC,IRXN)
               G%HGOMEGA  (2,IGSPEC,IZ,IX,IY)    =  G%HGOMEGA(2,IGSPEC,IZ,IX,IY) + G%OMEGAGDJL(IGSPEC,IRXN,IZ,IX,IY)
            ENDIF
         ENDDO !IGSPEC 
 
         ! Net rate:      
         G%HGOMEGA(3,IGSPEC,IZ,IX,IY) = G%HGOMEGA(1,IGSPEC,IZ,IX,IY) - G%HGOMEGA(2,IGSPEC,IZ,IX,IY)
   
      ENDDO
   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDIF 

! *****************************************************************************
END SUBROUTINE SPECIES_SOURCE_TERMS
! *****************************************************************************


! *****************************************************************************
SUBROUTINE SOLID_SPECIES_SOLVER(TS,DTIME)
! *****************************************************************************
REAL(EB), INTENT(IN) :: DTIME
REAL(EB), POINTER, DIMENSION(:,:,:),INTENT(IN):: TS

INTEGER     :: ISPEC,IZ,IX,IY,NB_ITER_TO_ADJUST
INTEGER     :: NSSPEC,NCELLZ,NCELLX,NCELLY
REAL(EB)    :: T
LOGICAL     :: LOCALY_CONVERGED

NCELLZ = G%NCELLZ
NCELLX = G%NCELLX
NCELLY = G%NCELLY
NSSPEC = SPROP%NSSPEC


G%CONV_INFO%CONVERGED_YIS(:) = .TRUE.


! Continuity equation
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ,ISPEC,T,LOCALY_CONVERGED,NB_ITER_TO_ADJUST)&
!$OMP DEFAULT(SHARED) COLLAPSE(3)
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ
    IF (G%IMASK(IZ,IX,IY)) CYCLE
    T=TS(IZ,IX,IY)
    CALL SOLVE_SPECIES_CONSERVATION_IN_CELL(IZ,IX,IY,T,DTIME,LOCALY_CONVERGED)

    !Adjust if there is a local issue with the conservation of species
    IF (.NOT. LOCALY_CONVERGED) THEN
        NB_ITER_TO_ADJUST = 0
        DO ISPEC = 1, SPROP%NSSPEC
            CALL ADJUST_SPECIES_IN_CELL(ISPEC, IZ, IX, IY, DTIME, NB_ITER_TO_ADJUST)
        ENDDO
        IF (NB_ITER_TO_ADJUST .GT. 2 * NSSPEC) THEN
        G%CONV_INFO%CONVERGED_YIS(:) = .FALSE.
        CYCLE
        ENDIF
        CALL SOLVE_SPECIES_CONSERVATION_IN_CELL(IZ,IX,IY,T,DTIME,LOCALY_CONVERGED)
        IF (.NOT. LOCALY_CONVERGED) G%CONV_INFO%CONVERGED_YIS(:) = .FALSE.
    ENDIF

ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO


! None iterative solver 
DO ISPEC = 0, NSSPEC
   IF (G%CONV_INFO%CONVERGED_YIS(ISPEC) .AND. G%CONV_INFO%ITER_YIS(ISPEC) .EQ. 0) THEN
      G%CONV_INFO%ITER_YIS(ISPEC) = G%CONV_INFO%ITER
   ENDIF
ENDDO


! *****************************************************************************
END SUBROUTINE SOLID_SPECIES_SOLVER
! *****************************************************************************
 
! *****************************************************************************
SUBROUTINE SOLVE_SPECIES_CONSERVATION_IN_CELL(IZ,IX,IY,T,DTIME,CONVERGED)
! *****************************************************************************
LOGICAL, INTENT(OUT) :: CONVERGED
REAL(EB), INTENT(IN) :: DTIME
INTEGER, INTENT(IN)  :: IZ,IX,IY
REAL(EB), INTENT(IN) :: T
REAL(EB)    :: SUM_YI_O_RHOI, RHOI,YI_SUM
INTEGER     :: ISPEC,NSSPEC
LOGICAL     :: NO_MORE_MASS
REAL(EB)    :: MASSE
REAL(EB), DIMENSION(SPROP%NSSPEC) :: MASSE_I
CHARACTER(LEN=512) :: MESSAGE

NSSPEC = SPROP%NSSPEC


! By default converged = TRUE
CONVERGED = .TRUE.

!If no reaction occurs in the cell, Yi, rho, etc., keep the same value from the previous iteration.
IF (.NOT. G%IS_REACTING(0,IZ,IX,IY)) THEN
   G%RDLTZN(IZ,IX,IY) = G%RDLTZ(IZ,IX,IY)
   G%RPN   (IZ,IX,IY) = G%RP   (IZ,IX,IY)
   DO ISPEC = 1, NSSPEC
      G%RYIDZPN(ISPEC,IZ,IX,IY) = G%RYIDZP(ISPEC,IZ,IX,IY)
      G%YIN    (ISPEC,IZ,IX,IY) = G%YI    (ISPEC,IZ,IX,IY)
      G%XIN    (ISPEC,IZ,IX,IY) = G%XI    (ISPEC,IZ,IX,IY)
   ENDDO
   RETURN
ENDIF


! Calculate the Rho*Dz (The total masse in the cell)
G%RDLTZN(IZ,IX,IY) = -(DTIME * G%GOMEGA(3,0,IZ,IX,IY) * G%DLTZ(IZ,IX,IY)) + G%RDLTZ(IZ,IX,IY)
MASSE=G%RDLTZN(IZ,IX,IY)
! check if mass remains in the cell
NO_MORE_MASS = .FALSE.
IF ( MASSE .LT. EPSILON_EB) NO_MORE_MASS = .TRUE.

SUM_YI_O_RHOI  = 0D0
YI_SUM         = 0D0 ! Help to check Final convergence (coverged if YI_SUM=1)

DO ISPEC = 1, NSSPEC
   ! CLACLUATE Rho*Yi*Dz (The masse of species i)
   G%RYIDZPN(ISPEC,IZ,IX,IY) = DTIME *G%SOMEGA(3,ISPEC,IZ,IX,IY) * G%DLTZ(IZ,IX,IY) &
                              + G%RYIDZP(ISPEC,IZ,IX,IY)
   MASSE_I(ISPEC) = G%RYIDZPN(ISPEC,IZ,IX,IY)

   ! Check that no surestimated mass was destroyed.
   IF (G%RYIDZPN(ISPEC,IZ,IX,IY) .LT. - EPSILON_EB) THEN
      ! If it occurs, we will call adjust ADJUST_SPECIES_IN_CELL
      CONVERGED = .FALSE.
      RETURN
   ELSEIF (G%RYIDZPN(ISPEC,IZ,IX,IY) .LT. 0D0 .AND.  G%RYIDZPN(ISPEC,IZ,IX,IY) .GT. - EPSILON_EB )THEN
      ! This is if the numerical approximation give a very small numerical negatif value
      G%RYIDZPN(ISPEC,IZ,IX,IY) = 0D0 
   ENDIF

   IF (NO_MORE_MASS) CYCLE

   ! Determine the masse fraction of species i (Y_i)

    G%YIN(ISPEC,IZ,IX,IY) = G%RYIDZPN(ISPEC,IZ,IX,IY) / G%RDLTZN(IZ,IX,IY)

   IF (G%YIN(ISPEC,IZ,IX,IY).NE.G%YIN(ISPEC,IZ,IX,IY) .OR. &
      G%YIN(ISPEC,IZ,IX,IY).EQ.GPG%POSINF .OR. G%YIN(ISPEC,IZ,IX,IY).EQ.GPG%NEGINF) THEN
      GPG%NAN = .TRUE.
   ENDIF

   YI_SUM = YI_SUM + G%YIN(ISPEC,IZ,IX,IY)

   ! To determine the bulk density
   IF (G%YIN(ISPEC,IZ,IX,IY) .EQ. 0) CYCLE ! No more species Yi in the cell
   RHOI  = RHOOFT(ISPEC,T)
   SUM_YI_O_RHOI  = SUM_YI_O_RHOI + G%YIN(ISPEC,IZ,IX,IY)/RHOI


ENDDO

IF (NO_MORE_MASS) THEN
   
   ! If the cell is fully depleted and no mass remains, it creates a critical issue for the solver,
   ! as handling zero-volume cells is currently not supported.

   ! As a temporary workaround, when the mass in a cell falls below EPSILON_EB,
   ! chemical reactions are halted in that cell and a minimal residual mass is retained
   ! to prevent the volume from becoming zero.

   ! If this condition is triggered, it means that EPSILON_EB was set too low.
   ! The current strategy has proven to yield accurate results if EPSILON_EB is sufficiently small,
   ! such that the residual mass is negligible in terms of physical impact,
   ! but still high enough to avoid the creation of empty cells.

   !$OMP CRITICAL
   IF (.NOT. GPG%SHUTDOWN_ALREADY_TRIGGERED) THEN
      GPG%SHUTDOWN_ALREADY_TRIGGERED = .TRUE.
   
      WRITE(MESSAGE,'(A)') NEW_LINE('A') // &
                           " WARNING: Mass in cell has dropped to zero, causing convergence issues."        // NEW_LINE('A') // &
                           " Please increase the value of  &GPYRO_GENERAL%EPS to avoid zero-mass cells."    // NEW_LINE('A') // &
                           " Note: &GPYRO_GENERAL%EPS is the minimum residual mass in a cell divided by its initial mass " // &
                           "that prevents divergence problems." // &
                           NEW_LINE('A')
      
      WRITE(MESSAGE(LEN_TRIM(MESSAGE)+1:), '(A,E10.1)') " Current value of EPS is ", GPG%EPS
      IF (IGPYRO_TYPE .EQ. 1) CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF
   !$OMP END CRITICAL
ENDIF

! Determine the bulk density of the cell
G%RPN(IZ,IX,IY) = 1D0 / SUM_YI_O_RHOI

! Calculate volume fraction:
DO ISPEC = 1, NSSPEC
   G%XIN(ISPEC,IZ,IX,IY)=G%RPN(IZ,IX,IY)*G%YIN(ISPEC,IZ,IX,IY) / RHOOFT(ISPEC,T) 
ENDDO

! Convergence tolerance = 0.01%
IF (YI_SUM .LT. 0.9999 .OR. YI_SUM .GT. 1.0001) THEN
   CONVERGED = .FALSE.
ENDIF

! *****************************************************************************
END SUBROUTINE SOLVE_SPECIES_CONSERVATION_IN_CELL
! *****************************************************************************


! *****************************************************************************
RECURSIVE SUBROUTINE ADJUST_SPECIES_IN_CELL(SPECIES,IZ,IX,IY,TIME_STEP,ITERATION)
! *****************************************************************************
INTEGER, INTENT(IN)    :: SPECIES, IZ, IX, IY
REAL(EB), INTENT(IN)   :: TIME_STEP
INTEGER, INTENT(INOUT) :: ITERATION

REAL(EB) :: MASS_DESTROYED, MASS_AVAILABLE, MASS_REMAINING
REAL(EB) :: EXCESS_MASS_DESTROYED
REAL(EB) :: OVERESTIMATED_ISPEC_DESTRUCTION_RATE
REAL(EB) :: ADJUSTMENT_FACTOR
REAL(EB) :: OLD_OMEGASFBK, OVERESTIMATED_SOLID_B_FORMATION
REAL(EB) :: OLD_OMEGASFGK, OVERESTIMATED_TOTAL_GAS_FORMATION

REAL(EB) :: OLD_OMEGASFJK, OVERESTIMATED_GAS_ISPEC_FORMATION
REAL(EB) :: OLD_OMEGASDJK, OVERESTIMATED_GAS_ISPEC_DESTRUCTION 

INTEGER  :: NRXN, NGSPEC
INTEGER  :: IRXN, IFROM, ITOS, IGSPEC

MASS_DESTROYED = -TIME_STEP * G%SOMEGA(3, SPECIES, IZ, IX, IY) * G%DLTZ(IZ, IX, IY)
MASS_AVAILABLE = G%RYIDZP(SPECIES, IZ, IX, IY)
MASS_REMAINING = MASS_AVAILABLE - MASS_DESTROYED

IF (MASS_REMAINING .GE. 0D0) RETURN ! Mass balance is already satisfied

NRXN   = SPROP%NRXN
NGSPEC = GPROP%NGSPEC

ITERATION = ITERATION + 1

! Stop if too many iterations are needed
IF (ITERATION .GT. 2 * SPROP%NSSPEC) THEN
   RETURN
ENDIF

! Excess mass destroyed beyond availability
EXCESS_MASS_DESTROYED = -MASS_REMAINING

! Adjust net production rate to ensure mass balance
G%SOMEGA(3, SPECIES, IZ, IX, IY) = -MASS_AVAILABLE / (TIME_STEP * G%DLTZ(IZ, IX, IY))

! Adjust destruction rate to avoid exceeding available mass
OVERESTIMATED_ISPEC_DESTRUCTION_RATE = G%SOMEGA(2, SPECIES, IZ, IX, IY)
G%SOMEGA(2, SPECIES, IZ, IX, IY) = OVERESTIMATED_ISPEC_DESTRUCTION_RATE &
                               - EXCESS_MASS_DESTROYED / (TIME_STEP * G%DLTZ(IZ, IX, IY))
ADJUSTMENT_FACTOR = G%SOMEGA(2, SPECIES, IZ, IX, IY) / OVERESTIMATED_ISPEC_DESTRUCTION_RATE

! Adjust reaction rates to maintain conservation
DO IRXN = 1, NRXN
   IFROM = RXN(IRXN)%IFROM
   IF (IFROM .NE. SPECIES) CYCLE 
   IF (.NOT. G%IS_REACTING(IRXN,IZ,IX,IY)) CYCLE

   ! Adjust reaction rates for solid species transformation
   ITOS = RXN(IRXN)%ITOS
   G%OMEGASDAK(IRXN, IZ, IX, IY) = ADJUSTMENT_FACTOR * G%OMEGASDAK(IRXN, IZ, IX, IY)

   OLD_OMEGASFBK = G%OMEGASFBK(IRXN, IZ, IX, IY)
   G%OMEGASFBK(IRXN, IZ, IX, IY) = ADJUSTMENT_FACTOR * G%OMEGASFBK(IRXN, IZ, IX, IY)
   OVERESTIMATED_SOLID_B_FORMATION = OLD_OMEGASFBK - G%OMEGASFBK(IRXN, IZ, IX, IY)



   IF (ITOS .NE. 0) THEN
      G%SOMEGA(1, ITOS, IZ, IX, IY) = G%SOMEGA(1, ITOS, IZ, IX, IY) - OVERESTIMATED_SOLID_B_FORMATION
      G%SOMEGA(3, ITOS, IZ, IX, IY) = G%SOMEGA(3, ITOS, IZ, IX, IY) - OVERESTIMATED_SOLID_B_FORMATION

      CALL ADJUST_SPECIES_IN_CELL(ITOS, IZ, IX, IY, TIME_STEP, ITERATION)
      IF (ITERATION .GT. 2 * SPROP%NSSPEC) RETURN
   ENDIF

   ! Adjust gas formation rates
   OLD_OMEGASFGK = G%OMEGASFGK(IRXN, IZ, IX, IY)
   G%OMEGASFGK(IRXN, IZ, IX, IY) = ADJUSTMENT_FACTOR * OLD_OMEGASFGK 
   OVERESTIMATED_TOTAL_GAS_FORMATION = OLD_OMEGASFGK - G%OMEGASFGK(IRXN, IZ, IX, IY)

   G%OMEGASFG(IZ, IX, IY) = G%OMEGASFG(IZ, IX, IY) - OVERESTIMATED_TOTAL_GAS_FORMATION
   G%GOMEGA(3, 0, IZ, IX, IY) = G%GOMEGA(3, 0, IZ, IX, IY) - OVERESTIMATED_TOTAL_GAS_FORMATION

   ! Adjust individual gas species rates
   DO IGSPEC = 1, NGSPEC
      IF (ABS(GPROP%YIELDS(IGSPEC,IRXN)) .LT. EPSILON_FB) CYCLE
      IF (GPROP%YIELDS(IGSPEC,IRXN) .GT. 0D0) THEN !Formation of gaseous species
            OLD_OMEGASFJK = G%OMEGASFJK(IGSPEC,IRXN,IZ,IX,IY)
            G%OMEGASFJK(IGSPEC,IRXN,IZ,IX,IY) = ADJUSTMENT_FACTOR * G%OMEGASFJK(IGSPEC,IRXN,IZ,IX,IY)
            OVERESTIMATED_GAS_ISPEC_FORMATION = OLD_OMEGASFJK - G%OMEGASFJK(IGSPEC,IRXN,IZ,IX,IY)
            G%GOMEGA   (1   ,IGSPEC,IZ,IX,IY) = G%GOMEGA(1,IGSPEC,IZ,IX,IY) - OVERESTIMATED_GAS_ISPEC_FORMATION
            G%GOMEGA   (3   ,IGSPEC,IZ,IX,IY) = G%GOMEGA(3,IGSPEC,IZ,IX,IY) - OVERESTIMATED_GAS_ISPEC_FORMATION
            
      ELSE !Destruction of gaseous species
            OLD_OMEGASDJK = G%OMEGASDJK(IGSPEC,IRXN,IZ,IX,IY)
            G%OMEGASDJK(IGSPEC,IRXN,IZ,IX,IY) = ADJUSTMENT_FACTOR * G%OMEGASDJK(IGSPEC,IRXN,IZ,IX,IY)
            OVERESTIMATED_GAS_ISPEC_DESTRUCTION = OLD_OMEGASDJK - G%OMEGASDJK(IGSPEC,IRXN,IZ,IX,IY)
            G%GOMEGA(2,IGSPEC,IZ,IX,IY) = G%GOMEGA(2,IGSPEC,IZ,IX,IY) - OVERESTIMATED_GAS_ISPEC_DESTRUCTION
            G%GOMEGA(3,IGSPEC,IZ,IX,IY) = G%GOMEGA(3,IGSPEC,IZ,IX,IY) + OVERESTIMATED_GAS_ISPEC_DESTRUCTION
            
      ENDIF
   ENDDO

ENDDO

! *****************************************************************************
END SUBROUTINE ADJUST_SPECIES_IN_CELL
! *****************************************************************************


! *****************************************************************************
SUBROUTINE UPDATE_MESH(IMESH)
! *****************************************************************************
INTEGER, INTENT(IN)  :: IMESH
INTEGER              :: IX,IY,IZ
INTEGER              :: NCELLX,NCELLY,NCELLZ
INTEGER              :: IOR, ICOUNT


NCELLX= G%NCELLX
NCELLY= G%NCELLY
NCELLZ= G%NCELLZ

IF (G%DIMENSION .EQ. 0) THEN  !OD
   G%DLTZN (1,1,1) = G%RDLTZN(1,1,1) / G%RPN(1,1,1)
ENDIF


IF (G%DIMENSION .EQ. 1) THEN  !1D
   IX=1 ; IY=1

   ! Get new Dz from rho and rho*Dz
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IZ) SHARED(G,NCELLZ,IX,IY)
   DO IZ = 1, NCELLZ
      ! IF (G%RDLTZN(IZ,IX,IY) .LT. EPSILON_EB) THEN
      !    G%DLTZN (IZ,IX,IY) = 0D0
      !    G%IMASK(IZ,IX,IY)  = .TRUE.
      !    IF (G%NEEDSBCB(IZ,IX,IY)) THEN
      !       G%NEEDSBCB(IZ,IX,IY) = .FALSE.
      !       G%NEEDSBCB(IZ-1,IX,IY) = .TRUE.
      !    ENDIF
      !    IF (G%NEEDSBCT(IZ,IX,IY)) THEN
      !       G%NEEDSBCT(IZ,IX,IY) = .FALSE.
      !       G%NEEDSBCT(IZ+1,IX,IY) = .TRUE.
      !    ENDIF
      !    CYCLE
      ! ENDIF
      G%DLTZN (IZ,IX,IY) = G%RDLTZN(IZ,IX,IY) / G%RPN(IZ,IX,IY)
   ENDDO
   !$OMP END PARALLEL DO


   DO ICOUNT = 1, NGPYRO_FACES_NEEDING_BCS(IMESH)

      IF(.NOT. GP_BOUDARYS(IMESH)%COMPLETE_CELL_AT_BC(ICOUNT)) CYCLE
   
      IZ  = GP_BOUDARYS(IMESH)%IZ_GPYRO (ICOUNT)
      IX  = GP_BOUDARYS(IMESH)%IX_GPYRO (ICOUNT)
      IY  = GP_BOUDARYS(IMESH)%IY_GPYRO (ICOUNT)
      IOR = GP_BOUDARYS(IMESH)%IOR_GPYRO(ICOUNT)
   
      SELECT CASE(IOR)
         CASE( 3) ; G%DLTZN(IZ-1,IX,IY) = G%DLTZN(IZ,IX,IY)
         CASE(-3) ; G%DLTZN(IZ+1,IX,IY) = G%DLTZN(IZ,IX,IY)
         CASE( 1) ; G%DLTZN(IZ,IX+1,IY) = G%DLTZN(IZ,IX,IY)
         CASE(-1) ; G%DLTZN(IZ,IX-1,IY) = G%DLTZN(IZ,IX,IY)
         CASE( 2) ; G%DLTZN(IZ,IX,IY+1) = G%DLTZN(IZ,IX,IY)
         CASE(-2) ; G%DLTZN(IZ,IX,IY-1) = G%DLTZN(IZ,IX,IY)
      END SELECT
   
   END DO
   

   ! Calculate new z values from new Dz
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IZ) SHARED(G,NCELLZ,IX,IY)
   DO IZ = 3, NCELLZ-1
      G%DZT(IZ,IX,IY) = 0.5D0 * (G%DLTZN(IZ,IX,IY) + G%DLTZN(IZ-1,IX,IY))
   ENDDO
   !$OMP END PARALLEL DO

   G%DZT(2     ,IX,IY) = G%DLTZN(1     ,IX,IY) + 0.5D0 * G%DLTZN(2       ,IX,IY)
   G%DZT(NCELLZ,IX,IY) = G%DLTZN(NCELLZ,IX,IY) + 0.5D0 * G%DLTZN(NCELLZ-1,IX,IY)
   G%DZT(1     ,IX,IY) = G%DZT  (2     ,IX,IY)

   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IZ) SHARED(G,NCELLZ,IX,IY)
   DO IZ = 1, NCELLZ - 1
      G%DZB(IZ,IX,IY) = G%DZT(IZ+1,IX,IY)
   ENDDO
   !$OMP END PARALLEL DO

   G%DZB(NCELLZ,IX,IY) = G%DZB(NCELLZ-1,IX,IY)

   G%DV => G%DLTZN

ENDIF

IF (G%DIMENSION .GT. 1 ) THEN  ! 2D/3D
   IF (.NOT. G%DEFORMATION) THEN
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) &
      !$OMP  COLLAPSE(3)
      DO IY = 1, NCELLY
      DO IX = 1, NCELLX
      DO IZ = 1, NCELLZ
         G%DLTZN (IZ,IX,IY) = G%DLTZ(IZ,IX,IY) 
         G%RPN(IZ,IX,IY)    =  G%RDLTZN(IZ,IX,IY) / G%DLTZN (IZ,IX,IY)
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   ENDIF
ENDIF
!! WARNING : no update mesh in 2D and 3D cases, to fix

! *****************************************************************************
END SUBROUTINE UPDATE_MESH
! *****************************************************************************

! *****************************************************************************
SUBROUTINE CALCULATE_INTERFACE_MESH_FACTORS()
! *****************************************************************************
! Calculate weighting factors

INTEGER :: NCELLZ,NCELLX,NCELLY,IZ,IX,IY


NCELLZ = G%NCELLZ
NCELLX = G%NCELLX
NCELLY = G%NCELLY

IF (NCELLZ .GT. 1) THEN
    !z-direction:
    G%FB(NCELLZ,:,:) = 0D0
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IZ) SHARED(G,NCELLZ)
    DO IZ = 1, NCELLZ-1
    G%FB(IZ,:,:) = G%DLTZN(IZ+1,:,:) / (G%DLTZN(IZ+1,:,:) + G%DLTZN(IZ,:,:))
    ENDDO
    !$OMP END PARALLEL DO 

    G%FT(1,:,:) = 0D0
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IZ) SHARED(G,NCELLZ)
    DO IZ = 2, NCELLZ
    G%FT(IZ,:,:) = G%DLTZN(IZ-1,:,:) / (G%DLTZN(IZ-1,:,:) + G%DLTZN(IZ,:,:))
    ENDDO
    !$OMP END PARALLEL DO 
ENDIF

!x-direction:
! As in 2D/3D ther is no deformation in x only calculate at the first time step
IF ((G%NTIMESTEPS .EQ. 0) .AND. (NCELLX .GT. 1)) THEN
   G%FE(:,NCELLX,:) = 0D0
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX) SHARED(G,NCELLX)
   DO IX = 1, NCELLX-1
      G%FE(:,IX,:) = G%DLTX (:,IX+1,:) / (G%DLTX (:,IX+1,:) + G%DLTX (:,IX,:))
   ENDDO
   !$OMP END PARALLEL DO 


   G%FW(:,1,:) = 0D0
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX) SHARED(G,NCELLX)
   DO IX = 2, NCELLX
      G%FW(:,IX,:) = G%DLTX (:,IX-1,:) / (G%DLTX (:,IX-1,:) + G%DLTX (:,IX,:))
   ENDDO
   !$OMP END PARALLEL DO 

ENDIF

!y-direction:
! As in 2D/3D ther is no deformation only calculate at the first time step
IF ((G%NTIMESTEPS .EQ. 0) .AND. (NCELLY .GT. 1 )) THEN
   G%FN(:,:,NCELLY) = 0D0
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IY) SHARED(G,NCELLY)
   DO IY = 1, NCELLY-1
      G%FN(:,:,IY) = G%DLTY (:,:,IY+1) / (G%DLTY (:,:,IY+1) + G%DLTY (:,:,IY))
   ENDDO
   !$OMP END PARALLEL DO 

   G%FS(:,:,1) = 0D0
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX) SHARED(G,NCELLX)
   DO IY = 2, NCELLY
      G%FS(:,:,IY) = G%DLTY (:,:,IY-1) / (G%DLTY (:,:,IY-1) + G%DLTY (:,:,IY))
   ENDDO
   !$OMP END PARALLEL DO 
ENDIF

! *****************************************************************************
END SUBROUTINE CALCULATE_INTERFACE_MESH_FACTORS
! *****************************************************************************

! *****************************************************************************
SUBROUTINE SOLVE_POROSITY
! *****************************************************************************
!Get the pure solid density and the porosity

REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YI
INTEGER :: NCELLZ,NCELLX,NCELLY,NSSPEC
INTEGER :: IZ,IX,IY,ISPEC
REAL(EB):: SUM_YI_O_RHOSI


NCELLZ = G%NCELLZ
NCELLX = G%NCELLX
NCELLY = G%NCELLY
NSSPEC = SPROP%NSSPEC

!The Solide masse fraction used is the new one that is alredy calculated
YI => G%YIN

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(SUM_YI_O_RHOSI,IX,IY,IZ,ISPEC) &
!$OMP SHARED(YI,G,NSSPEC,NCELLX,NCELLY,NCELLZ) COLLAPSE(3)
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ
   IF (G%IMASK(IZ,IX,IY)) CYCLE
   SUM_YI_O_RHOSI = 0
   DO ISPEC = 1, NSSPEC
      IF (YI(ISPEC,IZ,IX,IY) .EQ. 0) CYCLE ! No more species Yi in the cell
      SUM_YI_O_RHOSI =SUM_YI_O_RHOSI + YI(ISPEC,IZ,IX,IY) / SPROP%RS0(ISPEC)
   ENDDO
   G%RSPN(IZ,IX,IY) = 1D0 / SUM_YI_O_RHOSI
   G%POROSSN(IZ,IX,IY) = 1D0 - G%RPN(IZ,IX,IY) / G%RSPN(IZ,IX,IY)
ENDDO
ENDDO
ENDDO

! *****************************************************************************
END SUBROUTINE SOLVE_POROSITY
! *****************************************************************************


! *****************************************************************************
REAL(EB) FUNCTION RHOOFT(ISPEC,TMP) 
! *****************************************************************************

INTEGER,  INTENT(IN) :: ISPEC 
REAL(EB), INTENT(IN) :: TMP

IF (SPROP%NR(ISPEC) .EQ. 0D0) THEN
   RHOOFT = SPROP%R0(ISPEC)
ELSE
   RHOOFT = SPROP%R0(ISPEC) * (TMP/GPG%TREF)**SPROP%NR(ISPEC)
ENDIF

! *****************************************************************************
END FUNCTION RHOOFT
! *****************************************************************************

! *****************************************************************************
END MODULE GPYRO_MASS
! *****************************************************************************
