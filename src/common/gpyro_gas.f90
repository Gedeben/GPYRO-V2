! *****************************************************************************
MODULE GPYRO_GAS
! *****************************************************************************

USE PREC
USE GPYRO_VARS
USE GPYRO_FUNCS
USE GPYRO_BC, ONLY : SET_BC_FIXED_VALUE
IMPLICIT NONE

CONTAINS


! *****************************************************************************
SUBROUTINE PRESSURE_SOLVER(IMESH,DTIME,ALPHA_P)
! *****************************************************************************
INTEGER, INTENT(IN) :: IMESH
REAL(EB), INTENT(IN) :: DTIME,ALPHA_P
REAL(EB), POINTER, DIMENSION(:,:,:) :: AP0N,AP0,AP,AB,AT,AE,AW,AN,AS,B,PNOLD,PTR,BX,BY !,BZ
REAL(EB), POINTER, DIMENSION(:,:,:) :: TG, TGN

REAL(EB) :: DPMAX, DS_OR
INTEGER :: NCELLZ,NCELLX,NCELLY
INTEGER :: IZ,IX,IY, IOR, ICOUNT

REAL(EB), PARAMETER :: EPS_P = 0.01D0
LOGICAL :: ALL_CELLS_CONVERGED


IF (GPG%THERMAL_EQUILIBRIUM) THEN
   TG  => G%TPN
   TGN => G%TPN
ELSE
   TG  => G%TG
   TGN => G%TGN
ENDIF


AP0N   => G%RWORK01
AP0    => G%RWORK02
AP     => G%RWORK03
AB     => G%RWORK04
AT     => G%RWORK05
AE     => G%RWORK06
AW     => G%RWORK07
AN     => G%RWORK08
AS     => G%RWORK09
B      => G%RWORK10

PNOLD  => G%RWORK13
BX     => G%RWORK14
BY     => G%RWORK16
!BZ     => G%RWORK18


G=>GPM(IMESH)
GPBCP(1:,1:,1:,-3:)=>G%GPYRO_BOUNDARY_CONDITION(1:,1:,1:,-3:)
NCELLZ = G%NCELLZ
NCELLX = G%NCELLX
NCELLY = G%NCELLY


CALL GET_MEAN_MOLECULARE_WEIGHT()
CALL CALCULATE_WEIGHTED_POINT_PERMEABILITY(IMESH)
CALL CALCULATE_GAS_INTERFACE_QUANTITIES(IMESH,'PERMONU  ')
IF ( G%GZ .NE. 0D0) CALL CALCULATE_GAS_INTERFACE_QUANTITIES(IMESH,'RG       ') 
CALL CALCULATE_GAS_INTERFACE_QUANTITIES(IMESH,'PSIRGD   ')




!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) SHARED(AP0,AP0N,DTIME,G,GPROP,NCELLX,NCELLY,NCELLZ) COLLAPSE(3)
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ
   AP0 (IZ,IX,IY) = (G%M (IZ,IX,IY)/R) * G%POROSS (IZ,IX,IY) * G%DLTZ (IZ,IX,IY) * G%DLTX (IZ,IX,IY) * G%DLTY (IZ,IX,IY) / (TG (IZ,IX,IY) * DTIME)
   AP0N(IZ,IX,IY) = (G%MN(IZ,IX,IY)/R) * G%POROSSN(IZ,IX,IY) * G%DLTZN(IZ,IX,IY) * G%DLTX (IZ,IX,IY) * G%DLTY (IZ,IX,IY) / (TGN(IZ,IX,IY) * DTIME)
ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO


! Specify coefficients for TDMA

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) DEFAULT(SHARED) COLLAPSE(3)
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ
   IF (G%IMASK(IZ,IX,IY)) THEN
      AB(IZ,IX,IY) = 0D0; AT(IZ,IX,IY) = 0D0
      AE(IZ,IX,IY) = 0D0; AW(IZ,IX,IY) = 0D0
      AN(IZ,IX,IY) = 0D0; AS(IZ,IX,IY) = 0D0
      CYCLE
   ENDIF
   AB(IZ,IX,IY) = G%PERMONUB(IZ,IX,IY) * G%DXDY(IZ,IX,IY) / G%DZB(IZ,IX,IY)
   AT(IZ,IX,IY) = G%PERMONUT(IZ,IX,IY) * G%DXDY(IZ,IX,IY) / G%DZT(IZ,IX,IY)
   AE(IZ,IX,IY) = G%PERMONUE(IZ,IX,IY) * G%DYDZ(IZ,IX,IY) / G%DXE(IZ,IX,IY)
   AW(IZ,IX,IY) = G%PERMONUW(IZ,IX,IY) * G%DYDZ(IZ,IX,IY) / G%DXW(IZ,IX,IY)
   AN(IZ,IX,IY) = G%PERMONUN(IZ,IX,IY) * G%DXDZ(IZ,IX,IY) / G%DYN(IZ,IX,IY)
   AS(IZ,IX,IY) = G%PERMONUS(IZ,IX,IY) * G%DXDZ(IZ,IX,IY) / G%DYS(IZ,IX,IY)
ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO 

! Set AP and B
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) DEFAULT(SHARED) COLLAPSE(3)
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ
   IF (G%IMASK(IZ,IX,IY)) THEN
      AP(IZ,IX,IY)=1D0
      B (IZ,IX,IY)=GPG%P0
      CYCLE
   ENDIF
   AP(IZ,IX,IY) = AP0(IZ,IX,IY) + AB(IZ,IX,IY) + AT(IZ,IX,IY) + &
                                    AE(IZ,IX,IY) + AW(IZ,IX,IY) + &
                                    AN(IZ,IX,IY) + AS(IZ,IX,IY)
   B (IZ,IX,IY) = AP0(IZ,IX,IY) * G%P(IZ,IX,IY) + &
                  (AP0(IZ,IX,IY) - AP0N(IZ,IX,IY)) * G%PN(IZ,IX,IY) + &
                  G%OMEGASFG(IZ,IX,IY) * G%DLTZN(IZ,IX,IY) * G%DLTX (IZ,IX,IY) * G%DLTY (IZ,IX,IY)

ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO

! Apply boundary conditions

DO ICOUNT = 1, NGPYRO_FACES_NEEDING_BCS(IMESH)
   IF (GP_BOUDARYS(IMESH)%IMESH_GPYRO(ICOUNT) .NE. IMESH) CYCLE
   IZ        = GP_BOUDARYS(IMESH)%IZ_GPYRO (ICOUNT)
   IX        = GP_BOUDARYS(IMESH)%IX_GPYRO (ICOUNT)
   IY        = GP_BOUDARYS(IMESH)%IY_GPYRO (ICOUNT)
   IOR       = GP_BOUDARYS(IMESH)%IOR_GPYRO(ICOUNT)

   IF (ABS(IOR)==3)  DS_OR = G%DXDY(IZ,IX,IY)
   IF (ABS(IOR)==1)  DS_OR = G%DYDZ(IZ,IX,IY)
   IF (ABS(IOR)==2)  DS_OR = G%DXDZ(IZ,IX,IY)

   IF (GPBCP(IZ,IX,IY,IOR)%PRES .LT. 0D0) THEN
      B(IZ,IX,IY) = B(IZ,IX,IY) + 1D-3*GPBCP(IZ,IX,IY,IOR)%MFLUX*DS_OR
   ELSE 
      CALL SET_BC_FIXED_VALUE(AP(IZ,IX,IY),AB(IZ,IX,IY),AT(IZ,IX,IY),AE(IZ,IX,IY),AW(IZ,IX,IY),AN(IZ,IX,IY),AS(IZ,IX,IY),B(IZ,IX,IY),GPBCP(IZ,IX,IY,IOR)%PRES)
   ENDIF
ENDDO



! Add terms for 2D/3D:
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) DEFAULT(SHARED) COLLAPSE(3)
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ
   BX(IZ,IX,IY) = 0D0; BY(IZ,IX,IY) = 0D0
   IF (G%IMASK(IZ,IX,IY)) CYCLE
   IF ((.NOT. G%NEEDSBCW(IZ,IX,IY)) .AND. IX .NE. 1     ) BX(IZ,IX,IY) = BX(IZ,IX,IY) + AW(IZ,IX,IY)*G%PN(IZ,IX-1,IY)
   IF ((.NOT. G%NEEDSBCE(IZ,IX,IY)) .AND. IX .NE. NCELLX) BX(IZ,IX,IY) = BX(IZ,IX,IY) + AE(IZ,IX,IY)*G%PN(IZ,IX+1,IY)
   IF ((.NOT. G%NEEDSBCS(IZ,IX,IY)) .AND. IY .NE. 1     ) BY(IZ,IX,IY) = BY(IZ,IX,IY) + AS(IZ,IX,IY)*G%PN(IZ,IX,IY-1)
   IF ((.NOT. G%NEEDSBCN(IZ,IX,IY)) .AND. IY .NE. NCELLY) BY(IZ,IX,IY) = BY(IZ,IX,IY) + AN(IZ,IX,IY)*G%PN(IZ,IX,IY+1)
ENDDO
ENDDO
ENDDO
!$omp END PARALLEL DO

PNOLD(:,:,:) = G%PN(:,:,:)
! z-direction sweep
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) DEFAULT(SHARED) COLLAPSE(3)
DO IX = 1, NCELLX
DO IY = 1, NCELLY
DO IZ = 1, NCELLZ
   AP(IZ,IX,IY) = AP(IZ,IX,IY) / ALPHA_P
   B (IZ,IX,IY) = B (IZ,IX,IY) + BX(IZ,IX,IY) + BY(IZ,IX,IY) + AP(IZ,IX,IY) * PNOLD(IZ,IX,IY) * (1D0 - ALPHA_P) 
ENDDO
ENDDO
ENDDO
!$omp END PARALLEL DO

PTR => G%PN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY) DEFAULT(SHARED) COLLAPSE(2)
DO IX = 1, NCELLX
DO IY = 1, NCELLY
   CALL TDMA_SOLVER_GENERAL(NCELLZ,PTR(:,IX,IY),AP(:,IX,IY),AB(:,IX,IY),AT(:,IX,IY),B(:,IX,IY))
ENDDO
ENDDO
!$OMP END PARALLEL DO


! Check convergence - RELATIVE TOLERANCE

G%CONV_INFO%CONVERGED_P   = .TRUE.
ALL_CELLS_CONVERGED       = .TRUE.
DPMAX = MAX(MAXVAL(ABS(PNOLD(:,:,:) - GPG%P0)), EPS_P)


IF (G%CONV_INFO%ITER .EQ. 1) THEN
   ALL_CELLS_CONVERGED = .FALSE.
ELSE

   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) DEFAULT(SHARED) COLLAPSE(3)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
      IF (G%IMASK(IZ,IX,IY)) THEN
         G%RESIDUAL_P(IZ,IX,IY)= 0D0
         CYCLE
      ENDIF   
      G%RESIDUAL_P(IZ,IX,IY) = ABS(G%PN(IZ,IX,IY) - PNOLD(IZ,IX,IY)) / DPMAX
      IF (G%RESIDUAL_P(IZ,IX,IY) .GT. GPG%PTOL*ALPHA_P) THEN
         ALL_CELLS_CONVERGED = .FALSE.
      ENDIF
      IF (G%PN(IZ,IX,IY).NE.G%PN(IZ,IX,IY) .OR. G%PN(IZ,IX,IY).EQ.GPG%POSINF .OR. G%PN(IZ,IX,IY).EQ.GPG%NEGINF) THEN
         GPG%NAN = .TRUE.
         ALL_CELLS_CONVERGED = .FALSE.
      ENDIF
   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDIF
IF (.NOT. ALL_CELLS_CONVERGED) THEN
   G%CONV_INFO%CONVERGED_P   = .FALSE.
   G%CONV_INFO%CONVERGED_ALL = .FALSE.
ENDIF
IF(G%CONV_INFO%CONVERGED_P) G%CONV_INFO%ITER_P = G%CONV_INFO%ITER

! *****************************************************************************
END SUBROUTINE PRESSURE_SOLVER
! *****************************************************************************

! *****************************************************************************
SUBROUTINE INSTANTANEOUS_GAS_RELASE(IMESH)
! *****************************************************************************
! Calcuates z-direction gas mass flux as instantaneous relase 

INTEGER, INTENT(IN) :: IMESH
INTEGER :: IZ_IN,IZ,IX,IY,IGSPEC,IOR,ICOUNT
INTEGER :: NCELLZ,NCELLX,NCELLY

G=>GPM(IMESH)
NCELLZ = G%NCELLZ
NCELLX = G%NCELLX
NCELLY = G%NCELLY

GPBCP(1:,1:,1:,-3:)=>G%GPYRO_BOUNDARY_CONDITION(1:,1:,1:,-3:)

G%MDOTPPZ(:,:,:,:) = 0D0

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ,IZ_IN,IOR,ICOUNT,IGSPEC) &
!$OMP DEFAULT(SHARED)
DO ICOUNT = 1, NGPYRO_FACES_NEEDING_BCS(IMESH)
   IOR = GP_BOUDARYS(IMESH)%IOR_GPYRO(ICOUNT)
   IF (IOR .NE. -3) CYCLE ! start at the cell at the bottom surface
   IZ_IN = GP_BOUDARYS(IMESH)%IZ_GPYRO (ICOUNT)
   IX    = GP_BOUDARYS(IMESH)%IX_GPYRO (ICOUNT)
   IY    = GP_BOUDARYS(IMESH)%IY_GPYRO (ICOUNT)

   ! Account for mass flux due to bottom boundary condtion
   IF (GPBCP(IZ_IN,IX,IY,IOR)%MFLUX .GT. 0D0) THEN
      DO IGSPEC = 1, GPROP%NGSPEC
         G%MDOTPPZ(IGSPEC,IZ_IN,IX,IY) = 1D-3*GPBCP(IZ_IN,IX,IY,IOR)%MFLUX * GPBCP(IZ_IN,IX,IY,IOR)%YJINF(IGSPEC) 
      ENDDO
   ENDIF

   ! Initialise flux from bottom
   DO IGSPEC = 1, GPROP%NGSPEC
      G%MDOTPPZ(IGSPEC,IZ_IN,IX,IY) = G%MDOTPPZ(IGSPEC,IZ_IN,IX,IY) + G%GOMEGA (3,IGSPEC,IZ_IN,IX,IY) * G%DLTZ(IZ_IN,IX,IY)
      G%MDOTPPZ(0     ,IZ_IN,IX,IY) = G%MDOTPPZ(0     ,IZ_IN,IX,IY) + G%MDOTPPZ(  IGSPEC,IZ_IN,IX,IY)
   ENDDO

   DO IZ = IZ_IN-1, 1, -1
      IF (G%IMASK(IZ,IX,IY)) EXIT
      DO IGSPEC = 1, GPROP%NGSPEC
         G%MDOTPPZ(IGSPEC,IZ,IX,IY) = G%MDOTPPZ(IGSPEC,IZ+1,IX,IY) + G%GOMEGA (3,IGSPEC,IZ,IX,IY) * G%DLTZ(IZ,IX,IY)
         G%MDOTPPZ(0     ,IZ,IX,IY) = G%MDOTPPZ(0     ,IZ  ,IX,IY) + G%MDOTPPZ(  IGSPEC,IZ,IX,IY) 
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

! *****************************************************************************
END SUBROUTINE INSTANTANEOUS_GAS_RELASE
! *****************************************************************************

! *****************************************************************************
SUBROUTINE DARCIAN_MASS_FLUX(IMESH,NCELLZ,NCELLX,NCELLY,DIRECTION)
! *****************************************************************************
! Calculate flow rate from pressure gradient

INTEGER, INTENT(IN) :: IMESH,NCELLZ,NCELLX,NCELLY
CHARACTER(1), INTENT(IN) :: DIRECTION
INTEGER :: IZ,IX,IY
REAL(EB) :: DPDZMRG,DPDXMRG,DPDYMRG


G=>GPM(IMESH)

SELECT CASE(DIRECTION) 

   CASE ('Z')
      !$OMP PARALLEL DO PRIVATE(IX,IY,IZ,DPDZMRG) SHARED(G,NCELLX,NCELLY,NCELLZ) SCHEDULE(STATIC) COLLAPSE(2)
      DO IY = 1, NCELLY
      DO IX = 1, NCELLX
         DO IZ = 1, NCELLZ-1 !Calculate z-direction mass flux (bottom)
            DPDZMRG = (G%PN(IZ+1,IX,IY) - G%PN(IZ,IX,IY)) / G%DZB(IZ,IX,IY) !- G%RGNB(IZ,IX,IY)*G%GZ  
            G%MDOTPPDARCYB(IZ,IX,IY) = -G%PERMONUB(IZ,IX,IY) * DPDZMRG
         ENDDO
         
         DO IZ = 2, NCELLZ !Set z-direction mass flux (top) from bottom values:
            G%MDOTPPDARCYT(IZ,IX,IY) = G%MDOTPPDARCYB(IZ-1,IX,IY)
         ENDDO
         
         !here, we still need the top value for cell 1 and the bottom value for cell NCELLZ
         !For now, make this approximation:
         G%MDOTPPDARCYT(1     ,IX,IY) = G%MDOTPPDARCYB(1     ,IX,IY)         
         G%MDOTPPDARCYB(NCELLZ,IX,IY) = G%MDOTPPDARCYT(NCELLZ,IX,IY)
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO
      
   CASE('X')
      !$OMP PARALLEL DO PRIVATE(IX,IY,IZ,DPDXMRG) SHARED(G,NCELLX,NCELLY,NCELLZ) SCHEDULE(STATIC) COLLAPSE(2)
      DO IY = 1, NCELLY
      DO IZ = 1, NCELLZ
         DO IX = 1, NCELLX-1 !Calculate x-direction mass flux (east)
            DPDXMRG = (G%PN(IZ,IX+1,IY) - G%PN(IZ,IX,IY)) / G%DXE(IZ,IX,IY) !- G%RGNE(IZ,IX,IY)*G%GX  
            G%MDOTPPDARCYE(IZ,IX,IY) = -G%PERMONUE(IZ,IX,IY) * DPDXMRG
         ENDDO
         
         DO IX = 2, NCELLX !Set x-direction mass flux (west) from east values:
            G%MDOTPPDARCYW(IZ,IX,IY) = G%MDOTPPDARCYE(IZ,IX-1,IY)
         ENDDO
         
         !here, we still need the west value for cell 1 and the east value for cell NCELLX
         !For now, make this approximation:
         G%MDOTPPDARCYW(IZ,1     ,IY) = G%MDOTPPDARCYE(IZ,1     ,IY)         
         G%MDOTPPDARCYE(IZ,NCELLX,IY) = G%MDOTPPDARCYW(IZ,NCELLX,IY)
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO

   CASE('Y')
      !$OMP PARALLEL DO PRIVATE(IX,IY,IZ,DPDYMRG) SHARED(G,NCELLX,NCELLY,NCELLZ) SCHEDULE(STATIC) COLLAPSE(2)
      DO IX = 1, NCELLX
      DO IZ = 1, NCELLZ
         DO IY = 1, NCELLY-1 !Calculate y-direction mass flux (north)
            DPDYMRG = (G%PN(IZ,IX,IY+1) - G%PN(IZ,IX,IY)) / G%DYN(IZ,IX,IY) !- G%RGNN(IZ,IX,IY)*G%GY  
            G%MDOTPPDARCYN(IZ,IX,IY) = -G%PERMONUN(IZ,IX,IY) * DPDYMRG
         ENDDO
         
         DO IY = 2, NCELLY !Set y-direction mass flux (south) from north values:
            G%MDOTPPDARCYS(IZ,IX,IY) = G%MDOTPPDARCYN(IZ,IX,IY-1)
         ENDDO
         
         !here, we still need the south value for cell 1 and the north value for cell NCELLY
         !For now, make this approximation:
         G%MDOTPPDARCYS(IZ,IX,1     ) = G%MDOTPPDARCYN(IZ,IX,1     )
         G%MDOTPPDARCYN(IZ,IX,NCELLY) = G%MDOTPPDARCYS(IZ,IX,NCELLY)
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO

END SELECT

! *****************************************************************************
END SUBROUTINE DARCIAN_MASS_FLUX
! *****************************************************************************

! *****************************************************************************
SUBROUTINE GET_MEAN_MOLECULARE_WEIGHT
! *****************************************************************************
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YJG
INTEGER :: NCELLZ,NCELLX,NCELLY,NGSPEC
INTEGER :: IZ,IX,IY,IGSPEC


NCELLZ = G%NCELLZ
NCELLX = G%NCELLX
NCELLY = G%NCELLY
NGSPEC = GPROP%NGSPEC

!The Gas masse fraction used is that at the last iteration
! as the new one is not calculated yet
YJG => G%YJG

IF (GPG%SOLVE_GAS_YJ) THEN
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ,IGSPEC) &
   !$OMP SHARED(YJG,G,NGSPEC,NCELLX,NCELLY,NCELLZ) COLLAPSE(3)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
      G%MN(IZ,IX,IY)=0D0
      IF (G%IMASK(IZ,IX,IY)) CYCLE
      DO IGSPEC = 1, NGSPEC
         G%MN(IZ,IX,IY)=G%MN(IZ,IX,IY)+ YJG(IGSPEC,IZ,IX,IY) * GPROP%M(IGSPEC)
      !MORN(IZ,IX,IY) = MORN(IZ,IX,IY) + G%YJGN(IGSPEC,IZ,IX,IY) * GPROP%M(IGSPEC) / 8314.
      !MOR (IZ,IX,IY) = MOR (IZ,IX,IY) + G%YJG (IGSPEC,IZ,IX,IY) * GPROP%M(IGSPEC) / 8314.
      ENDDO
   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ELSE
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) &
   !$OMP SHARED(GPROP,NCELLX,NCELLY,NCELLZ) COLLAPSE(3)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
      IF (G%IMASK(IZ,IX,IY)) CYCLE
      G%MN(IZ,IX,IY) = GPROP%M(1) 
      !MORN(IZ,IX,IY) = GPROP%M(1) / 8314.
      !MOR (IZ,IX,IY) = MORN (IZ,IX,IY)
   ENDDO
   ENDDO
   ENDDO

ENDIF
! *****************************************************************************
END SUBROUTINE GET_MEAN_MOLECULARE_WEIGHT
! *****************************************************************************











! *****************************************************************************
SUBROUTINE CONVECTIVE_DIFFUSIVE_SOLVER(IMESH,NGSPEC,ISCHEME,DTIME, &
           ALPHA_YJG,UPDATE_COEFFS)
! *****************************************************************************

INTEGER, INTENT(IN) :: IMESH,NGSPEC,ISCHEME
REAL(EB), INTENT(IN) :: DTIME, ALPHA_YJG
LOGICAL, INTENT(IN) :: UPDATE_COEFFS
LOGICAL :: SWEEPZ, SWEEPX, SWEEPY
REAL(EB) :: RESIDABS,OMALPHA_YJG,RALPHA_YJG,RDTIME,YJSUM,SOURCE
REAL(EB), PARAMETER :: EPS = 1D-6
REAL(EB), POINTER, DIMENSION(:,:,:) :: DT,DB,PT,PB,BIGAT,BIGAB,MDOTPPT,MDOTPPB,AP0,AP,AB,AT,B, &
                             MDOTPPE,MDOTPPW,DE,DW,PW,PE,AW,AE,BIGAW,BIGAE, &
                             AN,AS,MDOTPPN,MDOTPPS,DN,DS,PS,PN,BIGAS,BIGAN,AWIXW,AEIXE,ASIYS,ANIYN, &
                             ATIZT,ABIZB,PTR,BZ,BX,BY,BBC,AP1,APTOT,BTOT
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: YJGNOLD,AP0YJG,GOMEGA3DXDYDZ
LOGICAL, DIMENSION(0:NGSPEC) :: ALL_CELLS_CONVERGED
INTEGER :: I,IZ,IX,IY,IGSPEC,IOR
INTEGER:: NCELLZ,NCELLX,NCELLY



G=>GPM(IMESH)

NCELLX = G%NCELLX
NCELLY = G%NCELLY
NCELLZ = G%NCELLZ
GPBCP(1:,1:,1:,-3:)=>G%GPYRO_BOUNDARY_CONDITION(1:,1:,1:,-3:)
OMALPHA_YJG = 1D0 - ALPHA_YJG
RALPHA_YJG  = 1D0 / ALPHA_YJG
RDTIME      = 1D0 / DTIME

PT      => G%RWORK01
PB      => G%RWORK02
BIGAT   => G%RWORK03
BIGAB   => G%RWORK04
B       => G%RWORK05
MDOTPPE => G%RWORK06
MDOTPPW => G%RWORK07
DE      => G%RWORK08
DW      => G%RWORK09
PW      => G%RWORK10
PE      => G%RWORK11
BIGAW   => G%RWORK12
BIGAE   => G%RWORK13
MDOTPPN => G%RWORK14
MDOTPPS => G%RWORK15
DN      => G%RWORK16
DS      => G%RWORK17
PS      => G%RWORK18
PN      => G%RWORK19
BIGAS   => G%RWORK20
BIGAN   => G%RWORK21
BZ      => G%RWORK22
BX      => G%RWORK23
BY      => G%RWORK24
BBC     => G%RWORK25
APTOT   => G%RWORK26
BTOT    => G%RWORK27

AP0     => G%RWORK28
AP      => G%RWORK29
AB      => G%RWORK30
AT      => G%RWORK31
AW      => G%RWORK32
AE      => G%RWORK33
AN      => G%RWORK34
AS      => G%RWORK35

AWIXW   => G%RWORK36
AEIXE   => G%RWORK37
ASIYS   => G%RWORK38
ANIYN   => G%RWORK39
ATIZT   => G%RWORK40
ABIZB   => G%RWORK41
AP1     => G%RWORK42

DT      => G%RWORK43
DB      => G%RWORK44
MDOTPPT => G%RWORK45
MDOTPPB => G%RWORK46

YJGNOLD         => G%RWORK110
AP0YJG          => G%RWORK111
GOMEGA3DXDYDZ   => G%RWORK112

! Vary sweep direction across iterations 
CALL SET_SWEEP_DIRECTION(SWEEPZ,SWEEPX,SWEEPY)

! Override this and just sweep in z:
SWEEPZ=.TRUE.
SWEEPX=.FALSE.
SWEEPY=.FALSE.

IF (GPG%SOLVE_PRESSURE) THEN
   MDOTPPT(:,:,:) =  G%MDOTPPDARCYT(:,:,:)
   MDOTPPB(:,:,:) =  G%MDOTPPDARCYB(:,:,:)

   IF (NCELLX .GT. 1) THEN
      MDOTPPE(:,:,:) = G%MDOTPPDARCYE(:,:,:)
      MDOTPPW(:,:,:) = G%MDOTPPDARCYW(:,:,:)
   ENDIF

   IF (NCELLX .GT. 1 .AND. NCELLY .GT. 1) THEN
      MDOTPPN(:,:,:) = G%MDOTPPDARCYN(:,:,:)
      MDOTPPS(:,:,:) = G%MDOTPPDARCYS(:,:,:)
   ENDIF

ELSE
   MDOTPPT(:,:,:) = -G%MDOTPPZ(0,:,:,:)
   DO IZ = 1, NCELLZ-1
      MDOTPPB(IZ,:,:) = MDOTPPT(IZ+1,:,:)
   ENDDO   
   MDOTPPB(NCELLZ,:,:) = MDOTPPB(NCELLZ-1,:,:) !THIS APPROXIMATION SHOULD BE FIXED

   IF (NCELLX .GT. 1) THEN
      MDOTPPE(:,:,:) = 0D0
      MDOTPPW(:,:,:) = 0D0
   ENDIF

   IF (NCELLX .GT. 1 .AND. NCELLY .GT. 1) THEN
      MDOTPPN(:,:,:) = 0D0
      MDOTPPS(:,:,:) = 0D0
   ENDIF
ENDIF         

! On first iteration, calculate quantities that depend on old values (AP0 and AP0*Yjg)
IF (G%CONV_INFO%ITER .EQ. 1) THEN
   !$OMP PARALLEL DO PRIVATE(IX,IY,IZ,IGSPEC) SHARED(G,AP,AB,AT,AE,AW,AN,AS,AP0,AP0YJG,RDTIME,NGSPEC,NCELLX,NCELLY,NCELLZ) SCHEDULE(STATIC) COLLAPSE(3)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
      IF (G%IMASK(IZ,IX,IY)) THEN
         AP(IZ,IX,IY) = 1D0
         AB(IZ,IX,IY) = 0D0 ; AT(IZ,IX,IY) = 0D0 
         AE(IZ,IX,IY) = 0D0 ; AW(IZ,IX,IY) = 0D0 
         AN(IZ,IX,IY) = 0D0 ; AS(IZ,IX,IY) = 0D0 
      ELSE
         AP0(IZ,IX,IY) = G%RG(IZ,IX,IY)*G%DV(IZ,IX,IY)*G%POROSS(IZ,IX,IY) * RDTIME
         DO IGSPEC = 1, NGSPEC
            AP0YJG(IGSPEC,IZ,IX,IY) = AP0(IZ,IX,IY)*G%YJG(IGSPEC,IZ,IX,IY)
         ENDDO
      ENDIF
   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDIF

! Set the effective front and back face diffusion coefficients
! Note that HM0 is set appropriately in subroutine GET_ALL_BOUNDARY_CONDITIONS

IF (NCELLZ .GT. 1 .AND. NCELLX .EQ. 1 .AND. NCELLY .EQ. 1) THEN
   IX=1; IY=1
   IF (UPDATE_COEFFS) THEN
      !$OMP PARALLEL DO PRIVATE(IZ,IGSPEC,IOR) DEFAULT(SHARED) SCHEDULE(STATIC)
      DO IZ = 1, NCELLZ
         IF (G%IMASK(IZ,IX,IY)) CYCLE

         ! Calculate "D", equal to diffusion coefficient divided by grid spacing
         DT(IZ,IX,IY) = G%PSIRGDT(IZ,IX,IY) / G%DZT(IZ,IX,IY)
         DB(IZ,IX,IY) = G%PSIRGDB(IZ,IX,IY) / G%DZB(IZ,IX,IY)
         IF (G%NEEDSBCT(IZ,IX,IY)) THEN
            IOR=3
            IF (GPG%USE_TORTUOSITY_FACTOR_FOR_FLUX) THEN
               DT(IZ,IX,IY) = GPG%TORTUOSITY_FACTOR*GPBCP(IZ,IX,IY,IOR)%HM0 !Set front-face boundary condition
            ELSE
               DT(IZ,IX,IY) = GPBCP(IZ,IX,IY,IOR)%HM0 !Set front-face boundary condition
            ENDIF
         ENDIF
         IF (G%NEEDSBCB(IZ,IX,IY)) THEN
            IOR=-3
            IF (GPG%USE_TORTUOSITY_FACTOR_FOR_FLUX) THEN
               DB(IZ,IX,IY) = GPG%TORTUOSITY_FACTOR*GPBCP(IZ,IX,IY,IOR)%HM0 !Set back-face boundary condition
            ELSE
               DB(IZ,IX,IY) = GPBCP(IZ,IX,IY, IOR)%HM0 !Set back-face boundary condition
            ENDIF
         ENDIF
         IF (DT(IZ,IX,IY) .LT. 1D-50) DT(IZ,IX,IY) = 1D-50
         IF (DB(IZ,IX,IY) .LT. 1D-50) DB(IZ,IX,IY) = 1D-50

         ! Calculate cell Peclet number
         PT(IZ,IX,IY) = MDOTPPT(IZ,IX,IY) / DT(IZ,IX,IY)
         PB(IZ,IX,IY) = MDOTPPB(IZ,IX,IY) / DB(IZ,IX,IY)
         CALL GET_A_SINGLE(ISCHEME,PT(IZ,IX,IY),PB(IZ,IX,IY),BIGAT(IZ,IX,IY),BIGAB(IZ,IX,IY)) ! Calculate A(|P|)
         AT(IZ,IX,IY) = (DT(IZ,IX,IY) * BIGAT(IZ,IX,IY) + MAX ( MDOTPPT(IZ,IX,IY), 0D0)) * G%DXDY(IZ,IX,IY)
         AB(IZ,IX,IY) = (DB(IZ,IX,IY) * BIGAB(IZ,IX,IY) + MAX (-MDOTPPB(IZ,IX,IY), 0D0)) * G%DXDY(IZ,IX,IY)

         IF (G%NEEDSBCT(IZ,IX,IY)) THEN
            ATIZT(IZ,IX,IY) = AT(IZ,IX,IY)
            AT   (IZ,IX,IY) = 0D0
         ENDIF

         IF (G%NEEDSBCB(IZ,IX,IY)) THEN
            ABIZB(IZ,IX,IY) = AB(IZ,IX,IY)
            AB   (IZ,IX,IY) = 0D0
         ENDIF
      ENDDO
      !$OMP END PARALLEL DO
   ENDIF !UPDATE_COEFFS

   !$OMP PARALLEL DO PRIVATE(IZ,IGSPEC) SHARED(NCELLZ,AP1,AP0,G,AT,AB,NGSPEC,GOMEGA3DXDYDZ) SCHEDULE(STATIC)
   DO IZ = 1, NCELLZ
      IF (G%IMASK(IZ,IX,IY)) CYCLE
      AP1(IZ,IX,IY) = AP0(IZ,IX,IY) + G%OMEGASFG(IZ,IX,IY)*G%DV(IZ,IX,IY) + AT(IZ,IX,IY) + AB(IZ,IX,IY)
      DO IGSPEC = 1, NGSPEC
         GOMEGA3DXDYDZ(IGSPEC,IZ,IX,IY) = G%GOMEGA(3,IGSPEC,IZ,IX,IY) * G%DV(IZ,IX,IY)
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDIF

IF (NCELLZ .GT. 1 .AND. NCELLX .GT. 1 .AND. NCELLY .EQ. 1) THEN
   IY = 1
   IF (UPDATE_COEFFS) THEN
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IOR) DEFAULT(SHARED) COLLAPSE(2)
      DO IX = 1, NCELLX
      DO IZ = 1, NCELLZ

         IF (G%IMASK(IZ,IX,IY)) CYCLE

         ! Calculate "D", equal to diffusion coefficient divided by grid spacing
         DT(IZ,IX,IY) = G%PSIRGDT(IZ,IX,IY) / G%DZT(IZ,IX,IY)
         DB(IZ,IX,IY) = G%PSIRGDB(IZ,IX,IY) / G%DZB(IZ,IX,IY)
         IF (G%NEEDSBCT(IZ,IX,IY)) THEN
            IOR=3
            IF (GPG%USE_TORTUOSITY_FACTOR_FOR_FLUX) THEN
               DT(IZ,IX,IY) = GPG%TORTUOSITY_FACTOR*GPBCP(IZ,IX,IY,IOR)%HM0 !Set front-face boundary condition
            ELSE
               DT(IZ,IX,IY) = GPBCP(IZ,IX,IY,IOR)%HM0 !Set front-face boundary condition
            ENDIF
         ENDIF
         IF (G%NEEDSBCB(IZ,IX,IY)) THEN
            IOR=-3
            IF (GPG%USE_TORTUOSITY_FACTOR_FOR_FLUX) THEN
               DB(IZ,IX,IY) = GPG%TORTUOSITY_FACTOR*GPBCP(IZ,IX,IY,IOR)%HM0 !Set back-face boundary condition
            ELSE
               DB(IZ,IX,IY) = GPBCP(IZ,IX,IY,IOR)%HM0 !Set back-face boundary condition
            ENDIF
         ENDIF
         IF (DT(IZ,IX,IY) .LT. 1D-50) DT(IZ,IX,IY) = 1D-50
         IF (DB(IZ,IX,IY) .LT. 1D-50) DB(IZ,IX,IY) = 1D-50

         ! Calculate cell Peclet number
         PT(IZ,IX,IY) = MDOTPPT(IZ,IX,IY) / DT(IZ,IX,IY)
         PB(IZ,IX,IY) = MDOTPPB(IZ,IX,IY) / DB(IZ,IX,IY)
         CALL GET_A_SINGLE(ISCHEME,PT(IZ,IX,IY),PB(IZ,IX,IY),BIGAT(IZ,IX,IY),BIGAB(IZ,IX,IY)) ! Calculate A(|P|)
         AT(IZ,IX,IY) = (DT(IZ,IX,IY) * BIGAT(IZ,IX,IY) + MAX ( MDOTPPT(IZ,IX,IY), 0D0)) * G%DXDY(IZ,IX,IY)
         AB(IZ,IX,IY) = (DB(IZ,IX,IY) * BIGAB(IZ,IX,IY) + MAX (-MDOTPPB(IZ,IX,IY), 0D0)) * G%DXDY(IZ,IX,IY)

         IF (G%NEEDSBCT(IZ,IX,IY)) THEN
            ATIZT(IZ,IX,IY) = AT(IZ,IX,IY)
            AT   (IZ,IX,IY) = 0D0
         ENDIF

         IF (G%NEEDSBCB(IZ,IX,IY)) THEN
            ABIZB(IZ,IX,IY) = AB(IZ,IX,IY)
            AB   (IZ,IX,IY) = 0D0
         ENDIF

         ! Calculate "D", equal to diffusion coefficient divided by grid spacing
         DW(IZ,IX,IY) = G%PSIRGDW(IZ,IX,IY) / G%DXW(IZ,IX,IY)
         DE(IZ,IX,IY) = G%PSIRGDE(IZ,IX,IY) / G%DXE(IZ,IX,IY)
         IF (G%NEEDSBCW(IZ,IX,IY)) THEN
            IOR=-1
            IF (GPG%USE_TORTUOSITY_FACTOR_FOR_FLUX) THEN
               DW(IZ,IX,IY) = GPG%TORTUOSITY_FACTOR*GPBCP(IZ,IX,IY,IOR)%HM0
            ELSE
               DW(IZ,IX,IY) = GPBCP(IZ,IX,IY,IOR)%HM0
            ENDIF
         ENDIF
         IF (G%NEEDSBCE(IZ,IX,IY)) THEN
            IOR= 1
            IF (GPG%USE_TORTUOSITY_FACTOR_FOR_FLUX) THEN
               DE(IZ,IX,IY) = GPG%TORTUOSITY_FACTOR*GPBCP(IZ,IX,IY,IOR)%HM0
            ELSE
               DE(IZ,IX,IY) = GPBCP(IZ,IX,IY,IOR)%HM0
            ENDIF
         ENDIF
         IF (DW(IZ,IX,IY) .LT. 1D-50) DW(IZ,IX,IY) = 1D-50
         IF (DE(IZ,IX,IY) .LT. 1D-50) DE(IZ,IX,IY) = 1D-50

         ! Calculate cell Peclet number
         PW(IZ,IX,IY) = MDOTPPW(IZ,IX,IY) / DW(IZ,IX,IY)
         PE(IZ,IX,IY) = MDOTPPE(IZ,IX,IY) / DE(IZ,IX,IY)
         CALL GET_A_SINGLE(ISCHEME,PW(IZ,IX,IY),PE(IZ,IX,IY),BIGAW(IZ,IX,IY),BIGAE(IZ,IX,IY)) ! Calculate A(|P|)
         AW(IZ,IX,IY) = (DW(IZ,IX,IY) * BIGAW(IZ,IX,IY) + MAX ( MDOTPPW(IZ,IX,IY), 0D0) ) * G%DYDZ(IZ,IX,IY)
         AE(IZ,IX,IY) = (DE(IZ,IX,IY) * BIGAE(IZ,IX,IY) + MAX (-MDOTPPE(IZ,IX,IY), 0D0) ) * G%DYDZ(IZ,IX,IY)

         IF (G%NEEDSBCW(IZ,IX,IY)) THEN 
            AWIXW(IZ,IX,IY) = AW(IZ,IX,IY)
            AW(IZ,IX,IY)    = 0D0
         ENDIF

         IF (G%NEEDSBCE(IZ,IX,IY)) THEN
            AEIXE(IZ,IX,IY) = AE(IZ,IX,IY)
            AE(IZ,IX,IY)    = 0D0
         ENDIF

      ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   ENDIF !UPDATE_COEFFS

   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IGSPEC) &
   !$OMP SHARED(G,AP1,AP0,AT,AB,AE,AW,GOMEGA3DXDYDZ,NCELLX,NCELLZ,NGSPEC) COLLAPSE(2)
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
      IF (G%IMASK(IZ,IX,IY)) CYCLE

      AP1(IZ,IX,IY) = AP0(IZ,IX,IY) + G%OMEGASFG(IZ,IX,IY)*G%DV(IZ,IX,IY) + AT(IZ,IX,IY) + AB(IZ,IX,IY) + AE(IZ,IX,IY) + AW(IZ,IX,IY)
      DO IGSPEC = 1, NGSPEC
         GOMEGA3DXDYDZ(IGSPEC,IZ,IX,IY) = G%GOMEGA(3,IGSPEC,IZ,IX,IY) * G%DV(IZ,IX,IY)
      ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

ENDIF

IF (NCELLZ .GT. 1 .AND. NCELLX .GT. 1 .AND. NCELLY .GT. 1) THEN
   IF (UPDATE_COEFFS) THEN      
      !$OMP PARALLEL DO PRIVATE(IX,IY,IZ,IGSPEC,IOR) DEFAULT(SHARED) SCHEDULE(STATIC) COLLAPSE(3)
      DO IY = 1, NCELLY
      DO IX = 1, NCELLX
      DO IZ = 1, NCELLZ

         IF (G%IMASK(IZ,IX,IY)) CYCLE

         ! Calculate "D", equal to diffusion coefficient divided by grid spacing
         DT(IZ,IX,IY) = G%PSIRGDT(IZ,IX,IY) / G%DZT(IZ,IX,IY)
         DB(IZ,IX,IY) = G%PSIRGDB(IZ,IX,IY) / G%DZB(IZ,IX,IY)
         IF (G%NEEDSBCT(IZ,IX,IY)) THEN
            IOR= 3
            IF (GPG%USE_TORTUOSITY_FACTOR_FOR_FLUX) THEN
               DT(IZ,IX,IY) = GPG%TORTUOSITY_FACTOR*GPBCP(IZ,IX,IY,IOR)%HM0 !Set front-face boundary condition
            ELSE
               DT(IZ,IX,IY) = GPBCP(IZ,IX,IY,IOR)%HM0 !Set front-face boundary condition
            ENDIF
         ENDIF
         IF (G%NEEDSBCB(IZ,IX,IY)) THEN
            IOR=-3
            IF (GPG%USE_TORTUOSITY_FACTOR_FOR_FLUX) THEN
               DB(IZ,IX,IY) = GPG%TORTUOSITY_FACTOR*GPBCP(IZ,IX,IY,IOR)%HM0 !Set back-face boundary condition
            ELSE
               DB(IZ,IX,IY) = GPBCP(IZ,IX,IY,IOR)%HM0 !Set back-face boundary condition
            ENDIF
         ENDIF
         IF (DT(IZ,IX,IY) .LT. 1D-50) DT(IZ,IX,IY) = 1D-50
         IF (DB(IZ,IX,IY) .LT. 1D-50) DB(IZ,IX,IY) = 1D-50

         ! Calculate cell Peclet number
         PT(IZ,IX,IY) = MDOTPPT(IZ,IX,IY) / DT(IZ,IX,IY)
         PB(IZ,IX,IY) = MDOTPPB(IZ,IX,IY) / DB(IZ,IX,IY)
         CALL GET_A_SINGLE(ISCHEME,PT(IZ,IX,IY),PB(IZ,IX,IY),BIGAT(IZ,IX,IY),BIGAB(IZ,IX,IY)) ! Calculate A(|P|)
         AT(IZ,IX,IY) = (DT(IZ,IX,IY) * BIGAT(IZ,IX,IY) + MAX ( MDOTPPT(IZ,IX,IY), 0D0)) * G%DXDY(IZ,IX,IY)
         AB(IZ,IX,IY) = (DB(IZ,IX,IY) * BIGAB(IZ,IX,IY) + MAX (-MDOTPPB(IZ,IX,IY), 0D0)) * G%DXDY(IZ,IX,IY)

         IF (G%NEEDSBCT(IZ,IX,IY)) THEN
            ATIZT(IZ,IX,IY) = AT(IZ,IX,IY)
            AT   (IZ,IX,IY) = 0D0
         ENDIF

         IF (G%NEEDSBCB(IZ,IX,IY)) THEN
            ABIZB(IZ,IX,IY) = AB(IZ,IX,IY)
            AB   (IZ,IX,IY) = 0D0
         ENDIF

         ! Calculate "D", equal to diffusion coefficient divided by grid spacing
         DW(IZ,IX,IY) = G%PSIRGDW(IZ,IX,IY) / G%DXW(IZ,IX,IY)
         DE(IZ,IX,IY) = G%PSIRGDE(IZ,IX,IY) / G%DXE(IZ,IX,IY)
         IF (G%NEEDSBCW(IZ,IX,IY)) THEN
            IOR=-1
            IF (GPG%USE_TORTUOSITY_FACTOR_FOR_FLUX) THEN
               DW(IZ,IX,IY) = GPG%TORTUOSITY_FACTOR*GPBCP(IZ,IX,IY,IOR)%HM0
            ELSE
               DW(IZ,IX,IY) = GPBCP(IZ,IX,IY,IOR)%HM0
            ENDIF
         ENDIF
         IF (G%NEEDSBCE(IZ,IX,IY)) THEN
            IOR= 1
            IF (GPG%USE_TORTUOSITY_FACTOR_FOR_FLUX) THEN
               DE(IZ,IX,IY) = GPG%TORTUOSITY_FACTOR*GPBCP(IZ,IX,IY,IOR)%HM0
            ELSE
               DE(IZ,IX,IY) = GPBCP(IZ,IX,IY,IOR)%HM0
            ENDIF
         ENDIF
         IF (DW(IZ,IX,IY) .LT. 1D-50) DW(IZ,IX,IY) = 1D-50
         IF (DE(IZ,IX,IY) .LT. 1D-50) DE(IZ,IX,IY) = 1D-50

         ! Calculate cell Peclet number
         PW(IZ,IX,IY) = MDOTPPW(IZ,IX,IY) / DW(IZ,IX,IY)
         PE(IZ,IX,IY) = MDOTPPE(IZ,IX,IY) / DE(IZ,IX,IY)
         CALL GET_A_SINGLE(ISCHEME,PW(IZ,IX,IY),PE(IZ,IX,IY),BIGAW(IZ,IX,IY),BIGAE(IZ,IX,IY)) ! Calculate A(|P|)
         AW(IZ,IX,IY) = (DW(IZ,IX,IY) * BIGAW(IZ,IX,IY) + MAX ( MDOTPPW(IZ,IX,IY), 0D0) ) * G%DYDZ(IZ,IX,IY)
         AE(IZ,IX,IY) = (DE(IZ,IX,IY) * BIGAE(IZ,IX,IY) + MAX (-MDOTPPE(IZ,IX,IY), 0D0) ) * G%DYDZ(IZ,IX,IY)

         IF (G%NEEDSBCW(IZ,IX,IY)) THEN 
            AWIXW(IZ,IX,IY) = AW(IZ,IX,IY)
            AW(IZ,IX,IY)    = 0D0
         ENDIF

         IF (G%NEEDSBCE(IZ,IX,IY)) THEN
            AEIXE(IZ,IX,IY) = AE(IZ,IX,IY)
            AE(IZ,IX,IY)    = 0D0
         ENDIF

         DS(IZ,IX,IY) = G%PSIRGDS(IZ,IX,IY) / G%DYS(IZ,IX,IY)
         DN(IZ,IX,IY) = G%PSIRGDN(IZ,IX,IY) / G%DYN(IZ,IX,IY)
         IF (G%NEEDSBCS(IZ,IX,IY)) THEN
            IOR=-2
            IF (GPG%USE_TORTUOSITY_FACTOR_FOR_FLUX) THEN
               DS(IZ,IX,IY) = GPG%TORTUOSITY_FACTOR*GPBCP(IZ,IX,IY,IOR)%HM0
            ELSE
               DS(IZ,IX,IY) = GPBCP(IZ,IX,IY,IOR)%HM0
            ENDIF
         ENDIF
         IF (G%NEEDSBCN(IZ,IX,IY)) THEN
            IOR= 2
            IF (GPG%USE_TORTUOSITY_FACTOR_FOR_FLUX) THEN
               DN(IZ,IX,IY) = GPG%TORTUOSITY_FACTOR*GPBCP(IZ,IX,IY,IOR)%HM0
            ELSE
               DN(IZ,IX,IY) = GPBCP(IZ,IX,IY,IOR)%HM0
            ENDIF
         ENDIF
         IF (DS(IZ,IX,IY) .LT. 1D-50) DS(IZ,IX,IY) = 1D-50
         IF (DN(IZ,IX,IY) .LT. 1D-50) DN(IZ,IX,IY) = 1D-50

         ! Calculate cell Peclet number
         PS(IZ,IX,IY) = MDOTPPS(IZ,IX,IY) / DS(IZ,IX,IY)
         PN(IZ,IX,IY) = MDOTPPN(IZ,IX,IY) / DN(IZ,IX,IY)
         CALL GET_A_SINGLE(ISCHEME,PS(IZ,IX,IY),PN(IZ,IX,IY),BIGAS(IZ,IX,IY),BIGAN(IZ,IX,IY)) ! Calculate A(|P|)
         AS(IZ,IX,IY) = (DS(IZ,IX,IY) * BIGAS(IZ,IX,IY) + MAX ( MDOTPPS(IZ,IX,IY), 0D0) ) * G%DXDZ(IZ,IX,IY)
         AN(IZ,IX,IY) = (DN(IZ,IX,IY) * BIGAN(IZ,IX,IY) + MAX (-MDOTPPN(IZ,IX,IY), 0D0) ) * G%DXDZ(IZ,IX,IY)

         IF (G%NEEDSBCS(IZ,IX,IY)) THEN 
            ASIYS(IZ,IX,IY)=AS(IZ,IX,IY)
            AS   (IZ,IX,IY)=0D0
         ENDIF

         IF (G%NEEDSBCN(IZ,IX,IY)) THEN
            ANIYN(IZ,IX,IY)=AN(IZ,IX,IY)
            AN   (IZ,IX,IY)=0D0
         ENDIF

      ENDDO
      ENDDO
      ENDDO
      !$omp END PARALLEL DO
   ENDIF !UPDATE_COEFFS

   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ,IGSPEC) SHARED(G,AP1,AP0,AT,AB,AE,AW,AN,AS,NCELLX,NCELLY,NCELLZ,NGSPEC,GOMEGA3DXDYDZ) COLLAPSE(2)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
      IF (G%IMASK(IZ,IX,IY)) CYCLE
      AP1(IZ,IX,IY) = AP0(IZ,IX,IY) + G%OMEGASFG(IZ,IX,IY)*G%DV(IZ,IX,IY) + AT(IZ,IX,IY) + AB(IZ,IX,IY) + AE(IZ,IX,IY) + AW(IZ,IX,IY) + AN(IZ,IX,IY) + AS(IZ,IX,IY)
      DO IGSPEC = 1, NGSPEC
         GOMEGA3DXDYDZ(IGSPEC,IZ,IX,IY) = G%GOMEGA(3,IGSPEC,IZ,IX,IY) * G%DV(IZ,IX,IY)
      ENDDO
   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDIF

YJGNOLD(:,:,:,:) = G%YJGN(:,:,:,:)

! Loop over all gaseous species
DO IGSPEC = 1, GPROP%NGSPEC

   !IF (ITER .GT. 1 .AND. MAXVAL(G%RESIDUAL_YJG(IGSPEC,:,:,:)*REAL(G%IMASK(:,:,:)))*RALPHA_YJG .LT. 1D-30) CYCLE

   PTR => G%YJGN(IGSPEC,:,:,:) !This is to avoid creating a temporary array in TDMA_SOLVER_GENERAL

   ! Apply boundary conditions:
   IF (NCELLZ .GT. 1 .AND. NCELLX .EQ. 1 .AND. NCELLY .EQ. 1) THEN !1D
      IX=1; IY=1
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IZ) &
      !$OMP SHARED(IGSPEC,IX,IY,BBC,ATIZT,ABIZB,GPBCP,G,NCELLZ)
      DO IZ = 1, NCELLZ
         IF (G%IMASK(IZ,IX,IY)) CYCLE
         BBC(IZ,IX,IY) = 0D0
         IF (G%NEEDSBCT(IZ,IX,IY)) BBC (IZ,IX,IY) = BBC (IZ,IX,IY) + ATIZT(IZ,IX,IY)*(GPBCP(IZ,IX,IY, 3)%YJINF(IGSPEC)-G%YJGN(IGSPEC,IZ,IX,IY))
         IF (G%NEEDSBCB(IZ,IX,IY)) BBC (IZ,IX,IY) = BBC (IZ,IX,IY) + ABIZB(IZ,IX,IY)*(GPBCP(IZ,IX,IY,-3)%YJINF(IGSPEC)-G%YJGN(IGSPEC,IZ,IX,IY))
         !IF (G%NEEDSBCT(IZ,IX,IY)) BBC (IZ,IX,IY) = BBC (IZ,IX,IY) + ATIZT(IZ,IX,IY)*(GPBCP(IZ,IX,IY, 3)%YJINF(IGSPEC)-G%YJG (IGSPEC,IZ,IX,IY))
         !IF (G%NEEDSBCB(IZ,IX,IY)) BBC (IZ,IX,IY) = BBC (IZ,IX,IY) + ABIZB(IZ,IX,IY)*(GPBCP(IZ,IX,IY,-3)%YJINF(IGSPEC)-G%YJG (IGSPEC,IZ,IX,IY))

      ENDDO
      !$OMP END PARALLEL DO
   ENDIF

   IF (NCELLZ .GT. 1 .AND. NCELLX .GT. 1 .AND. NCELLY .EQ. 1) THEN !2D
      IY=1
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IZ) &
      !$OMP SHARED(IGSPEC,IY,BBC,ATIZT,ABIZB,AWIXW,AEIXE,GPBCP,G,NCELLX,NCELLZ) COLLAPSE(2)
      DO IX = 1, NCELLX
      DO IZ = 1, NCELLZ
         IF (G%IMASK(IZ,IX,IY)) CYCLE
         BBC(IZ,IX,IY) = 0D0
         IF (G%NEEDSBCT(IZ,IX,IY)) BBC (IZ,IX,IY) = BBC (IZ,IX,IY) + ATIZT(IZ,IX,IY)*(GPBCP(IZ,IX,IY, 3)%YJINF(IGSPEC)-G%YJGN(IGSPEC,IZ,IX,IY))
         IF (G%NEEDSBCB(IZ,IX,IY)) BBC (IZ,IX,IY) = BBC (IZ,IX,IY) + ABIZB(IZ,IX,IY)*(GPBCP(IZ,IX,IY,-3)%YJINF(IGSPEC)-G%YJGN(IGSPEC,IZ,IX,IY))
         IF (G%NEEDSBCW(IZ,IX,IY)) BBC (IZ,IX,IY) = BBC (IZ,IX,IY) + AWIXW(IZ,IX,IY)*(GPBCP(IZ,IX,IY,-1)%YJINF(IGSPEC)-G%YJGN(IGSPEC,IZ,IX,IY))
         IF (G%NEEDSBCE(IZ,IX,IY)) BBC (IZ,IX,IY) = BBC (IZ,IX,IY) + AEIXE(IZ,IX,IY)*(GPBCP(IZ,IX,IY, 1)%YJINF(IGSPEC)-G%YJGN(IGSPEC,IZ,IX,IY))
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   ENDIF

   IF (NCELLZ .GT. 1 .AND. NCELLX .GT. 1 .AND. NCELLY .GT. 1) THEN !3D
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) &
      !$OMP SHARED(IGSPEC,BBC,ATIZT,ABIZB,AWIXW,AEIXE,ASIYS,ANIYN,GPBCP,G,NCELLX,NCELLY,NCELLZ) COLLAPSE(3)
      DO IY = 1, NCELLY
      DO IX = 1, NCELLX
      DO IZ = 1, NCELLZ
         IF (G%IMASK(IZ,IX,IY)) CYCLE
         BBC(IZ,IX,IY) = 0D0
         IF (G%NEEDSBCT(IZ,IX,IY)) BBC (IZ,IX,IY) = BBC (IZ,IX,IY) + ATIZT(IZ,IX,IY)*(GPBCP(IZ,IX,IY, 3)%YJINF(IGSPEC)-G%YJGN(IGSPEC,IZ,IX,IY))
         IF (G%NEEDSBCB(IZ,IX,IY)) BBC (IZ,IX,IY) = BBC (IZ,IX,IY) + ABIZB(IZ,IX,IY)*(GPBCP(IZ,IX,IY,-3)%YJINF(IGSPEC)-G%YJGN(IGSPEC,IZ,IX,IY))
         IF (G%NEEDSBCW(IZ,IX,IY)) BBC (IZ,IX,IY) = BBC (IZ,IX,IY) + AWIXW(IZ,IX,IY)*(GPBCP(IZ,IX,IY,-1)%YJINF(IGSPEC)-G%YJGN(IGSPEC,IZ,IX,IY))
         IF (G%NEEDSBCE(IZ,IX,IY)) BBC (IZ,IX,IY) = BBC (IZ,IX,IY) + AEIXE(IZ,IX,IY)*(GPBCP(IZ,IX,IY, 1)%YJINF(IGSPEC)-G%YJGN(IGSPEC,IZ,IX,IY))
         IF (G%NEEDSBCS(IZ,IX,IY)) BBC (IZ,IX,IY) = BBC (IZ,IX,IY) + ASIYS(IZ,IX,IY)*(GPBCP(IZ,IX,IY,-2)%YJINF(IGSPEC)-G%YJGN(IGSPEC,IZ,IX,IY))
         IF (G%NEEDSBCN(IZ,IX,IY)) BBC (IZ,IX,IY) = BBC (IZ,IX,IY) + ANIYN(IZ,IX,IY)*(GPBCP(IZ,IX,IY, 2)%YJINF(IGSPEC)-G%YJGN(IGSPEC,IZ,IX,IY))
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   ENDIF


   ! Calculate diffusion terms:
   IF (NCELLZ .GT. 1 .AND. NCELLX .GT. 1 .AND. NCELLY .EQ. 1) THEN !2D
      IY=1

      IF (SWEEPZ .AND. SWEEPX) THEN
         DO IX = 1, NCELLX
         DO IZ = 1, NCELLZ
            IF (G%IMASK(IZ,IX,IY)) CYCLE
            BZ (IZ,IX,IY) = 0D0; BX (IZ,IX,IY) = 0D0
            IF ( (.NOT. G%NEEDSBCT(IZ,IX,IY)) .AND. IZ .NE. 1     ) BZ (IZ,IX,IY) = BZ (IZ,IX,IY) + AT(IZ,IX,IY) * G%YJGN(IGSPEC,IZ-1,IX,IY)
            IF ( (.NOT. G%NEEDSBCB(IZ,IX,IY)) .AND. IZ .NE. NCELLZ) BZ (IZ,IX,IY) = BZ (IZ,IX,IY) + AB(IZ,IX,IY) * G%YJGN(IGSPEC,IZ+1,IX,IY)
            IF ( (.NOT. G%NEEDSBCW(IZ,IX,IY)) .AND. IX .NE. 1     ) BX (IZ,IX,IY) = BX (IZ,IX,IY) + AW(IZ,IX,IY) * G%YJGN(IGSPEC,IZ,IX-1,IY)
            IF ( (.NOT. G%NEEDSBCE(IZ,IX,IY)) .AND. IX .NE. NCELLX) BX (IZ,IX,IY) = BX (IZ,IX,IY) + AE(IZ,IX,IY) * G%YJGN(IGSPEC,IZ,IX+1,IY)
         ENDDO
         ENDDO
      ENDIF

      IF (SWEEPZ .AND. (.NOT. SWEEPX)) THEN
         DO IX = 1, NCELLX
         DO IZ = 1, NCELLZ
            IF (G%IMASK(IZ,IX,IY)) CYCLE
            BZ (IZ,IX,IY) = 0D0; BX (IZ,IX,IY) = 0D0
            IF ( (.NOT. G%NEEDSBCW(IZ,IX,IY)) .AND. IX .NE. 1     ) BX (IZ,IX,IY) = BX (IZ,IX,IY) + AW(IZ,IX,IY) * G%YJGN(IGSPEC,IZ,IX-1,IY)
            IF ( (.NOT. G%NEEDSBCE(IZ,IX,IY)) .AND. IX .NE. NCELLX) BX (IZ,IX,IY) = BX (IZ,IX,IY) + AE(IZ,IX,IY) * G%YJGN(IGSPEC,IZ,IX+1,IY)
         ENDDO
         ENDDO
      ENDIF

      IF (SWEEPX .AND. (.NOT. SWEEPZ)) THEN
         DO IX = 1, NCELLX
         DO IZ = 1, NCELLZ
            IF (G%IMASK(IZ,IX,IY)) CYCLE
            BZ (IZ,IX,IY) = 0D0; BX (IZ,IX,IY) = 0D0
            IF ( (.NOT. G%NEEDSBCT(IZ,IX,IY)) .AND. IZ .NE. 1     ) BZ (IZ,IX,IY) = BZ (IZ,IX,IY) + AT(IZ,IX,IY) * G%YJGN(IGSPEC,IZ-1,IX,IY)
            IF ( (.NOT. G%NEEDSBCB(IZ,IX,IY)) .AND. IZ .NE. NCELLZ) BZ (IZ,IX,IY) = BZ (IZ,IX,IY) + AB(IZ,IX,IY) * G%YJGN(IGSPEC,IZ+1,IX,IY)
         ENDDO
         ENDDO
      ENDIF

      BY(:,:,:) = 0.

   ENDIF !End 2D diffusion terms

   IF (NCELLZ .GT. 1 .AND. NCELLX .GT. 1 .AND. NCELLY .GT. 1) THEN !3D

      IF (SWEEPX .AND. SWEEPY .AND. SWEEPZ) THEN
         DO IY = 1, NCELLY
         DO IX = 1, NCELLX
         DO IZ = 1, NCELLZ
            IF (G%IMASK(IZ,IX,IY)) CYCLE
            BZ (IZ,IX,IY) = 0D0; BX (IZ,IX,IY) = 0D0; ; BY (IZ,IX,IY) = 0D0
            IF ( (.NOT. G%NEEDSBCT(IZ,IX,IY)) .AND. IZ .NE. 1     ) BZ (IZ,IX,IY) = BZ (IZ,IX,IY) + AT(IZ,IX,IY) * G%YJGN(IGSPEC,IZ-1,IX,IY)
            IF ( (.NOT. G%NEEDSBCB(IZ,IX,IY)) .AND. IZ .NE. NCELLZ) BZ (IZ,IX,IY) = BZ (IZ,IX,IY) + AB(IZ,IX,IY) * G%YJGN(IGSPEC,IZ+1,IX,IY)
            IF ( (.NOT. G%NEEDSBCW(IZ,IX,IY)) .AND. IX .NE. 1     ) BX (IZ,IX,IY) = BX (IZ,IX,IY) + AW(IZ,IX,IY) * G%YJGN(IGSPEC,IZ,IX-1,IY)
            IF ( (.NOT. G%NEEDSBCE(IZ,IX,IY)) .AND. IX .NE. NCELLX) BX (IZ,IX,IY) = BX (IZ,IX,IY) + AE(IZ,IX,IY) * G%YJGN(IGSPEC,IZ,IX+1,IY)
            IF ( (.NOT. G%NEEDSBCS(IZ,IX,IY)) .AND. IY .NE. 1     ) BY (IZ,IX,IY) = BY (IZ,IX,IY) + AS(IZ,IX,IY) * G%YJGN(IGSPEC,IZ,IX,IY-1)
            IF ( (.NOT. G%NEEDSBCN(IZ,IX,IY)) .AND. IY .NE. NCELLY) BY (IZ,IX,IY) = BY (IZ,IX,IY) + AN(IZ,IX,IY) * G%YJGN(IGSPEC,IZ,IX,IY+1)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (SWEEPX .AND. (.NOT. SWEEPY)  .AND. (.NOT. SWEEPZ)) THEN
         DO IY = 1, NCELLY
         DO IX = 1, NCELLX
         DO IZ = 1, NCELLZ
            IF (G%IMASK(IZ,IX,IY)) CYCLE
            BZ (IZ,IX,IY) = 0D0; BX (IZ,IX,IY) = 0D0; ; BY (IZ,IX,IY) = 0D0
            IF ( (.NOT. G%NEEDSBCT(IZ,IX,IY)) .AND. IZ .NE. 1     ) BZ (IZ,IX,IY) = BZ (IZ,IX,IY) + AT(IZ,IX,IY) * G%YJGN(IGSPEC,IZ-1,IX,IY)
            IF ( (.NOT. G%NEEDSBCB(IZ,IX,IY)) .AND. IZ .NE. NCELLZ) BZ (IZ,IX,IY) = BZ (IZ,IX,IY) + AB(IZ,IX,IY) * G%YJGN(IGSPEC,IZ+1,IX,IY)
            IF ( (.NOT. G%NEEDSBCS(IZ,IX,IY)) .AND. IY .NE. 1     ) BY (IZ,IX,IY) = BY (IZ,IX,IY) + AS(IZ,IX,IY) * G%YJGN(IGSPEC,IZ,IX,IY-1)
            IF ( (.NOT. G%NEEDSBCN(IZ,IX,IY)) .AND. IY .NE. NCELLY) BY (IZ,IX,IY) = BY (IZ,IX,IY) + AN(IZ,IX,IY) * G%YJGN(IGSPEC,IZ,IX,IY+1)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF ((.NOT. SWEEPX) .AND. SWEEPY .AND. (.NOT. SWEEPZ)) THEN
         DO IY = 1, NCELLY
         DO IX = 1, NCELLX
         DO IZ = 1, NCELLZ
            IF (G%IMASK(IZ,IX,IY)) CYCLE
            BZ (IZ,IX,IY) = 0D0; BX (IZ,IX,IY) = 0D0; ; BY (IZ,IX,IY) = 0D0
            IF ( (.NOT. G%NEEDSBCT(IZ,IX,IY)) .AND. IZ .NE. 1     ) BZ (IZ,IX,IY) = BZ (IZ,IX,IY) + AT(IZ,IX,IY) * G%YJGN(IGSPEC,IZ-1,IX,IY)
            IF ( (.NOT. G%NEEDSBCB(IZ,IX,IY)) .AND. IZ .NE. NCELLZ) BZ (IZ,IX,IY) = BZ (IZ,IX,IY) + AB(IZ,IX,IY) * G%YJGN(IGSPEC,IZ+1,IX,IY)
            IF ( (.NOT. G%NEEDSBCW(IZ,IX,IY)) .AND. IX .NE. 1     ) BX (IZ,IX,IY) = BX (IZ,IX,IY) + AW(IZ,IX,IY) * G%YJGN(IGSPEC,IZ,IX-1,IY)
            IF ( (.NOT. G%NEEDSBCE(IZ,IX,IY)) .AND. IX .NE. NCELLX) BX (IZ,IX,IY) = BX (IZ,IX,IY) + AE(IZ,IX,IY) * G%YJGN(IGSPEC,IZ,IX+1,IY)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF ((.NOT. SWEEPX) .AND. (.NOT. SWEEPY) .AND. SWEEPZ) THEN
         DO IY = 1, NCELLY
         DO IX = 1, NCELLX
         DO IZ = 1, NCELLZ
            IF (G%IMASK(IZ,IX,IY)) CYCLE
            BZ (IZ,IX,IY) = 0D0; BX (IZ,IX,IY) = 0D0; ; BY (IZ,IX,IY) = 0D0
            IF ( (.NOT. G%NEEDSBCW(IZ,IX,IY)) .AND. IX .NE. 1     ) BX (IZ,IX,IY) = BX (IZ,IX,IY) + AW(IZ,IX,IY) * G%YJGN(IGSPEC,IZ,IX-1,IY)
            IF ( (.NOT. G%NEEDSBCE(IZ,IX,IY)) .AND. IX .NE. NCELLX) BX (IZ,IX,IY) = BX (IZ,IX,IY) + AE(IZ,IX,IY) * G%YJGN(IGSPEC,IZ,IX+1,IY)
            IF ( (.NOT. G%NEEDSBCS(IZ,IX,IY)) .AND. IY .NE. 1     ) BY (IZ,IX,IY) = BY (IZ,IX,IY) + AS(IZ,IX,IY) * G%YJGN(IGSPEC,IZ,IX,IY-1)
            IF ( (.NOT. G%NEEDSBCN(IZ,IX,IY)) .AND. IY .NE. NCELLY) BY (IZ,IX,IY) = BY (IZ,IX,IY) + AN(IZ,IX,IY) * G%YJGN(IGSPEC,IZ,IX,IY+1)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

   ENDIF !End 3D diffusion terms



   IF (SWEEPX) THEN
      !$omp PARALLEL DO SCHEDULE(DYNAMIC,4) PRIVATE(IZ,IX,IY)
      DO IZ = 1, NCELLZ
      DO IY = 1, NCELLY
         DO IX = 1, NCELLX
            IF (G%IMASK(IZ,IX,IY)) THEN
               B(IZ,IX,IY) = G%YJG(IGSPEC,IZ,IX,IY)
               CYCLE
            ENDIF

            APTOT(IZ,IX,IY) = 0D0
            BTOT (IZ,IX,IY) = 0D0

            DO I = 1, 4
               IF (I .EQ. 1) SOURCE = BZ           (       IZ,IX,IY)
               IF (I .EQ. 2) SOURCE = BY           (       IZ,IX,IY)
               IF (I .EQ. 3) SOURCE = BBC          (       IZ,IX,IY)
               IF (I .EQ. 4) SOURCE = GOMEGA3DXDYDZ(IGSPEC,IZ,IX,IY)

               IF (SOURCE .LT. 0D0) THEN
                  APTOT(IZ,IX,IY) = APTOT(IZ,IX,IY) - SOURCE / MAX(G%YJGN(IGSPEC,IZ,IX,IY),EPS)
               ELSE
                  BTOT(IZ,IX,IY) = BTOT(IZ,IX,IY) + SOURCE
               ENDIF
            ENDDO

            AP(IZ,IX,IY) = (AP1(IZ,IX,IY) + APTOT(IZ,IX,IY)) * RALPHA_YJG
            B (IZ,IX,IY) = AP0YJG(IGSPEC,IZ,IX,IY) + BTOT(IZ,IX,IY) + AP(IZ,IX,IY)*YJGNOLD(IGSPEC,IZ,IX,IY) * OMALPHA_YJG
         ENDDO
         CALL TDMA_SOLVER_GENERAL(NCELLX,PTR(IZ,:,IY),AP(IZ,:,IY),AE(IZ,:,IY),AW(IZ,:,IY),B(IZ,:,IY))
      ENDDO
      ENDDO
      !$omp END PARALLEL DO
   ENDIF 

   IF (SWEEPY) THEN

      !$omp PARALLEL DO SCHEDULE(DYNAMIC,2) PRIVATE(IZ,IX,IY)
      DO IZ = 1, NCELLZ
      DO IX = 1, NCELLX
         DO IY = 1, NCELLY
            IF (G%IMASK(IZ,IX,IY)) THEN
               B(IZ,IX,IY) = G%YJG(IGSPEC,IZ,IX,IY)
               CYCLE
            ENDIF

            APTOT(IZ,IX,IY) = 0D0
            BTOT (IZ,IX,IY) = 0D0

            DO I = 1, 4
               IF (I .EQ. 1) SOURCE = BZ           (       IZ,IX,IY)
               IF (I .EQ. 2) SOURCE = BX           (       IZ,IX,IY)
               IF (I .EQ. 3) SOURCE = BBC          (       IZ,IX,IY)
               IF (I .EQ. 4) SOURCE = GOMEGA3DXDYDZ(IGSPEC,IZ,IX,IY)

               IF (SOURCE .LT. 0D0) THEN
                  APTOT(IZ,IX,IY) = APTOT(IZ,IX,IY) - SOURCE / MAX(G%YJGN(IGSPEC,IZ,IX,IY),EPS)
               ELSE
                  BTOT(IZ,IX,IY) = BTOT(IZ,IX,IY) + SOURCE
               ENDIF
            ENDDO

            AP(IZ,IX,IY) = (AP1(IZ,IX,IY) + APTOT(IZ,IX,IY)) * RALPHA_YJG
            B (IZ,IX,IY) = AP0YJG(IGSPEC,IZ,IX,IY) + BTOT(IZ,IX,IY) + AP(IZ,IX,IY)*YJGNOLD(IGSPEC,IZ,IX,IY) * OMALPHA_YJG
         ENDDO
         CALL TDMA_SOLVER_GENERAL(NCELLY,PTR(IZ,IX,:),AP(IZ,IX,:),AN(IZ,IX,:),AS(IZ,IX,:),B(IZ,IX,:))
      ENDDO
      ENDDO
      !$omp END PARALLEL DO
   ENDIF

   IF(SWEEPZ) THEN !SWEEPZ
      DO IX = 1, NCELLX
      DO IY = 1, NCELLY
         DO IZ = 1, NCELLZ

            IF (G%IMASK(IZ,IX,IY)) THEN
               B(IZ,IX,IY) = G%YJG(IGSPEC,IZ,IX,IY)
               CYCLE
            ENDIF

            APTOT(IZ,IX,IY) = 0D0
            BTOT (IZ,IX,IY) = 0D0

            DO I = 1, 4
               IF (I .EQ. 1) SOURCE = BX           (       IZ,IX,IY)
               IF (I .EQ. 2) SOURCE = BY           (       IZ,IX,IY)
               IF (I .EQ. 3) SOURCE = BBC          (       IZ,IX,IY)
               IF (I .EQ. 4) SOURCE = GOMEGA3DXDYDZ(IGSPEC,IZ,IX,IY)

               IF (SOURCE .LT. 0D0) THEN
                  APTOT(IZ,IX,IY) = APTOT(IZ,IX,IY) - SOURCE / MAX(G%YJGN(IGSPEC,IZ,IX,IY),EPS)
               ELSE
                  BTOT(IZ,IX,IY) = BTOT(IZ,IX,IY) + SOURCE
               ENDIF
            ENDDO
                  
            AP(IZ,IX,IY) = (AP1(IZ,IX,IY) + APTOT(IZ,IX,IY)) * RALPHA_YJG
            B (IZ,IX,IY) = AP0YJG(IGSPEC,IZ,IX,IY) + BTOT(IZ,IX,IY) + AP(IZ,IX,IY)*YJGNOLD(IGSPEC,IZ,IX,IY) * OMALPHA_YJG
         ENDDO
         CALL TDMA_SOLVER_GENERAL(NCELLZ,PTR(:,IX,IY),AP(:,IX,IY),AB(:,IX,IY),AT(:,IX,IY),B(:,IX,IY))
      ENDDO
      ENDDO
   ENDIF


ENDDO !IGSPEC ENDDO


! Check convergence (relative tolerance)

! This is how convergence was checked in 0.700:
!EPS = 1D-5
!YJGDIFF(:,:) = ABS(YJGN(:,:)-YJGNOLD(:,:)) 
!WHERE(YJGNOLD(:,:) .LE. EPS) YJGDIFF(:,:) = YJGDIFF(:,:) / EPS
!WHERE(YJGNOLD(:,:) .GT. EPS) YJGDIFF(:,:) = YJGDIFF(:,:) / YJGNOLD(:,:)
!MAXDIFF = MAXVAL(YJGDIFF(:,:))

IF (G%CONV_INFO%ITER .EQ. 1) THEN
   G%CONV_INFO%CONVERGED_YJG (:) = .FALSE.
   G%RESIDUAL_YJG(:,:,:,:) = 1.
ELSE
   ALL_CELLS_CONVERGED       (:) = .TRUE.
   G%CONV_INFO%CONVERGED_YJG (:) = .TRUE.
   !$omp PARALLEL DO SCHEDULE(STATIC) PRIVATE(RESIDABS)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
      IF (G%IMASK(IZ,IX,IY)) THEN
         G%RESIDUAL_YJG(IGSPEC,IZ,IX,IY)= 0D0
         CYCLE
      ENDIF
      YJSUM = 0D0
      DO IGSPEC = 1, GPROP%NGSPEC
         YJSUM = YJSUM + G%YJGN(IGSPEC,IZ,IX,IY)
         
         RESIDABS = ABS(G%YJGN(IGSPEC,IZ,IX,IY)-YJGNOLD(IGSPEC,IZ,IX,IY))
         G%RESIDUAL_YJG(IGSPEC,IZ,IX,IY) =  RESIDABS / MAX(YJGNOLD(IGSPEC,IZ,IX,IY),GPG%EPS_YJG)
         IF (G%YJGN(IGSPEC,IZ,IX,IY).NE.G%YJGN(IGSPEC,IZ,IX,IY) .OR. G%YJGN(IGSPEC,IZ,IX,IY).EQ.GPG%POSINF .OR. G%YJGN(IGSPEC,IZ,IX,IY).EQ.GPG%NEGINF) THEN
            G%RESIDUAL_YJG(IGSPEC,IZ,IX,IY)=9D9
            GPG%NAN = .TRUE.
         ENDIF
         IF (G%RESIDUAL_YJG(IGSPEC,IZ,IX,IY) .GT. GPG%YJTOL*ALPHA_YJG) THEN
            ALL_CELLS_CONVERGED(IGSPEC) = .FALSE.
            ALL_CELLS_CONVERGED(0     ) = .FALSE.
            
         ENDIF
         IF (YJSUM .LT. 0.999 .OR. YJSUM .GT. 1.001) THEN
            G%RESIDUAL_YJG(:,IZ,IX,IY) = 10.
            ALL_CELLS_CONVERGED(IGSPEC) = .FALSE.
            ALL_CELLS_CONVERGED(0     ) = .FALSE.
         ENDIF
   ENDDO

   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   DO IGSPEC = 0, GPROP%NGSPEC
      IF (.NOT. ALL_CELLS_CONVERGED(IGSPEC)) THEN
         G%CONV_INFO%CONVERGED_YJG(IGSPEC)= .FALSE.
         G%CONV_INFO%CONVERGED_ALL = .FALSE.
      ENDIF
      IF (G%CONV_INFO%CONVERGED_YJG(IGSPEC) .AND. (G%CONV_INFO%ITER_YJG(IGSPEC) .EQ. 0)) THEN
         G%CONV_INFO%ITER_YJG(IGSPEC) = G%CONV_INFO%ITER
      ENDIF
   ENDDO
ENDIF

! *****************************************************************************
END SUBROUTINE CONVECTIVE_DIFFUSIVE_SOLVER
! *****************************************************************************


! *****************************************************************************
SUBROUTINE CALCULATE_GAS_DENSITY(IMESH)
! *****************************************************************************

INTEGER, INTENT(IN) :: IMESH
INTEGER :: IZ,IX,IY,IGSPEC
REAL(EB), POINTER, DIMENSION (:,:,:) :: MWN
INTEGER :: NCELLX,NCELLY,NCELLZ,NGSPEC


G=>GPM(IMESH)

MWN => G%RWORK01


NCELLX = G%NCELLX
NCELLY = G%NCELLY
NCELLZ = G%NCELLZ
NGSPEC = GPROP%NGSPEC

! Calculate new gas density
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ,IGSPEC) &
!$OMP SHARED(MWN,G,GPROP,GPG,NCELLX,NCELLY,NCELLZ,NGSPEC) COLLAPSE(3)
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ
   MWN(IZ,IX,IY) = 0D0
   DO IGSPEC = 1, NGSPEC
      MWN(IZ,IX,IY) = MWN(IZ,IX,IY) + G%YJGN(IGSPEC,IZ,IX,IY) / GPROP%M(IGSPEC)
   ENDDO
   MWN(IZ,IX,IY) = 1D0 / MWN(IZ,IX,IY)
   IF (GPG%SOLVE_PRESSURE) THEN
      G%RGN(IZ,IX,IY) = RHOGOFT(G%PN(IZ,IX,IY), MWN(IZ,IX,IY), G%TP(IZ,IX,IY))
   ELSE
      G%RGN(IZ,IX,IY) = RHOGOFT(GPG%P0, MWN(IZ,IX,IY), G%TP(IZ,IX,IY))
   ENDIF
ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO


! *****************************************************************************
END SUBROUTINE CALCULATE_GAS_DENSITY
! *****************************************************************************



! *****************************************************************************
SUBROUTINE GAS_ENERGY_SOURCE_TERMS(IMESH)
! *****************************************************************************
INTEGER, INTENT(IN) :: IMESH
INTEGER :: IGSPEC,IZ,IX,IY,IRXN
INTEGER :: NCELLZ, NCELLX, NCELLY

G => GPM(IMESH)

NCELLZ = G%NCELLZ
NCELLX = G%NCELLX
NCELLY = G%NCELLY
! Gas-phase energy source terms:
IF (GPG%SOLVE_GAS_ENERGY .AND. (.NOT. GPG%THERMAL_EQUILIBRIUM)) THEN
   IF (GPG%HCV .LT. 0D0) THEN
      !$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(3) DEFAULT(SHARED) PRIVATE(IX,IY,IZ)
      DO IY = 1, NCELLY
      DO IX = 1, NCELLX
      DO IZ = 1, NCELLZ
         G%SHGP(IZ,IX,IY) = G%HCV(IZ,IX,IY)*G%TPN(IZ,IX,IY)
         G%SHGM(IZ,IX,IY) = G%HCV(IZ,IX,IY)*G%TGN(IZ,IX,IY)
      ENDDO
      ENDDO
      ENDDO
   !$OMP END PARALLEL DO
   ELSE
      !$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(3) DEFAULT(SHARED) PRIVATE(IX,IY,IZ)
      DO IY = 1, NCELLY
      DO IX = 1, NCELLX
      DO IZ = 1, NCELLZ
         G%SHGP(IZ,IX,IY) = GPG%HCV*G%TPN(IZ,IX,IY)
         G%SHGM(IZ,IX,IY) = GPG%HCV*G%TGN(IZ,IX,IY)
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   ENDIF

   !These terms cancel with -omegasfg * hg only if gases are produced at Tgas (not Tsolid)
   IF (GPG%GASES_PRODUCED_AT_TSOLID) THEN
      !$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(3) DEFAULT(SHARED) PRIVATE(IX,IY,IZ,IGSPEC)
      DO IY = 1, NCELLY
      DO IX = 1, NCELLX
      DO IZ = 1, NCELLZ
      IF (.NOT. G%IS_REACTING(0,IZ,IX,IY)) CYCLE
      DO IGSPEC = 1, GPROP%NGSPEC
         G%SHGP(IZ,IX,IY) = G%SHGP(IZ,IX,IY) + GPROP%CPG * (G%TPN(IZ,IX,IY) - GPG%TDATUM) * G%GOMEGA (1,IGSPEC,IZ,IX,IY) ! 1
         G%SHGM(IZ,IX,IY) = G%SHGM(IZ,IX,IY) + GPROP%CPG * (G%TPN(IZ,IX,IY) - GPG%TDATUM) * G%GOMEGA (2,IGSPEC,IZ,IX,IY) ! 3
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   ENDIF

    
   DO IRXN = 1, GPROP%NHGRXN
      IF (HGRXN(IRXN)%DH .LT. 0D0) THEN !Exothermic
         !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IRXN,IX,IY,IZ) DEFAULT(SHARED) COLLAPSE(3) 
         DO IY = 1, NCELLY
         DO IX = 1, NCELLX
         DO IZ = 1, NCELLZ        
            G%SHGP(IZ,IX,IY) = G%SHGP(IZ,IX,IY) - G%HGRR(IRXN,IZ,IX,IY) * HGRXN(IRXN)%DH
         ENDDO
         ENDDO
         ENDDO
      !$OMP END PARALLEL DO
      ELSE
         !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IRXN,IX,IY,IZ) DEFAULT(SHARED) COLLAPSE(3) 
         DO IY = 1, NCELLY
         DO IX = 1, NCELLX
         DO IZ = 1, NCELLZ
            G%SHGM(IZ,IX,IY) = G%SHGM(IZ,IX,IY) + G%HGRR(IRXN,IZ,IX,IY) * HGRXN(IRXN)%DH
         ENDDO
         ENDDO
         ENDDO
         !$OMP END PARALLEL DO
      ENDIF         
   ENDDO
ENDIF


! *****************************************************************************
END SUBROUTINE GAS_ENERGY_SOURCE_TERMS
! *****************************************************************************


!This generalized routine calculates interface quantities:
! *****************************************************************************
SUBROUTINE CALCULATE_GAS_INTERFACE_QUANTITIES(IMESH,CQUANTITY)
! *****************************************************************************

INTEGER, INTENT(IN) :: IMESH
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
 
   CASE('PERMONU')

      IF (NCELLZ .GT. 1 .AND. NCELLX .LE. 1) THEN
         IX = 1; IY = 1
         !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IZ) SHARED(G,GAMPZ,NCELLZ)
         DO IZ = 1, NCELLZ
            IF (G%IMASK(IZ,IX,IY)) THEN; GAMPZ(IZ,IX,IY)=0D0; CYCLE; ENDIF
            GAMPZ(IZ,IX,IY) = G%PERMZ(IZ,IX,IY) / G%D12(IZ,IX,IY) !Gamma-z (point)
            IF(G%IMASK(IZ,IX,IY)) GAMPZ(IZ,IX,IY) = GAMPZ(IZ,IX,IY) * 1D6
         ENDDO
         !$OMP END PARALLEL DO
         GAMT => G%PERMONUT
         GAMB => G%PERMONUB
      ENDIF

      IF (NCELLZ .GT. 1  .AND. NCELLX .GT. 1 .AND. NCELLY .LE. 1) THEN
         IY = 1
         !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IZ) SHARED(G,GAMPZ,GAMPX,GAMPY,NCELLX,NCELLZ) COLLAPSE(2)
         DO IX = 1, NCELLX
         DO IZ = 1, NCELLZ
            GAMPZ(IZ,IX,IY) = G%PERMZ(IZ,IX,IY) / G%D12(IZ,IX,IY) !Gamma-z (point)
            GAMPX(IZ,IX,IY) = G%PERMX(IZ,IX,IY) / G%D12(IZ,IX,IY) !Gamma-x (point)
            IF(G%IMASK(IZ,IX,IY)) THEN 
               GAMPZ(IZ,IX,IY) = GAMPZ(IZ,IX,IY) * 1D6
               GAMPX(IZ,IX,IY) = GAMPX(IZ,IX,IY) * 1D6
            ENDIF
         ENDDO
         ENDDO
         !$OMP END PARALLEL DO
         GAMT => G%PERMONUT
         GAMB => G%PERMONUB
         GAME => G%PERMONUE
         GAMW => G%PERMONUW
      ENDIF

      IF (NCELLZ .GT. 1 .AND. NCELLX .GT. 1 .AND. NCELLY .GT. 1) THEN
         !$omp PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) SHARED(G,GAMPZ,GAMPX,GAMPY,NCELLX,NCELLY,NCELLZ) COLLAPSE(3)
         DO IY = 1, NCELLY
         DO IX = 1, NCELLX
         DO IZ = 1, NCELLZ
            GAMPZ(IZ,IX,IY) = G%PERMZ(IZ,IX,IY) / G%D12(IZ,IX,IY) !Gamma-z (point)
            GAMPX(IZ,IX,IY) = G%PERMX(IZ,IX,IY) / G%D12(IZ,IX,IY) !Gamma-x (point)
            GAMPY(IZ,IX,IY) = G%PERMY(IZ,IX,IY) / G%D12(IZ,IX,IY) !Gamma-y (point)
            IF(G%IMASK(IZ,IX,IY)) THEN 
               GAMPZ(IZ,IX,IY) = GAMPZ(IZ,IX,IY) * 1D6
               GAMPX(IZ,IX,IY) = GAMPX(IZ,IX,IY) * 1D6
               GAMPY(IZ,IX,IY) = GAMPY(IZ,IX,IY) * 1D6
            ENDIF
         ENDDO
         ENDDO
         ENDDO
         !$omp END PARALLEL DO
         GAMT => G%PERMONUT
         GAMB => G%PERMONUB
         GAME => G%PERMONUE
         GAMW => G%PERMONUW
         GAMN => G%PERMONUN
         GAMS => G%PERMONUS
      ENDIF

   CASE('PSIRGD')
      IF (NCELLZ .GT. 1) THEN
         !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) SHARED(G,GPG,GAMPZ,NCELLX,NCELLY,NCELLZ) COLLAPSE(3)
         DO IY = 1, NCELLY
         DO IX = 1, NCELLX
         DO IZ = 1, NCELLZ
            GAMPZ(IZ,IX,IY) = GPG%TORTUOSITY_FACTOR * G%POROSS(IZ,IX,IY) * G%RGN(IZ,IX,IY) * G%D12(IZ,IX,IY) !Gamma-z (point)
         ENDDO
         ENDDO
         ENDDO
         !$OMP END PARALLEL DO
         GAMT => G%PSIRGDT
         GAMB => G%PSIRGDB
      ENDIF

      IF (NCELLX .GT. 1) THEN
         GAMPX = GAMPZ
         GAME => G%PSIRGDE
         GAMW => G%PSIRGDW
      ENDIF

      IF (NCELLY .GT. 1) THEN
         GAMPY = GAMPZ
         GAMN => G%PSIRGDN
         GAMS => G%PSIRGDS
      ENDIF

   CASE('RG')
      IF (NCELLZ .GT. 1) THEN
         GAMPZ = G%RGN !Gamma-z (point)
         GAMT => G%RGNT
         GAMB => G%RGNB
      ENDIF

      IF (NCELLX .GT. 1) THEN   
         GAMPX = GAMPZ !G%RGN !Gamma-x (point)
         GAME => G%RGNE
         GAMW => G%RGNW
      ENDIF

      IF (NCELLY .GT. 1) THEN   
         GAMPY = GAMPZ !G%RGN !Gamma-y (point)
         GAMN => G%RGNN
         GAMS => G%RGNS
      ENDIF

   CASE DEFAULT

      CONTINUE
          
END SELECT


CALL CALCULATE_INTERFACE_QUANTITIES(IMESH,GAMPZ,GAMPX,GAMPY, &
                                    GAMT,GAMB,GAME,GAMW,GAMN,GAMS)

! *****************************************************************************
END SUBROUTINE CALCULATE_GAS_INTERFACE_QUANTITIES
! *****************************************************************************


! *****************************************************************************
SUBROUTINE CALCULATE_WEIGHTED_POINT_PERMEABILITY(IMESH)
! *****************************************************************************
INTEGER, INTENT(IN) :: IMESH
INTEGER :: IZ,IX,IY,ISPEC
REAL(EB) ::  NUMER, DENOM, OMPOROS, FPOROS
INTEGER ::  NCELLX, NCELLY, NCELLZ,NSSPEC

G => GPM(IMESH)

NCELLX = G%NCELLX
NCELLY = G%NCELLY
NCELLZ = G%NCELLZ
NSSPEC = SPROP%NSSPEC

IF (GPG%KOZENY_CARMAN) THEN
   !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(PRIVATE) SHARED(G,SPROP,NCELLX,NCELLY,NCELLZ,NSSPEC) COLLAPSE(3)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ

      G%PERMZ(IZ,IX,IY) = 0D0
      G%PERMX(IZ,IX,IY) = 0D0
      G%PERMY(IZ,IX,IY) = 0D0

      IF (G%IMASK(IZ,IX,IY)) CYCLE

      NUMER   = G%POROSSN(IZ,IX,IY) * G%POROSSN(IZ,IX,IY) * G%POROSSN(IZ,IX,IY)
      OMPOROS = 1D0 - G%POROSSN(IZ,IX,IY)
      DENOM   = MAX(OMPOROS*OMPOROS, 1D-20)
      FPOROS  = NUMER / DENOM

      DO ISPEC = 1, NSSPEC
         IF (G%XIN(ISPEC,IZ,IX,IY) .LT. EPSILON_FB) CYCLE
         G%PERMZ(IZ,IX,IY) = G%PERMZ(IZ,IX,IY) + G%XIN(ISPEC,IZ,IX,IY) * SPROP%PERMZ(ISPEC) * FPOROS
         G%PERMX(IZ,IX,IY) = G%PERMX(IZ,IX,IY) + G%XIN(ISPEC,IZ,IX,IY) * SPROP%PERMX(ISPEC) * FPOROS
         G%PERMY(IZ,IX,IY) = G%PERMY(IZ,IX,IY) + G%XIN(ISPEC,IZ,IX,IY) * SPROP%PERMY(ISPEC) * FPOROS
      ENDDO

   ENDDO
   ENDDO
   ENDDO
   !$omp END PARALLEL DO


ELSE
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ,ISPEC) SHARED(G,SPROP,NCELLX,NCELLY,NCELLZ,NSSPEC) COLLAPSE(3)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ

      G%PERMZ(IZ,IX,IY) = 0D0
      G%PERMX(IZ,IX,IY) = 0D0
      G%PERMY(IZ,IX,IY) = 0D0

      IF (G%IMASK(IZ,IX,IY)) CYCLE

      DO ISPEC = 1, NSSPEC 
         IF (G%XIN(ISPEC,IZ,IX,IY) .LT. EPSILON_FB) CYCLE  
         G%PERMZ(IZ,IX,IY) = G%PERMZ(IZ,IX,IY) + G%XIN(ISPEC,IZ,IX,IY)*SPROP%PERMZ(ISPEC)   
         G%PERMX(IZ,IX,IY) = G%PERMX(IZ,IX,IY) + G%XIN(ISPEC,IZ,IX,IY)*SPROP%PERMX(ISPEC)   
         G%PERMY(IZ,IX,IY) = G%PERMY(IZ,IX,IY) + G%XIN(ISPEC,IZ,IX,IY)*SPROP%PERMY(ISPEC)   
      ENDDO

   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDIF



! *****************************************************************************
END SUBROUTINE CALCULATE_WEIGHTED_POINT_PERMEABILITY
! *****************************************************************************

! *****************************************************************************
SUBROUTINE GET_D(IMESH,IBG,IO2)
! *****************************************************************************

INTEGER, INTENT(IN) :: IMESH,IBG,IO2 !Background species, O2 index
REAL(EB), POINTER, DIMENSION(:,:,:) :: TSTAR,OMEGAD,NUMER,PRES
REAL(EB), SAVE :: EPSOK = -1D0, SIG2 = -1D0, NUMER1 = -1D0 
INTEGER :: IZ,IX,IY, NCELLX, NCELLY, NCELLZ


G=>GPM(IMESH)

NCELLX = G%NCELLX
NCELLY = G%NCELLY
NCELLZ = G%NCELLZ

TSTAR  => G%RWORK01
OMEGAD => G%RWORK02
NUMER  => G%RWORK03
PRES   => G%RWORK04

! It's assumed all species have the effective binary diffusion coefficient
! of oxygen into the "background species", usually fuel (pyrolysate)
IF (EPSOK .LT. 0D0) THEN
   EPSOK  = SQRT(GPROP%EPSOK(IO2)*GPROP%EPSOK(IBG))
   SIG2   = (0.5D0*(GPROP%SIGMA(IO2)+GPROP%SIGMA(IBG)))**2D0
   NUMER1 = 1D0/GPROP%M(IO2) + 1D0/GPROP%M(IBG)
ENDIF

IF (GPG%SOLVE_PRESSURE) THEN
   PRES(:,:,:) = G%PN(:,:,:)
ELSE
   PRES(:,:,:) = GPG%P0
ENDIF

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) DEFAULT(SHARED) COLLAPSE(3)
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ
   IF (G%IMASK(IZ,IX,IY)) THEN
      G%D12(IZ,IX,IY) = 1D-20
   ELSE
      TSTAR (IZ,IX,IY) = G%TPN(IZ,IX,IY) / EPSOK
      OMEGAD(IZ,IX,IY) = AD(1)*TSTAR(IZ,IX,IY)**AD(2) + AD(3)*EXP(TSTAR(IZ,IX,IY)*AD(4)) + AD(5)*EXP(AD(6)*TSTAR(IZ,IX,IY))
      NUMER (IZ,IX,IY) = SQRT(NUMER1 * G%TPN(IZ,IX,IY)*G%TPN(IZ,IX,IY)*G%TPN(IZ,IX,IY))
      G%D12 (IZ,IX,IY) = 0.0188129*NUMER(IZ,IX,IY) / (PRES(IZ,IX,IY)*SIG2*OMEGAD(IZ,IX,IY))
   ENDIF
ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO

! *****************************************************************************
END SUBROUTINE GET_D
! *****************************************************************************

! *****************************************************************************
SUBROUTINE CALC_HCV(IMESH,NCELLZ,NCELLX,NCELLY)
! *****************************************************************************
INTEGER, INTENT(IN) :: IMESH,NCELLZ,NCELLX,NCELLY
INTEGER :: IX,IY,ISPEC

REAL(EB), DIMENSION (1:NCELLZ) :: MDOTPP,DP


G=>GPM(IMESH)

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,ISPEC,DP,MDOTPP) SHARED(G,SPROP,GPG,NCELLZ) COLLAPSE(2)
DO IY = 1, NCELLY
DO IX = 1, NCELLX

   ! Calculate average pore diameter:
   DP(:) = 0D0
   DO ISPEC = 1, SPROP%NSSPEC
      DP(:) = DP(:) + G%YI(ISPEC,:,IX,IY) * SPROP%PORE_DIAMETER(ISPEC)
   ENDDO
   WHERE (DP(:) .LT. 1D-10) DP(:) = 1D-10

   ! Get mass flux:
   IF (GPG%SOLVE_PRESSURE) THEN
      MDOTPP(:) =  G%MDOTPPDARCYT(:,IX,IY)
   ELSE
      MDOTPP(:) = -G%MDOTPPZ(0,:,IX,IY)
   ENDIF

   IF (NCELLX .GT. 1) THEN
      MDOTPP(:) = SQRT(MDOTPP(:)**2D0 + G%MDOTPPDARCYE(:,IX,IY)**2D0)
   ELSE
      MDOTPP(:) = ABS(MDOTPP(:)) 
   ENDIF
      
   !Calculate Reynolds number:      
   G%RE(:,IX,IY) = MDOTPP(:) * DP(:) / (G%RGN(:,IX,IY)*G%D12(:,IX,IY)*MAX(G%POROSSN(:,IX,IY), 1D-10) )
      
   G%NU(:,IX,IY) = GPG%NU_A + GPG%NU_B * G%RE(:,IX,IY)**GPG%NU_C
      
   G%HCV(:,IX,IY) = G%NU(:,IX,IY) * G%RGN(:,IX,IY)*G%D12(:,IX,IY)*GPROP%CPG / (DP(:)**2D0)

ENDDO
ENDDO
!$OMP END PARALLEL DO


! *****************************************************************************
END SUBROUTINE CALC_HCV
! *****************************************************************************



! *****************************************************************************
END MODULE GPYRO_GAS
! *****************************************************************************
