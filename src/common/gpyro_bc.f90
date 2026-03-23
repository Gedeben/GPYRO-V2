! *****************************************************************************
MODULE GPYRO_BC
! *****************************************************************************

USE PREC
USE GPYRO_VARS
USE GPYRO_FUNCS

IMPLICIT NONE

CONTAINS

! *****************************************************************************
SUBROUTINE GET_ALL_BOUNDARY_CONDITIONS(IMESH,T)
! *****************************************************************************

USE GPYRO_VARS

INTEGER, INTENT(IN) :: IMESH
REAL(EB), INTENT(IN) :: T
INTEGER :: IL, IH, IZ, IX, IY,IOR_GPYRO, ICOUNT, IMESH_FDS
REAL(EB) :: F

G=>GPM(IMESH)


!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICOUNT, IMESH_FDS, IX, IY, IZ, IOR_GPYRO, IL, IH, F) &
!$OMP SHARED(NGPYRO_FACES_NEEDING_BCS, GP_BOUDARYS, IGPYRO_TYPE, IMESH, G )
DO ICOUNT = 1, NGPYRO_FACES_NEEDING_BCS(IMESH)
   IMESH_FDS = GP_BOUDARYS(IMESH)%IMESH_FDS(ICOUNT)
   IF ((IGPYRO_TYPE .EQ. 2) .AND. IMESH_FDS.NE. 0) CYCLE  !  If this is à FDS simulation, and this cell have FDS BC skip 
   IF (GP_BOUDARYS(IMESH)%IMESH_GPYRO(ICOUNT) .NE. IMESH) CYCLE
   IZ        = GP_BOUDARYS(IMESH)%IZ_GPYRO (ICOUNT)
   IX        = GP_BOUDARYS(IMESH)%IX_GPYRO (ICOUNT)
   IY        = GP_BOUDARYS(IMESH)%IY_GPYRO (ICOUNT)
   IOR_GPYRO = GP_BOUDARYS(IMESH)%IOR_GPYRO(ICOUNT)
   SELECT CASE(IOR_GPYRO)
   CASE ( 3); CALL GET_BC_F(G%SURF_IDX_BCT(IZ,IX,IY),T,IL,IH,F) !-z boundary condition: Top (it is upside-down)
   CASE (-3); CALL GET_BC_F(G%SURF_IDX_BCB(IZ,IX,IY),T,IL,IH,F) !+z boundary condition: Bottom
   CASE (-1); CALL GET_BC_F(G%SURF_IDX_BCW(IZ,IX,IY),T,IL,IH,F) !-x boundary condition
   CASE ( 1); CALL GET_BC_F(G%SURF_IDX_BCE(IZ,IX,IY),T,IL,IH,F) !+x boundary condition
   CASE (-2); CALL GET_BC_F(G%SURF_IDX_BCS(IZ,IX,IY),T,IL,IH,F) !-y boundary condition
   CASE ( 2); CALL GET_BC_F(G%SURF_IDX_BCN(IZ,IX,IY),T,IL,IH,F) !+y boundary condition
   END SELECT
   CALL INTERPOLATE_BCS(IL,IH,F,IZ,IX,IY, IOR_GPYRO)
ENDDO
!$omp END PARALLEL DO

GPBCP(1:,1:,1:,-3:)=>G%GPYRO_BOUNDARY_CONDITION(1:,1:,1:,-3:)
 

CONTAINS

! *****************************************************************************
SUBROUTINE GET_BC_F(SURF_IDX_TARG,T,IL,IH,F)
! *****************************************************************************

INTEGER, INTENT(IN) :: SURF_IDX_TARG
REAL(EB), INTENT(IN) :: T
INTEGER, INTENT(OUT) :: IL, IH
REAL(EB), INTENT(OUT) :: F

INTEGER :: I 

IL = -1
IH = -1
DO I = 1, GPG%NSURF_IDX
   IF (GPG%ALLBC(I)%SURF_IDX .EQ. SURF_IDX_TARG) THEN

      IF (GPG%ALLBC(I)%T .LE. T) THEN 
         IL = I
         IH = IL
         IF (I .LT. GPG%NSURF_IDX) THEN
            IF (GPG%ALLBC(I+1)%SURF_IDX .EQ. SURF_IDX_TARG) IH = IL + 1
         ENDIF
      ENDIF
      
   ENDIF
ENDDO

IF (IH .EQ. IL) THEN
   F  = 0.
ELSE
   F = (T - GPG%ALLBC(IL)%T) / (GPG%ALLBC(IH)%T - GPG%ALLBC(IL)%T)
ENDIF

! *****************************************************************************
END SUBROUTINE GET_BC_F
! *****************************************************************************

! *****************************************************************************
SUBROUTINE INTERPOLATE_BCS(IL,IH,F,IZ,IX,IY,IOR)
! *****************************************************************************

INTEGER, INTENT(IN) :: IL,IH,IZ,IX,IY,IOR
REAL(EB), INTENT(IN) :: F
INTEGER :: ISPEC
REAL(EB) :: HFIXEDLO, HFIXEDHI

! Linearly interpolate real quantities:
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%QE      = GPG%ALLBC(IL)%QE      + F * (GPG%ALLBC(IH)%QE      - GPG%ALLBC(IL)%QE      )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%HC0     = GPG%ALLBC(IL)%HC      + F * (GPG%ALLBC(IH)%HC      - GPG%ALLBC(IL)%HC      )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%NHC     = GPG%ALLBC(IL)%NHC     + F * (GPG%ALLBC(IH)%NHC     - GPG%ALLBC(IL)%NHC     )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%TINF    = GPG%ALLBC(IL)%TINF    + F * (GPG%ALLBC(IH)%TINF    - GPG%ALLBC(IL)%TINF    )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%TFIXED  = GPG%ALLBC(IL)%TFIXED  + F * (GPG%ALLBC(IH)%TFIXED  - GPG%ALLBC(IL)%TFIXED  )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%PRES    = GPG%ALLBC(IL)%PRES    + F * (GPG%ALLBC(IH)%PRES    - GPG%ALLBC(IL)%PRES    )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%MFLUX   = GPG%ALLBC(IL)%MDOTPP  + F * (GPG%ALLBC(IH)%MDOTPP  - GPG%ALLBC(IL)%MDOTPP  )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%HM0     = GPG%ALLBC(IL)%HM      + F * (GPG%ALLBC(IH)%HM      - GPG%ALLBC(IL)%HM      )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%QEG     = GPG%ALLBC(IL)%QEG     + F * (GPG%ALLBC(IH)%QEG     - GPG%ALLBC(IL)%QEG     )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%HC0G    = GPG%ALLBC(IL)%HCG     + F * (GPG%ALLBC(IH)%HCG     - GPG%ALLBC(IL)%HCG     )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%TINFG   = GPG%ALLBC(IL)%TINFG   + F * (GPG%ALLBC(IH)%TINFG   - GPG%ALLBC(IL)%TINFG   )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%TFIXEDG = GPG%ALLBC(IL)%TFIXEDG + F * (GPG%ALLBC(IH)%TFIXEDG - GPG%ALLBC(IL)%TFIXEDG )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%PRES    = GPG%ALLBC(IL)%PRES    + F * (GPG%ALLBC(IH)%PRES    - GPG%ALLBC(IL)%PRES    )
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%YJINF(1:GPROP%NGSPEC) = GPG%ALLBC(IL)%YJINF(1:GPROP%NGSPEC) + F * (GPG%ALLBC(IH)%YJINF(1:GPROP%NGSPEC) - GPG%ALLBC(IL)%YJINF(1:GPROP%NGSPEC))

! Reradiation is discrete so use value at IL:
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%RERAD = GPG%ALLBC(IL)%RERADIATION

!Maybe we should rethink all this subroutine
! Set solid enthalpy:
HFIXEDLO = 0.
HFIXEDHI = 0.
DO ISPEC = 1, SPROP%NSSPEC
   HFIXEDLO = HFIXEDLO + G%YI(ISPEC,IZ,IX,IY) * HOFT(ISPEC,GPG%ALLBC(IL)%TFIXED)
   HFIXEDHI = HFIXEDHI + G%YI(ISPEC,IZ,IX,IY) * HOFT(ISPEC,GPG%ALLBC(IH)%TFIXED)
ENDDO
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%HFIXED  = HFIXEDLO  + F * (HFIXEDHI - HFIXEDLO )

! Set gas enthalpy:
HFIXEDLO = GPROP%CPG * ( GPG%ALLBC(IL)%TFIXEDG - GPG%TDATUM)
HFIXEDHI = GPROP%CPG * ( GPG%ALLBC(IH)%TFIXEDG - GPG%TDATUM)
G%GPYRO_BOUNDARY_CONDITION(IZ,IX,IY,IOR)%HFIXEDG = HFIXEDLO + F * (HFIXEDHI - HFIXEDLO )

! *****************************************************************************
END SUBROUTINE INTERPOLATE_BCS
! *****************************************************************************

! *****************************************************************************
END SUBROUTINE GET_ALL_BOUNDARY_CONDITIONS
! *****************************************************************************

! *****************************************************************************
SUBROUTINE ADJUST_HC_FOR_BLOWING(IZ,IX,IY,NCELLZ,NCELLX,NCELLY,CTYPE,HC0,HC0BLOWING)
! *****************************************************************************

INTEGER, INTENT(IN) :: IZ, IX, IY, NCELLZ, NCELLX, NCELLY
CHARACTER(1), INTENT(IN) :: CTYPE
REAL(EB), INTENT(IN) :: HC0
REAL(EB), INTENT(OUT) :: HC0BLOWING
REAL(EB) :: MC

IF (GPG%SOLVE_PRESSURE) THEN
   SELECT CASE(CTYPE)
      CASE('X')
         IF (IX .EQ. 1     ) MC=ABS(G%MDOTPPDARCYW(IZ,IX,IY))
         IF (IX .EQ. NCELLX) MC=ABS(G%MDOTPPDARCYE(IZ,IX,IY))
      CASE('Y')
         IF (IY .EQ. 1     ) MC=ABS(G%MDOTPPDARCYS(IZ,IX,IY))
         IF (IY .EQ. NCELLY) MC=ABS(G%MDOTPPDARCYN(IZ,IX,IY))
      CASE('Z')
         IF (IZ .EQ. 1     ) MC=ABS(G%MDOTPPDARCYT(IZ,IX,IY))
         IF (IZ .EQ. NCELLZ) MC=ABS(G%MDOTPPDARCYB(IZ,IX,IY))
   END SELECT
ELSE
   SELECT CASE(CTYPE)
      CASE('X')
         WRITE(*,*) 'Cannot adjust x-direction hc for blowing with SOLVE_PRESSURE = F' 
      CASE('Y')
         WRITE(*,*) 'Cannot adjust y-direction hc for blowing with SOLVE_PRESSURE = F' 
      CASE('Z')
         MC=ABS(G%MDOTPPZ(0,IZ,IX,IY))     
   END SELECT
ENDIF

MC         = MAX(MC*GPROP%CPG,1D-10)
HC0BLOWING = MC / (EXP(MC/HC0) - 1D0)
 
! *****************************************************************************
END SUBROUTINE 
! *****************************************************************************


! *****************************************************************************
SUBROUTINE GET_BC_INFO(IX,IZ,IY,NBCPATCH,PATCH_BCLOC,PATCH_IX1,PATCH_IX2,PATCH_IZ1, &
           PATCH_IZ2,PATCH_IY1,PATCH_IY2,PATCH_IPNEXT,PATCH_T,BCLOC,TI,IP,IPNEXT,F)
! *****************************************************************************
INTEGER, INTENT(IN) :: IX,IZ,IY,NBCPATCH
REAL(EB), INTENT(IN) :: TI
CHARACTER(1), INTENT(IN) :: BCLOC
INTEGER, INTENT(IN), DIMENSION (:) :: PATCH_IX1,PATCH_IX2,PATCH_IZ1,PATCH_IZ2,PATCH_IY1,PATCH_IY2,PATCH_IPNEXT !prevents array temporary
REAL(EB), INTENT(IN), DIMENSION (:) :: PATCH_T
CHARACTER(1), INTENT(IN), DIMENSION (:) :: PATCH_BCLOC


INTEGER, INTENT(OUT) :: IP,IPNEXT
REAL(EB), INTENT(OUT) :: F

DO IP = 1, NBCPATCH
   IF (PATCH_BCLOC(IP) .NE. BCLOC) CYCLE
   IPNEXT = PATCH_IPNEXT(IP)

   IF (IPNEXT .NE. -1) THEN
      IF (TI .GT. PATCH_T(IPNEXT)) CYCLE
   ENDIF
         
   SELECT CASE (BCLOC)
      CASE ('T','B')
         IF (IX .LT. PATCH_IX1(IP)) CYCLE
         IF (IX .GT. PATCH_IX2(IP)) CYCLE
         IF (IY .LT. PATCH_IY1(IP)) CYCLE
         IF (IY .GT. PATCH_IY2(IP)) CYCLE         
      CASE ('E','W')
         IF (IZ .LT. PATCH_IZ1(IP)) CYCLE
         IF (IZ .GT. PATCH_IZ2(IP)) CYCLE
         IF (IY .LT. PATCH_IY1(IP)) CYCLE
         IF (IY .GT. PATCH_IY2(IP)) CYCLE
      CASE ('N','S')
         IF (IZ .LT. PATCH_IZ1(IP)) CYCLE
         IF (IZ .GT. PATCH_IZ2(IP)) CYCLE
         IF (IX .LT. PATCH_IX1(IP)) CYCLE
         IF (IX .GT. PATCH_IX2(IP)) CYCLE
   END SELECT

   IF (IPNEXT .EQ. -1) THEN
      F = 1D0
   ELSE
      F = 1D0 - (TI-PATCH_T(IP)) / (PATCH_T(IPNEXT)-PATCH_T(IP))
   ENDIF
                  
   RETURN
ENDDO 

! *****************************************************************************            
END SUBROUTINE GET_BC_INFO
! *****************************************************************************

! *****************************************************************************
SUBROUTINE SET_BC_FIXED_VALUE(AP,AT,AB,AE,AW,AN,AS,B,FIXEDVAL)
! *****************************************************************************
REAL(EB), INTENT(OUT) :: AP,AT,AB,AE,AW,AN,AS,B
REAL(EB), INTENT(IN)  :: FIXEDVAL
      
AP = 1D0
AB = 0D0
AT = 0D0
AE = 0D0
AW = 0D0
AN = 0D0
AS = 0D0
B  = FIXEDVAL

! *****************************************************************************
END SUBROUTINE SET_BC_FIXED_VALUE
! *****************************************************************************

! *****************************************************************************
END MODULE GPYRO_BC
! *****************************************************************************