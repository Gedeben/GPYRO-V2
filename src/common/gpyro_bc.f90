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
INTEGER :: IL, IH, ID_SURF, ICOUNT, IMESH_FDS
REAL(EB) :: F
TYPE (GPYRO_BOUNDARYS_INFORMATION), POINTER :: BOUNDARYS

G         => GPM(IMESH)
BOUNDARYS => GP_BOUDARYS(IMESH)

!!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICOUNT,IMESH_FDS,ID_SURF,IL,IH,F) &
!!$OMP SHARED(NGPYRO_FACES_NEEDING_BCS, BOUNDARYS, IGPYRO_TYPE, IMESH,G,T)
DO ICOUNT = 1, NGPYRO_FACES_NEEDING_BCS(IMESH)
   IMESH_FDS = BOUNDARYS%IMESH_FDS(ICOUNT)
   IF ((IGPYRO_TYPE .EQ. 2) .AND. IMESH_FDS.NE. 0) CYCLE  ! If this is à FDS simulation, and this cell have FDS BC skip 
   IF (BOUNDARYS%IMESH_GPYRO(ICOUNT) .NE. IMESH) CYCLE
   ID_SURF   = BOUNDARYS%SURF_IDX (ICOUNT)
   
   CALL GET_BC_F    (ID_SURF,T,IL,IH,F)
   CALL INTERPOLATE_BCS(ICOUNT,IL,IH,F)
   IF (GPG%BLOWING) CALL ADJUST_HC_FOR_BLOWING(ICOUNT)

ENDDO
!!$OMP END PARALLEL DO


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
SUBROUTINE INTERPOLATE_BCS(ICOUNT,IL,IH,F)
! *****************************************************************************

INTEGER, INTENT(IN) :: ICOUNT,IL,IH
REAL(EB), INTENT(IN) :: F
INTEGER :: ISPEC, IX,IY,IZ
REAL(EB) :: HFIXEDLO, HFIXEDHI

! Linearly interpolate real quantities:
BOUNDARYS%GPBC(ICOUNT)%QE      = GPG%ALLBC(IL)%QE      + F * (GPG%ALLBC(IH)%QE      - GPG%ALLBC(IL)%QE      )
BOUNDARYS%GPBC(ICOUNT)%HC0     = GPG%ALLBC(IL)%HC      + F * (GPG%ALLBC(IH)%HC      - GPG%ALLBC(IL)%HC      )
BOUNDARYS%GPBC(ICOUNT)%TINF    = GPG%ALLBC(IL)%TINF    + F * (GPG%ALLBC(IH)%TINF    - GPG%ALLBC(IL)%TINF    )
BOUNDARYS%GPBC(ICOUNT)%TFIXED  = GPG%ALLBC(IL)%TFIXED  + F * (GPG%ALLBC(IH)%TFIXED  - GPG%ALLBC(IL)%TFIXED  )
BOUNDARYS%GPBC(ICOUNT)%PRES    = GPG%ALLBC(IL)%PRES    + F * (GPG%ALLBC(IH)%PRES    - GPG%ALLBC(IL)%PRES    )
BOUNDARYS%GPBC(ICOUNT)%MFLUX   = GPG%ALLBC(IL)%MDOTPP  + F * (GPG%ALLBC(IH)%MDOTPP  - GPG%ALLBC(IL)%MDOTPP  )
BOUNDARYS%GPBC(ICOUNT)%QEG     = GPG%ALLBC(IL)%QEG     + F * (GPG%ALLBC(IH)%QEG     - GPG%ALLBC(IL)%QEG     )
BOUNDARYS%GPBC(ICOUNT)%HC0G    = GPG%ALLBC(IL)%HCG     + F * (GPG%ALLBC(IH)%HCG     - GPG%ALLBC(IL)%HCG     )
BOUNDARYS%GPBC(ICOUNT)%TINFG   = GPG%ALLBC(IL)%TINFG   + F * (GPG%ALLBC(IH)%TINFG   - GPG%ALLBC(IL)%TINFG   )
BOUNDARYS%GPBC(ICOUNT)%TFIXEDG = GPG%ALLBC(IL)%TFIXEDG + F * (GPG%ALLBC(IH)%TFIXEDG - GPG%ALLBC(IL)%TFIXEDG )
BOUNDARYS%GPBC(ICOUNT)%PRES    = GPG%ALLBC(IL)%PRES    + F * (GPG%ALLBC(IH)%PRES    - GPG%ALLBC(IL)%PRES    )
BOUNDARYS%GPBC(ICOUNT)%YJINF(1:GPROP%NGSPEC) = GPG%ALLBC(IL)%YJINF(1:GPROP%NGSPEC) + F * (GPG%ALLBC(IH)%YJINF(1:GPROP%NGSPEC) - GPG%ALLBC(IL)%YJINF(1:GPROP%NGSPEC))

! Reradiation is discrete so use value at IL:
BOUNDARYS%GPBC(ICOUNT)%RERAD = GPG%ALLBC(IL)%RERADIATION

BOUNDARYS%GPBC(ICOUNT)%HM0     = GPG%ALLBC(IL)%HM      + F * (GPG%ALLBC(IH)%HM      - GPG%ALLBC(IL)%HM      )

IF (GPG%USE_TORTUOSITY_FACTOR_FOR_FLUX) THEN
   BOUNDARYS%GPBC(ICOUNT)%HM0 = GPG%TORTUOSITY_FACTOR * BOUNDARYS%GPBC(ICOUNT)%HM0
ENDIF
IF (BOUNDARYS%GPBC(ICOUNT)%HM0 .LT. EPSILON_FB) BOUNDARYS%GPBC(ICOUNT)%HM0=EPSILON_FB


! *****************************************************************************
END SUBROUTINE INTERPOLATE_BCS
! *****************************************************************************


! *****************************************************************************
SUBROUTINE ADJUST_HC_FOR_BLOWING(ICOUNT)
! *****************************************************************************
INTEGER, INTENT(IN) :: ICOUNT
INTEGER :: IZ, IX, IY,IOR
REAL(EB) :: HC0, HC0BLOWING, MC

IZ  = BOUNDARYS%IZ_GPYRO (ICOUNT)
IX  = BOUNDARYS%IX_GPYRO (ICOUNT)
IY  = BOUNDARYS%IY_GPYRO (ICOUNT)
IOR = BOUNDARYS%IOR_GPYRO(ICOUNT)

IF (GPG%SOLVE_PRESSURE) THEN
   SELECT CASE(IOR)
      CASE( 3) ; MC=ABS(G%MDOTPPDARCYT(IZ,IX,IY))
      CASE(-3) ; MC=ABS(G%MDOTPPDARCYB(IZ,IX,IY))
      CASE( 1) ; MC=ABS(G%MDOTPPDARCYE(IZ,IX,IY))
      CASE(-1) ; MC=ABS(G%MDOTPPDARCYW(IZ,IX,IY))
      CASE( 2) ; MC=ABS(G%MDOTPPDARCYN(IZ,IX,IY))
      CASE(-2) ; MC=ABS(G%MDOTPPDARCYS(IZ,IX,IY))
   END SELECT

ELSE ! Instantaneous gas transport
   IF (IOR .NE. 3) RETURN
   MC=ABS(G%MDOTPPZ(0,IZ,IX,IY))
ENDIF

HC0        = BOUNDARYS%GPBC(ICOUNT)%HC0  ! Default value computed from boundary type
MC         = MAX(MC*GPROP%CPG,1D-10)
HC0BLOWING = MC / (EXP(MC/HC0) - 1D0)
BOUNDARYS%GPBC(ICOUNT)%HC0 = HC0BLOWING

! *****************************************************************************
END SUBROUTINE ADJUST_HC_FOR_BLOWING
! *****************************************************************************


! *****************************************************************************
END SUBROUTINE GET_ALL_BOUNDARY_CONDITIONS
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