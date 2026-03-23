! *****************************************************************************
MODULE GPYRO_PYRO
! *****************************************************************************

USE PREC
USE GPYRO_VARS
USE GPYRO_FUNCS
USE GPYRO_IO
USE GPYRO_BC
USE GPYRO_MASS
USE GPYRO_SOLID_THERMAL
USE GPYRO_GAS

IMPLICIT NONE

CONTAINS

! *****************************************************************************
SUBROUTINE GPYRO_PYROLYSIS(IMESH,ICASE,NCELLZ,NCELLX,NCELLY,TI,FAILED_TIMESTEP)
! *****************************************************************************

INTEGER,  INTENT(IN)  :: IMESH,ICASE,NCELLZ,NCELLX,NCELLY
REAL(EB), INTENT(INOUT) :: TI
LOGICAL, INTENT(OUT) :: FAILED_TIMESTEP

! Local variables
INTEGER :: ITER 
LOGICAL :: UPDATE_COEFFS 
REAL(EB) :: TMID1, TMID2
TYPE(SOLVER_CONVERGENCE_INFO), POINTER :: CI


G=>GPM(IMESH)
CI =>G%CONV_INFO

! Calculate timestep
G%DTIME_TMP   = TI - G%TLAST_TMP
G%DTIME_YIS   = TI - G%TLAST_YIS
G%DTIME_YJG   = TI - G%TLAST_YJG
G%DTIME_P     = TI - G%TLAST_P
G%DTIME_HG    = TI - G%TLAST_HG

G%TLAST_TMP = TI
G%TLAST_YIS = TI
G%TLAST_YJG = TI
G%TLAST_P   = TI
G%TLAST_HG  = TI

! Initialize misc variables
FAILED_TIMESTEP = .FALSE.
GPG%NAN         = .FALSE.

! Set initialy convergence = False for equation beeing solve.
CALL INITIALIZATION_SUBITERATION_CONVERGENCE_INFO(IMESH)

CALL GET_CPU_TIME(TMID1)
! Begin by getting boundary condition info. We only need to do this once
! per timestep (not every iteration). 
CALL GET_ALL_BOUNDARY_CONDITIONS(IMESH,TI)
CALL GET_CPU_TIME(TMID2) ; GPG%TUSED(54) = GPG%TUSED(54) + TMID2 - TMID1

! To prevent the solver from crashing.
! Stop when the solution has converged  or the maximum number of sub iterations is reached,
DO WHILE ((CI%ITER .LE. GPG%NTDMA_ITERATIONS) .AND. (.NOT.CI%CONVERGED_ALL) .AND. (.NOT. GPG%NAN))

   CI%CONVERGED_ALL = .TRUE. ! The solvers change this value to False if no convergence
   CI%ITER = CI%ITER + 1
   ITER = CI%ITER 
   UPDATE_COEFFS=.FALSE.;
   IF (MOD((ITER+1),GPG%NCOEFF_UPDATE_SKIP) .EQ. 0) UPDATE_COEFFS = .TRUE.
   IF (ITER .EQ. 1) UPDATE_COEFFS = .TRUE. 

   !===============================================================!
   !============ GET SPECIES CONCENTRATION EVOLUTION ==============!
   !===============================================================!
   
   CALL UPDATE_SOLID_SPECIES_EVOLUTION(IMESH)

   !===============================================================!
   !====================== GAS TRANSPORT ==========================!
   !===============================================================!

   CALL GET_CPU_TIME(TMID1)
   
   ! Get pyrolysate mass flux. This is not where the darcy mass flux is calculated
   IF (SPROP%NRXN .GT. 0 .AND. (.NOT. GPG%SOLVE_PRESSURE)) THEN
      CALL INSTANTANEOUS_GAS_RELASE(IMESH)
   ENDIF
   ! Get D
   IF (UPDATE_COEFFS .AND. GPG%SOLVE_PRESSURE .OR. GPG%SOLVE_GAS_YJ .OR. GPG%SOLVE_GAS_ENERGY) THEN 
      CALL GET_D(IMESH,MIN(GPROP%IBG,GPROP%NGSPEC),MIN(GPROP%IO2,GPROP%NGSPEC))
   ENDIF

   !Update gas density
   IF (UPDATE_COEFFS .AND. GPG%SOLVE_GAS_YJ) CALL CALCULATE_GAS_DENSITY(IMESH)


   IF (GPG%SOLVE_PRESSURE .AND. (.NOT. CI%CONVERGED_P)) THEN
      CALL PRESSURE_SOLVER(IMESH,G%DTIME_P,GPG%ALPHA_P)
      IF (NCELLZ .GT. 1) CALL DARCIAN_MASS_FLUX(IMESH,NCELLZ,NCELLX,NCELLY,'Z')
      IF (NCELLX .GT. 1) CALL DARCIAN_MASS_FLUX(IMESH,NCELLZ,NCELLX,NCELLY,'X')
      IF (NCELLY .GT. 1) CALL DARCIAN_MASS_FLUX(IMESH,NCELLZ,NCELLX,NCELLY,'Y')
   ENDIF

   CALL GET_CPU_TIME(TMID2);  GPG%TUSED(30) = GPG%TUSED(30) + TMID2 - TMID1

   !===============================================================!
   !==================== SOLID ENTHALPY SOLVER ====================!
   !===============================================================!

   CALL GET_CPU_TIME(TMID1)


   IF (.NOT. CI%CONVERGED_TMP)  THEN
      IF (GPG%FIX_DOMAIN_TEMPERATURE) THEN

         CALL GET_CPU_TIME(TMID2)
         G%TPN(:,:,:) = GPG%TAMB + GPG%BETA(ICASE) * TI / 60D0
         CI%CONVERGED_TMP = .TRUE.
      ELSE
         CALL SOLID_ENTHALPY_SOLVER(IMESH,G%DTIME_TMP)
      ENDIF
   ENDIF
   CALL GET_CPU_TIME(TMID2)
   GPG%TUSED(1) = GPG%TUSED(1) + TMID2 - TMID1 ! Total time for solid enthalpy solver

   !===============================================================!
   !================= GAS ENTHALPY SOLVER =========================!
   !===============================================================!
   
   CALL GET_CPU_TIME(TMID1)
   !Calculate volumetric heat transfer coefficient
   IF (UPDATE_COEFFS .AND. (GPG%SOLVE_GAS_ENERGY .AND. (.NOT. GPG%THERMAL_EQUILIBRIUM)) .AND. GPG%HCV .LT. 0D0) THEN
      CALL CALC_HCV(IMESH,NCELLZ,NCELLX,NCELLY)
   ENDIF

   CALL GAS_ENERGY_SOURCE_TERMS(IMESH)


   IF ((GPG%SOLVE_GAS_YJ .OR. GPG%SOLVE_GAS_ENERGY) .AND. ((.NOT. CI%CONVERGED_YJG(0)) &
       .OR. (.NOT. CI%CONVERGED_HG))) &
   CALL CONVECTIVE_DIFFUSIVE_SOLVER(IMESH,GPROP%NGSPEC,GPG%CONV_DIFF_SCHEME,&
                                       G%DTIME_YJG,GPG%ALPHA_YJG, UPDATE_COEFFS)

   CALL GET_CPU_TIME(TMID2);  GPG%TUSED(30) = GPG%TUSED(30) + TMID2 - TMID1
   !===============================================================!
ENDDO !ITER

CI%ITER_ALL = ITER

CALL CHECK_TIME_STEP_CONVERGENCE(IMESH,TI,FAILED_TIMESTEP)

! If this is  property estimation run, dont dump data
IF (IGPYRO_TYPE .EQ. 1) CALL DUMP_GPYRO(IMESH,ICASE,TI) !Standalone
IF (IGPYRO_TYPE .EQ. 2 .AND. (GPG%MYID .EQ. GPG%PROCESS_GPYRO(IMESH) .OR. (.NOT. GPG%USE_MPI))) &
CALL DUMP_GPYRO(IMESH,ICASE,TI)

CALL GET_CPU_TIME(TMID1)
CALL NEW_UPDATE_SOLUTION(IMESH)
CALL GET_CPU_TIME(TMID2);  GPG%TUSED(55) = GPG%TUSED(55) + TMID2 - TMID1
   

! *****************************************************************************
END SUBROUTINE GPYRO_PYROLYSIS
! *****************************************************************************

! *****************************************************************************
SUBROUTINE TG_DRIVER(ICASE,TI,FAILED_TIMESTEP)
! *****************************************************************************

! INTENT IN variables
INTEGER,  INTENT(IN) :: ICASE
REAL(EB), INTENT(INOUT) :: TI
LOGICAL, INTENT(OUT) :: FAILED_TIMESTEP
!REAL(EB), INTENT(IN) :: TMP

! Local variables
REAL(EB) :: DTIME
INTEGER :: ITERN,IMESH
LOGICAL :: CONVERGED 
IMESH = GPG%IMESH(ICASE)
G=>GPM(IMESH)


! Calculate timestep
DTIME       = TI - G%TLAST_YIS
G%DTIME_YIS = DTIME
G%TLAST_YIS = TI
GPG%NAN= .FALSE.


! ! Calculate temperature assuming linear increase
G%TPN(1,1,1) = GPG%TAMB + GPG%BETA(ICASE) * TI / 60D0

ITERN           = 0
CONVERGED       = .FALSE.
FAILED_TIMESTEP = .FALSE.

G%CONV_INFO%CONVERGED_YIS(:) = .FALSE.
G%CONV_INFO%CONVERGED_YJG(:) = .FALSE.

DO WHILE (ITERN .LT. GPG%NTDMA_ITERATIONS .AND. .NOT. CONVERGED)
   ITERN = ITERN + 1
   IF (GPG%NAN) CYCLE

   CALL UPDATE_SOLID_SPECIES_EVOLUTION(IMESH)

   IF (G%CONV_INFO%CONVERGED_YIS(0)) CONVERGED = .TRUE.
   IF (GPG%NAN                     ) CONVERGED = .FALSE.

ENDDO
      
G%CONV_INFO%ITER_ALL = ITERN
IF (G%RDLTZN(1,1,1) .LE. EPSILON_EB) G%CONSUMED(1,1,1) = .TRUE. 

IF (.NOT. CONVERGED) THEN

   CALL ROLLBACK_SOLUTION(1,1,IMESH)
   IF (GPG%DTNEXT .LE. 1D-4*GPG%DT0) THEN
      IF (IGPYRO_TYPE .EQ. 1) WRITE(*,*) "Failed timestep. Dt can't be reduced anymore. t = ",TI
      FAILED_TIMESTEP = .FALSE.
   ELSE
      IF (IGPYRO_TYPE .EQ. 1) WRITE(*,*) 'Failed timestep. Reducing Dt. t = ',TI
      FAILED_TIMESTEP = .TRUE.
      G%TLAST_YIS = TI - DTIME
      GPG%DTNEXT = 0.9D0 * GPG%DTNEXT
      GPG%DTNEXT = MAX(GPG%DTNEXT, 1D-4*GPG%DT0)
      TI = G%TLAST_YIS + GPG%DTNEXT
      G%CONV_INFO%CONVERGED_YIS(:)=.FALSE.

      G%TPN(1,1,1)     = G%TP(1,1,1)
      G%HPN(1,1,1)     = G%HP(1,1,1)

      G%YIN   (:,1,1,1) = G%YI(:,1,1,1)
      G%XIN   (:,1,1,1) = G%XI(:,1,1,1)
      G%RPN   (1,1,1)   = G%RP(1,1,1)
      G%DLTZN (1,1,1)   = G%DLTZ(1,1,1)
      G%RDLTZN(1,1,1)   = G%RP(1,1,1)*G%DLTZ(1,1,1)

      G%CONSUMED(1,1,1) = .FALSE.
      IF (G%RDLTZN(1,1,1) .LE. EPSILON_EB) G%CONSUMED(1,1,1) = .TRUE. 

      RETURN
   ENDIF
ENDIF

! If this is a standlone implementation, call dump routines
IF(IGPYRO_TYPE .EQ. 1 ) CALL DUMP_GPYRO(1,ICASE,TI)

! Update solution
CALL NEW_UPDATE_SOLUTION(IMESH)

! *****************************************************************************
END SUBROUTINE TG_DRIVER
! *****************************************************************************

! *****************************************************************************
SUBROUTINE NEW_UPDATE_SOLUTION(IMESH)
! *****************************************************************************
INTEGER, INTENT(IN) :: IMESH
INTEGER  ::NCELLZ,NCELLX,NCELLY
INTEGER :: IZ, IX, IY
INTEGER :: IOR, ICOUNT

G=>GPM(IMESH)
NCELLZ = G%NCELLZ
NCELLX = G%NCELLX
NCELLY = G%NCELLY
! Switch field pointers safely based on current association

IF (ASSOCIATED(G%TP, G%STORAGE%TP_W1)) THEN
   G%TP => G%STORAGE%TP_W2
   G%TPN => G%STORAGE%TP_W1
   ! G%STORAGE%TP_W1 = HUGE(0._EB) usfull for detecting beug
ELSE IF (ASSOCIATED(G%TP, G%STORAGE%TP_W2)) THEN
   G%TP => G%STORAGE%TP_W1
   G%TPN => G%STORAGE%TP_W2
   !G%STORAGE%TP_W2 = HUGE(0._EB)
END IF

IF (ASSOCIATED(G%HP, G%STORAGE%HP_W1)) THEN
   G%HP => G%STORAGE%HP_W2
   G%HPN => G%STORAGE%HP_W1
   !G%STORAGE%HP_W1 = -HUGE(0._EB)
ELSE IF (ASSOCIATED(G%HP, G%STORAGE%HP_W2)) THEN
   G%HP => G%STORAGE%HP_W1
   G%HPN => G%STORAGE%HP_W2
   !G%STORAGE%HP_W2 = -HUGE(0._EB)
END IF

IF (ASSOCIATED(G%RP, G%STORAGE%RP_W1)) THEN
   G%RP => G%STORAGE%RP_W2
   G%RPN => G%STORAGE%RP_W1
   !G%STORAGE%RP_W1 =-HUGE(0._EB)
ELSE IF (ASSOCIATED(G%RP, G%STORAGE%RP_W2)) THEN
   G%RP => G%STORAGE%RP_W1
   G%RPN => G%STORAGE%RP_W2
   !G%STORAGE%RP_W2 = -HUGE(0._EB)
END IF

IF (ASSOCIATED(G%RSP, G%STORAGE%RSP_W1)) THEN
   G%RSP => G%STORAGE%RSP_W2
   G%RSPN => G%STORAGE%RSP_W1
   !G%STORAGE%RSP_W1 = -HUGE(0._EB)
ELSE IF (ASSOCIATED(G%RSP, G%STORAGE%RSP_W2)) THEN
   G%RSP => G%STORAGE%RSP_W1
   G%RSPN => G%STORAGE%RSP_W2
   !G%STORAGE%RSP_W2 = -HUGE(0._EB)
END IF

IF (ASSOCIATED(G%DLTZ, G%STORAGE%DLTZ_W1)) THEN
   G%DLTZ => G%STORAGE%DLTZ_W2
   G%DLTZN => G%STORAGE%DLTZ_W1
   !G%STORAGE%DLTZ_W1 = -HUGE(0._EB)
ELSE IF (ASSOCIATED(G%DLTZ, G%STORAGE%DLTZ_W2)) THEN
   G%DLTZ => G%STORAGE%DLTZ_W1
   G%DLTZN => G%STORAGE%DLTZ_W2
   !G%STORAGE%DLTZ_W2 = -HUGE(0._EB)
END IF

IF (ASSOCIATED(G%RDLTZ, G%STORAGE%RDLTZ_W1)) THEN
   G%RDLTZ => G%STORAGE%RDLTZ_W2
   G%RDLTZN => G%STORAGE%RDLTZ_W1
   !G%STORAGE%RDLTZ_W1 = -HUGE(0._EB)
ELSE IF (ASSOCIATED(G%RDLTZ, G%STORAGE%RDLTZ_W2)) THEN
   G%RDLTZ => G%STORAGE%RDLTZ_W1
   G%RDLTZN => G%STORAGE%RDLTZ_W2
   !G%STORAGE%RDLTZ_W2 = -HUGE(0._EB)
END IF

IF (ASSOCIATED(G%TG, G%STORAGE%TG_W1)) THEN
   G%TG => G%STORAGE%TG_W2
   G%TGN => G%STORAGE%TG_W1
   !G%STORAGE%TG_W1 = HUGE(0._EB)
ELSE IF (ASSOCIATED(G%TG, G%STORAGE%TG_W2)) THEN
   G%TG => G%STORAGE%TG_W1
   G%TGN => G%STORAGE%TG_W2
   !G%STORAGE%TG_W2 = HUGE(0._EB)
END IF

IF (ASSOCIATED(G%HG, G%STORAGE%HG_W1)) THEN
   G%HG => G%STORAGE%HG_W2
   G%HGN => G%STORAGE%HG_W1
   !G%STORAGE%HG_W1 = HUGE(0._EB)
ELSE IF (ASSOCIATED(G%HG, G%STORAGE%HG_W2)) THEN
   G%HG => G%STORAGE%HG_W1
   G%HGN => G%STORAGE%HG_W2
   !G%STORAGE%HG_W2 = HUGE(0._EB)
END IF

IF (ASSOCIATED(G%P, G%STORAGE%P_W1)) THEN
   G%P => G%STORAGE%P_W2
   G%PN => G%STORAGE%P_W1
   !G%STORAGE%P_W1 = HUGE(0._EB)
ELSE IF (ASSOCIATED(G%P, G%STORAGE%P_W2)) THEN
   G%P => G%STORAGE%P_W1
   G%PN => G%STORAGE%P_W2
   !G%STORAGE%P_W2 = HUGE(0._EB)
END IF

IF (ASSOCIATED(G%M, G%STORAGE%M_W1)) THEN
   G%M => G%STORAGE%M_W2
   G%MN => G%STORAGE%M_W1
   !G%STORAGE%M_W1 = HUGE(0._EB)
ELSE IF (ASSOCIATED(G%M, G%STORAGE%M_W2)) THEN
   G%M => G%STORAGE%M_W1
   G%MN => G%STORAGE%M_W2
   !G%STORAGE%M_W2 = HUGE(0._EB)
END IF

IF (ASSOCIATED(G%POROSS, G%STORAGE%POROSS_W1)) THEN
   G%POROSS => G%STORAGE%POROSS_W2
   G%POROSSN => G%STORAGE%POROSS_W1
   !G%STORAGE%POROSS_W1 = HUGE(0._EB)
ELSE IF (ASSOCIATED(G%POROSS, G%STORAGE%POROSS_W2)) THEN
   G%POROSS => G%STORAGE%POROSS_W1
   G%POROSSN => G%STORAGE%POROSS_W2
   !G%STORAGE%POROSS_W2 = HUGE(0._EB)
END IF

IF (ASSOCIATED(G%YI, G%STORAGE%YI_W1)) THEN
   G%YI => G%STORAGE%YI_W2
   G%YIN => G%STORAGE%YI_W1
   !G%STORAGE%YI_W1 = HUGE(0._EB)
ELSE IF (ASSOCIATED(G%YI, G%STORAGE%YI_W2)) THEN
   G%YI => G%STORAGE%YI_W1
   G%YIN => G%STORAGE%YI_W2
   !G%STORAGE%YI_W2 = HUGE(0._EB)
END IF

IF (ASSOCIATED(G%XI, G%STORAGE%XI_W1)) THEN
   G%XI => G%STORAGE%XI_W2
   G%XIN => G%STORAGE%XI_W1
   !G%STORAGE%XI_W1 = HUGE(0._EB)
ELSE IF (ASSOCIATED(G%XI, G%STORAGE%XI_W2)) THEN
   G%XI => G%STORAGE%XI_W1
   G%XIN => G%STORAGE%XI_W2
   !G%STORAGE%XI_W2 = HUGE(0._EB)
END IF

IF (ASSOCIATED(G%RYIDZP, G%STORAGE%RYIDZP_W1)) THEN
   G%RYIDZP => G%STORAGE%RYIDZP_W2
   G%RYIDZPN => G%STORAGE%RYIDZP_W1
   !G%STORAGE%RYIDZP_W1 = HUGE(0._EB)
ELSE IF (ASSOCIATED(G%RYIDZP, G%STORAGE%RYIDZP_W2)) THEN
   G%RYIDZP => G%STORAGE%RYIDZP_W1
   G%RYIDZPN => G%STORAGE%RYIDZP_W2
   !G%STORAGE%RYIDZP_W2 = HUGE(0._EB)
END IF

IF (ASSOCIATED(G%RYIDZSIGMA, G%STORAGE%RYIDZSIGMA_W1)) THEN
   G%RYIDZSIGMA => G%STORAGE%RYIDZSIGMA_W2
   G%RYIDZSIGMAN => G%STORAGE%RYIDZSIGMA_W1
   !G%STORAGE%RYIDZSIGMA_W1 = HUGE(0._EB)
ELSE IF (ASSOCIATED(G%RYIDZSIGMA, G%STORAGE%RYIDZSIGMA_W2)) THEN
   G%RYIDZSIGMA => G%STORAGE%RYIDZSIGMA_W1
   G%RYIDZSIGMAN => G%STORAGE%RYIDZSIGMA_W2
   !G%STORAGE%RYIDZSIGMA_W2 = HUGE(0._EB)
END IF

IF (ASSOCIATED(G%YJG, G%STORAGE%YJG_W1)) THEN
   G%YJG => G%STORAGE%YJG_W2
   G%YJGN => G%STORAGE%YJG_W1
   !G%STORAGE%YJG_W1 = HUGE(0._EB)
ELSE IF (ASSOCIATED(G%YJG, G%STORAGE%YJG_W2)) THEN
   G%YJG => G%STORAGE%YJG_W1
   G%YJGN => G%STORAGE%YJG_W2
   !G%STORAGE%YJG_W2 = HUGE(0._EB)
END IF



!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ,IOR,ICOUNT) &
!$OMP SHARED(NGPYRO_FACES_NEEDING_BCS, IMESH, GP_BOUDARYS,GPBCP)
DO ICOUNT = 1, NGPYRO_FACES_NEEDING_BCS(IMESH)
   IZ  = GP_BOUDARYS(IMESH)%IZ_GPYRO (ICOUNT)
   IX  = GP_BOUDARYS(IMESH)%IX_GPYRO (ICOUNT)
   IY  = GP_BOUDARYS(IMESH)%IY_GPYRO (ICOUNT)
   IOR = GP_BOUDARYS(IMESH)%IOR_GPYRO (ICOUNT)

   GPBCP(IZ,IX,IY,IOR)%T_SURFACE_OLD = GPBCP(IZ,IX,IY,IOR)%T_SURFACE
END DO
!$OMP END PARALLEL DO

IF (NCELLX .EQ. 1 .AND. NCELLY .EQ. 1) THEN
   IF (G%THICKNESS(1,1) .GT. 0D0) THEN
      G%THICKNESS(1,1) = SUM(G%DLTZ(:,1,1))
   ENDIF
ENDIF


! If the cell is fully depleted and no mass remains, it creates a critical issue for the solver,
! as handling zero-volume cells is currently not supported.
! As a temporary workaround, when the mass in a cell/ mass initial in the cell falls below GPG%EPS,
! chemical reactions are halted in that cell and a minimal residual mass is retained
! to prevent the volume from becoming zero.

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) &
!$OMP DEFAULT(SHARED) COLLAPSE(3)
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ
   IF(G%CONSUMED(IZ,IX,IY) .OR. G%IMASK(IZ,IX,IY) .OR. (.NOT. G%IS_REACTING(0,IZ,IX,IY))) CYCLE
   IF (G%RDLTZ(IZ,IX,IY)/G%RYIDZ0(0,IZ,IX,IY) .LE. GPG%EPS) THEN
      G%CONSUMED(IZ,IX,IY) = .TRUE.
   ENDIF
ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO

! ***************************************************************************** 
 END SUBROUTINE NEW_UPDATE_SOLUTION
! *****************************************************************************



! *****************************************************************************
SUBROUTINE UPDATE_SOLUTION(IMESH)
! *****************************************************************************

INTEGER, INTENT(IN) :: IMESH

INTEGER  ::NCELLZ,NCELLX,NCELLY
INTEGER :: IZ, IX, IY, IGSPEC, ISSPEC
INTEGER :: IOR, ICOUNT

G=>GPM(IMESH)
NCELLZ = G%NCELLZ
NCELLX = G%NCELLX
NCELLY = G%NCELLY

! Update solution
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ,ISSPEC,IGSPEC) &
!$OMP SHARED(G,GPG,GPROP,NCELLY,NCELLX,NCELLZ) COLLAPSE(3)
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ

   G%TP    (IZ,IX,IY) = G%TPN   (IZ,IX,IY)
   G%HP    (IZ,IX,IY) = G%HPN   (IZ,IX,IY)
   G%RP    (IZ,IX,IY) = G%RPN   (IZ,IX,IY)
   G%RSP   (IZ,IX,IY) = G%RSPN  (IZ,IX,IY)
   G%DLTZ  (IZ,IX,IY) = G%DLTZN (IZ,IX,IY)
   G%RDLTZ (IZ,IX,IY) = G%RDLTZN(IZ,IX,IY)

   DO ISSPEC = 1, SPROP%NSSPEC
      G%YI        (ISSPEC,IZ,IX,IY)   = G%YIN        (ISSPEC,IZ,IX,IY)
      G%XI        (ISSPEC,IZ,IX,IY)   = G%XIN        (ISSPEC,IZ,IX,IY)
      G%RYIDZP    (ISSPEC,IZ,IX,IY)   = G%RYIDZPN    (ISSPEC,IZ,IX,IY)
      G%RYIDZSIGMA(ISSPEC,IZ,IX,IY)   = G%RYIDZSIGMAN(ISSPEC,IZ,IX,IY)
   ENDDO

   IF (GPG%SOLVE_GAS_YJ) THEN 
      G%RG(IZ,IX,IY) = G%RGN(IZ,IX,IY)
      DO IGSPEC = 1, GPROP%NGSPEC
         G%YJG(IGSPEC,IZ,IX,IY) = G%YJGN(IGSPEC,IZ,IX,IY)
      ENDDO 
   ENDIF

   IF (GPG%SOLVE_PRESSURE) G%P(IZ,IX,IY) = G%PN(IZ,IX,IY)

   IF (GPG%SOLVE_GAS_YJ .AND. GPG%SOLVE_GAS_ENERGY) THEN
      G%TG(IZ,IX,IY) = G%TGN(IZ,IX,IY)
      G%HG(IZ,IX,IY) = G%HGN(IZ,IX,IY)
   ENDIF

   IF (GPG%SOLVE_GAS_ENERGY .OR. GPG%SOLVE_GAS_YJ .OR. GPG%SOLVE_PRESSURE) THEN
      G%POROSS(IZ,IX,IY) = G%POROSSN(IZ,IX,IY)
   ENDIF

ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO


!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ,IOR,ICOUNT) &
!$OMP SHARED(NGPYRO_FACES_NEEDING_BCS, IMESH, GP_BOUDARYS,GPBCP)
DO ICOUNT = 1, NGPYRO_FACES_NEEDING_BCS(IMESH)
   IZ  = GP_BOUDARYS(IMESH)%IZ_GPYRO (ICOUNT)
   IX  = GP_BOUDARYS(IMESH)%IX_GPYRO (ICOUNT)
   IY  = GP_BOUDARYS(IMESH)%IY_GPYRO (ICOUNT)
   IOR = GP_BOUDARYS(IMESH)%IOR_GPYRO (ICOUNT)

   GPBCP(IZ,IX,IY,IOR)%T_SURFACE_OLD = GPBCP(IZ,IX,IY,IOR)%T_SURFACE
END DO
!$OMP END PARALLEL DO

IF (NCELLX .EQ. 1 .AND. NCELLY .EQ. 1) THEN
   IF (G%THICKNESS(1,1) .GT. 0D0) THEN
      G%THICKNESS(1,1) = SUM(G%DLTZN(:,1,1))
   ENDIF
ENDIF


! If the cell is fully depleted and no mass remains, it creates a critical issue for the solver,
! as handling zero-volume cells is currently not supported.
! As a temporary workaround, when the mass in a cell/ mass initial in the cell falls below GPG%EPS,
! chemical reactions are halted in that cell and a minimal residual mass is retained
! to prevent the volume from becoming zero.

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) &
!$OMP DEFAULT(SHARED) COLLAPSE(3)
DO IY = 1, NCELLY
DO IX = 1, NCELLX
DO IZ = 1, NCELLZ
   IF(G%CONSUMED(IZ,IX,IY) .OR. G%IMASK(IZ,IX,IY) .OR. (.NOT. G%IS_REACTING(0,IZ,IX,IY))) CYCLE
   IF (G%RDLTZN(IZ,IX,IY)/G%RYIDZ0(0,IZ,IX,IY) .LE. GPG%EPS) THEN
      G%CONSUMED(IZ,IX,IY) = .TRUE.
   ENDIF
ENDDO
ENDDO
ENDDO
!$OMP END PARALLEL DO

! *****************************************************************************
END SUBROUTINE UPDATE_SOLUTION
! *****************************************************************************

! *****************************************************************************
SUBROUTINE ROLLBACK_SOLUTION(IXP,IYP,IMESH)
! *****************************************************************************

INTEGER, INTENT(IN) :: IXP,IYP,IMESH
INTEGER :: IOR, ICOUNT, IZ,IX,IY

G=>GPM(IMESH)

G%TPN(:,IXP,IYP)     = G%TP(:,IXP,IYP)
G%HPN(:,IXP,IYP)     = G%HP(:,IXP,IYP)
IF (GPG%SOLVE_PRESSURE) G%PN(:,IXP,IYP) = G%P(:,IXP,IYP)
IF (GPG%SOLVE_GAS_YJ) THEN 
   G%YJGN(:,:,IXP,IYP) = G%YJG(:,:,IXP,IYP)
   G%RGN(:,IXP,IYP)    = G%RG(:,IXP,IYP)
   IF (GPG%SOLVE_GAS_ENERGY) THEN
      G%TGN(:,IXP,IYP)     = G%TG(:,IXP,IYP)
      G%HGN(:,IXP,IYP)     = G%HG(:,IXP,IYP)      
   ENDIF
ENDIF
      
IF (GPG%SOLVE_POROSITY) THEN
   G%RSPN   (:,IXP,IYP) = G%RSP(:,IXP,IYP)
   G%POROSSN(:,IXP,IYP) = G%POROSS(:,IXP,IYP)
ENDIF

G%YIN   (:,:,IXP,IYP) = G%YI(:,:,IXP,IYP)
G%XIN   (:,:,IXP,IYP) = G%XI(:,:,IXP,IYP)
G%RPN   (:,IXP,IYP)   = G%RP(:,IXP,IYP)
G%DLTZN (:,IXP,IYP)   = G%DLTZ(:,IXP,IYP)
G%RDLTZN(:,IXP,IYP)   = G%RP(:,IXP,IYP)*G%DLTZ(:,IXP,IYP)

IF (G%THICKNESS(IXP,IYP) .GT. 0D0) THEN
   G%THICKNESS(IXP,IYP) = SUM(G%DLTZ(:,IXP,IYP))
ENDIF


DO ICOUNT = 1, NGPYRO_FACES_NEEDING_BCS(IMESH)
   IZ  = GP_BOUDARYS(IMESH)%IZ_GPYRO (ICOUNT)
   IX  = GP_BOUDARYS(IMESH)%IX_GPYRO (ICOUNT)
   IY  = GP_BOUDARYS(IMESH)%IY_GPYRO (ICOUNT)
   IOR = GP_BOUDARYS(IMESH)%IOR_GPYRO (ICOUNT)
   GPBCP(IZ,IX,IY,IOR)%T_SURFACE = GPBCP(IZ,IX,IY,IOR)%T_SURFACE_OLD
END DO


G%CONSUMED(:,IXP,IYP) = .FALSE.
IF (GPG%SOLVE_PRESSURE) THEN !Should calculate this from pressure gradient, fix! 
   G%MDOTPPDARCYT(:,IXP,IYP) = 0D0 
   G%MDOTPPDARCYB(:,IXP,IYP) = 0D0
   IF (G%NCELLX .GT. 1) THEN 
      G%MDOTPPDARCYE(:,IXP,IYP) = 0D0
      G%MDOTPPDARCYW(:,IXP,IYP) = 0D0
   ENDIF
   IF (G%NCELLY .GT. 1) THEN 
      G%MDOTPPDARCYS(:,IXP,IYP) = 0D0
      G%MDOTPPDARCYN(:,IXP,IYP) = 0D0
   ENDIF

ENDIF

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) &
!$OMP DEFAULT(SHARED)
DO IZ = 1, G%NCELLZ
   IF(G%CONSUMED(IZ,IXP,IYP) .OR. G%IMASK(IZ,IXP,IYP) .OR. (.NOT. G%IS_REACTING(0,IZ,IXP,IYP))) CYCLE
   IF (G%RDLTZN(IZ,IXP,IYP)/G%RYIDZ0(0,IZ,IXP,IYP) .LE. GPG%EPS) THEN
      G%CONSUMED(IZ,IXP,IYP) = .TRUE.
   ENDIF
ENDDO
!$OMP END PARALLEL DO


! *****************************************************************************
END SUBROUTINE ROLLBACK_SOLUTION
! *****************************************************************************


! *****************************************************************************
SUBROUTINE INITIALIZATION_SUBITERATION_CONVERGENCE_INFO(IMESH)
! *****************************************************************************
INTEGER ,INTENT(IN) :: IMESH
TYPE(SOLVER_CONVERGENCE_INFO), POINTER :: CI


G=>GPM(IMESH)
CI =>G%CONV_INFO
   
! Set convergence variables according to which equations are being solved. 
! If an equation isn't being solved, we automatically set its convergence status to true.

CI%CONVERGED_ALL = .FALSE.
CI%CONVERGED_TMP = .FALSE.

CI%ITER_TMP = 0
CI%ITER_P   = 0
CI%ITER_HG  = 0
CI%ITER_YIS = 0
CI%ITER_YJG = 0
CI%ITER_ALL = 0
CI%ITER     = 0


IF (SPROP%NRXN .EQ. 0 .OR. GPG%NOCONSUMPTION) THEN 
   CI%CONVERGED_YIS(:) =.TRUE.
ELSE
   CI%CONVERGED_YIS(:) =.FALSE.
ENDIF

IF (GPG%SOLVE_GAS_YJ) THEN
   CI%CONVERGED_YJG(:)= .FALSE.
ELSE
   CI%CONVERGED_YJG(:) = .TRUE.
ENDIF

IF (GPG%SOLVE_PRESSURE) THEN 
   CI%CONVERGED_P   = .FALSE.
ELSE
   CI%CONVERGED_P   = .TRUE.
ENDIF

IF (GPG%SOLVE_GAS_ENERGY .AND. (.NOT. GPG%THERMAL_EQUILIBRIUM)) THEN
   CI%CONVERGED_HG = .FALSE.
ELSE 
   CI%CONVERGED_HG = .TRUE.
ENDIF
! *****************************************************************************
END SUBROUTINE INITIALIZATION_SUBITERATION_CONVERGENCE_INFO
! *****************************************************************************

! *****************************************************************************
SUBROUTINE CHECK_TIME_STEP_CONVERGENCE(IMESH,TI,FAILED_TIMESTEP)
! *****************************************************************************

INTEGER,  INTENT(IN)  :: IMESH
REAL(EB), INTENT(INOUT) :: TI
LOGICAL, INTENT(OUT) :: FAILED_TIMESTEP

! Local variables
INTEGER :: IX,IY,ISSPEC,IGSPEC,LUCONVG=1700
LOGICAL :: LOPEN 
CHARACTER(300) :: MESSAGE,FN
CHARACTER(2) :: TWO

CHARACTER(2000) :: WRITESTR, STR
TYPE(SOLVER_CONVERGENCE_INFO), POINTER :: CI
INTEGER :: NCELLZ, NCELLX, NCELLY

G=>GPM(IMESH)
CI =>G%CONV_INFO

NCELLZ = G%NCELLZ
NCELLX = G%NCELLX
NCELLY = G%NCELLY

FAILED_TIMESTEP = .FALSE.


80 FORMAT (A,F9.4, ' s')
81 FORMAT (A,F9.4, ' s')
82 FORMAT (A,F9.4,A,E12.6, ' s')

! Write convergence info:
IF (GPG%DUMP_DETAILED_CONVERGENCE) THEN

  INQUIRE(UNIT=LUCONVG+IMESH,OPENED=LOPEN)
  IF (.NOT. LOPEN) THEN
     WRITE(TWO,'(I2.2)') IMESH
     FN = 'iterations_mesh_' // TWO // '.csv'
     OPEN(UNIT=LUCONVG+IMESH,FILE=TRIM(FN),FORM='FORMATTED',STATUS='REPLACE')
     WRITESTR='t(s),Dt(s),iter,iter_T,'

     IF (GPG%SOLVE_PRESSURE) WRITESTR = TRIM(WRITESTR) // "iter_P,"

     DO ISSPEC = 1, SPROP%NSSPEC
        WRITE(TWO,'(I2.2)') ISSPEC
        WRITESTR = TRIM(WRITESTR) // "iter_YIS(" // TWO // ")," 
     ENDDO

     IF (GPG%SOLVE_GAS_YJ) THEN 
        DO IGSPEC = 1, GPROP%NGSPEC
           WRITE(TWO,'(I2.2)') IGSPEC
           WRITESTR = TRIM(WRITESTR) // "iter_YJG(" // TWO // ")," 
        ENDDO
     ENDIF 

     WRITE(LUCONVG+IMESH,'(A)') TRIM(WRITESTR)
  ENDIF

  WRITE(WRITESTR,'(2(F12.6,","),2(I3,","))') TI, G%DTIME_TMP, CI%ITER, CI%ITER_TMP

  IF (GPG%SOLVE_PRESSURE) THEN 
     WRITE(STR,'(I3,",")') CI%ITER_P
     WRITESTR = TRIM(WRITESTR) // TRIM(STR)
  ENDIF

  DO ISSPEC = 1, SPROP%NSSPEC
     WRITE(STR,'(I3,",")') CI%ITER_YIS(ISSPEC)
     WRITESTR = TRIM(WRITESTR) // TRIM(STR)
  ENDDO

  IF (GPG%SOLVE_GAS_YJ) THEN 
     DO IGSPEC = 1, GPROP%NGSPEC
        WRITE(STR,'(I3,",")') CI%ITER_YJG(IGSPEC)
        WRITESTR = TRIM(WRITESTR) // TRIM(STR)
     ENDDO
  ENDIF

  WRITE(LUCONVG+IMESH,'(A)') TRIM(WRITESTR)

ENDIF

! If All converged no need to do what following
IF (CI%CONVERGED_ALL) RETURN

IF (IGPYRO_TYPE .NE. 3) THEN
   IF (.NOT. CI%CONVERGED_YIS(0)) THEN
      WRITE(*,80) 'Solid Yi not converged @ time: ', TI
   ENDIF

   IF (GPG%SOLVE_GAS_YJ .AND. .NOT. CI%CONVERGED_YJG(0) )THEN
      WRITE(*,80) 'Gas Yj not converged @ time: ', TI
      WRITE(*,* ) 'Gas Yj max residual: ', MAXVAL(G%RESIDUAL_YJG(:,:,:,:))
   ENDIF

   IF (.NOT. CI%CONVERGED_HG) THEN
      WRITE(*,80) 'Gas h not converged @ time: ', TI
   ENDIF

   IF (GPG%SOLVE_PRESSURE .AND. (.NOT. CI%CONVERGED_P))THEN
      WRITE(*,80) 'Pressure not converged @ time: ', TI
      WRITE(*,* ) 'Pressure max residual: ', MAXVAL(G%RESIDUAL_P(:,:,:))
      CONTINUE
   ENDIF

   IF (.NOT. CI%CONVERGED_TMP) THEN
      WRITE(*,80) 'Solid temperature not converged @ time: ', TI
      WRITE(*,* ) 'Temperature max residual: ', MAXVAL(G%RESIDUAL_TMP(:,:,:))
   ENDIF
ENDIF

IF (GPG%DTNEXT .LE. GPG%DTMIN_KILL .AND. IGPYRO_TYPE .NE. 2) THEN
   IF (IGPYRO_TYPE .EQ. 1) THEN
      FAILED_TIMESTEP = .TRUE.
      WRITE(*,81) "Failed timestep.  Timestep can't be reduced anymore at t = ", TI
      MESSAGE='Shutting down because timestep is smaller than DTMIN_KILL.'
      CALL SHUTDOWN_GPYRO(MESSAGE)
   ENDIF
   ELSE
   IF (IGPYRO_TYPE .NE. 2) THEN
      IF (GPG%NAN .AND. IGPYRO_TYPE .EQ. 1) WRITE(*,*) 'NaN, +Infinity, or -Infinity encountered.' 
      IF (IGPYRO_TYPE .EQ. 1) THEN 
         WRITE(*,82) 'Failed timestep at t = ', TI, '. Reducing timestep to ', 0.5D0 * GPG%DTNEXT 
      ENDIF
      FAILED_TIMESTEP = .TRUE.
      G%TLAST_TMP = TI - G%DTIME_TMP
      G%TLAST_YIS = TI - G%DTIME_YIS
      G%TLAST_YJG = TI - G%DTIME_YJG
      G%TLAST_P   = TI - G%DTIME_P
      G%TLAST_HG  = TI - G%DTIME_HG

      GPG%DTNEXT = 0.5D0 * GPG%DTNEXT
      !GPG%DTNEXT = MAX(GPG%DTNEXT, 1D-4*GPG%DT0)
      TI = G%TLAST_TMP + GPG%DTNEXT

      DO IY = 1, NCELLY
      DO IX = 1, NCELLX
         CALL ROLLBACK_SOLUTION(IX,IY,IMESH)
      ENDDO
      ENDDO
   ELSE !FDS
   
      FAILED_TIMESTEP = .TRUE.
      G%TLAST_TMP = TI - G%DTIME_TMP
      G%TLAST_YIS = TI - G%DTIME_YIS
      G%TLAST_YJG = TI - G%DTIME_YJG
      G%TLAST_P   = TI - G%DTIME_P
      G%TLAST_HG  = TI - G%DTIME_HG
      
      DO IY = 1, NCELLY
      DO IX = 1, NCELLX
         CALL ROLLBACK_SOLUTION(IX,IY,IMESH)
      ENDDO
      ENDDO
   
   ENDIF
   
   RETURN
ENDIF



GPG%DTNEXT = MIN(1.001*GPG%DTNEXT,GPG%DT0)


! *****************************************************************************
END SUBROUTINE CHECK_TIME_STEP_CONVERGENCE
! *****************************************************************************

! *****************************************************************************
END MODULE GPYRO_PYRO
! *****************************************************************************
