PROGRAM GPYRO_STANDALONE

#ifndef GITHASH_PP
#define GITHASH_PP "unknown"
#endif
#ifndef GITDATE_PP
#define GITDATE_PP "unknown"
#endif
#ifndef BUILDDATE_PP
#define BUILDDATE_PP "unknown"
#endif

USE PREC
USE GPYRO_PYRO
USE GPYRO_VARS
USE GPYRO_IO
USE GPYRO_INIT
USE GPYRO_FUNCS
USE OMP_UTILS

IMPLICIT NONE

INTEGER  :: ICASE,IMESH,LO10
REAL(EB) :: TI
REAL(EB) :: TSTART,TEND, TMID1,TMID2, TREAD
LOGICAL  :: LOPEN, FAILED_TS, END_SIMULATION, LAST_TIME_STEP
CHARACTER(300) :: MESSAGE

WRITE(0,'(A)') ""
WRITE(0,'(A)') "        (      (   '  )      (      (         )(    "
WRITE(0,'(A)') "      ()/)   (()     ((       )    )/(      ()/(    "
WRITE(0,'(A)') "   /█(▌▌))  (█(███)  (█)    (▌) (████(▖)  (/ ██\\)  "
WRITE(0,'(A)') "  |█(  __   |█|__)█|  (██\_/█)  |█(__)▌|  |█|  |█|  "
WRITE(0,'(A)') "  |█( |██|  |█ ███/    \███/    |█████/   |█|  |█|  "
WRITE(0,'(A)') "  |█(__|█|  |█|         |█|     |█| \█\   |█|__|█|  "
WRITE(0,'(A)') "   \█████|  |█|          █|     |█|  \█\   \████/   "
WRITE(0,'(A)') "                                                    "
WRITE(0,'(A)') '             STANDALONE SIMULATION MODE             '
WRITE(0,'(A)') "                                                    "

WRITE(0,'(A,A)')      ' Revision         : ',TRIM(GITHASH_PP)
WRITE(0,'(A,A)')      ' Revision Date    : ',TRIM(GITDATE_PP)
WRITE(0,'(A,A)')      ' Compilation Date : ',TRIM(BUILDDATE_PP)
WRITE(0,'(1X,A)') 'https://github.com/reaxfire/gpyro'
WRITE(0,'(A)') ""
WRITE(0, '(A)') OMP_STATUS_MESSAGE()
WRITE(0,'(A)') ""


! Initialize system clock
CALL SYSTEM_CLOCK(COUNT_RATE=CLOCK_COUNT_RATE)
CALL SYSTEM_CLOCK(COUNT_MAX=CLOCK_COUNT_MAX)

GPG%GPYRO_WALL_CLOCK_START = WALL_CLOCK_TIME_THANKS_FDS()

CALL GET_CPU_TIME(TSTART)

IGPYRO_TYPE = 1 !Standalone implementation

CALL GET_CPU_TIME(TMID1)
CALL READ_GPYRO
CALL CHECK_GPYRO(0) !Checks for potential conditions that could lead to seg fault or cause other problems
CALL GET_CPU_TIME(TMID2) ; TREAD = TMID2 - TMID1

DO ICASE = 1, GPG%NCASES

   CALL GET_CPU_TIME(TMID1)
   IMESH = GPG%IMESH(ICASE) ! One cases can only have one mesh for now
   CALL ALLOCATE_GPYRO(IMESH)
   G=>GPM(IMESH)
   GPBCP(1:,1:,1:,-3:)=>G%GPYRO_BOUNDARY_CONDITION(1:,1:,1:,-3:)
   CALL INIT_GPYRO(ICASE,IMESH)
   GPG%TUSED(51) = TREAD
   !If summary file is open, close it:
   INQUIRE(UNIT=LUSUM,OPENED=LOPEN) ; IF(LOPEN) CLOSE(LUSUM)
   CALL GET_CPU_TIME(TMID2) ; GPG%TUSED(52) = GPG%TUSED(52) + TMID2 - TMID1

   CALL DUMP_GPYRO(IMESH,ICASE,0D0)
   IF (GPG%FDSMODE) CALL DUMP_FDS_INPUT_FILE(ICASE)

   GPG%DT     = GPG%DT0
   TI         = GPG%DT0
   G%NTIMESTEPS = 0

   WRITE(0,'(1X,A)')
   WRITE(0,'(1X,A)') 'Entering main timestepping loop.' 

   END_SIMULATION = .FALSE.
   LAST_TIME_STEP = .FALSE.

   DO WHILE (.NOT. END_SIMULATION)
      IF (LAST_TIME_STEP) END_SIMULATION = .TRUE.
      FAILED_TS = .TRUE.
      DO WHILE (FAILED_TS)
         GPG%DT = GPG%DTNEXT
         IF (GPG%ZEROD(ICASE)) THEN
            CALL TG_DRIVER(ICASE,TI,FAILED_TS)
         ELSE
            CALL GPYRO_PYROLYSIS(IMESH,ICASE,G%NCELLZ,G%NCELLX,G%NCELLY,TI,FAILED_TS)
            IF(G%THICKNESS(1,1) .LT. 5D-4*G%ZDIM) THEN 
               TI = GPG%TSTOP(ICASE) + 1.
               CONTINUE
            ENDIF
         ENDIF
      ENDDO
     
      G%NTIMESTEPS = G%NTIMESTEPS + 1
      LO10 = LOG10(REAL(MAX(1,ABS(G%NTIMESTEPS)),EB)) !Borrow a little trick from FDS
      IF (MOD(G%NTIMESTEPS,10**LO10) .EQ. 0 .OR. MOD(G%NTIMESTEPS,1000)==0 .OR. LAST_TIME_STEP) THEN
         WRITE(0,'(1X,A,I7,A,F10.4,A)') 'Timestep #: ', G%NTIMESTEPS, '  Time: ', TI, ' s'
         CALL WRITE_DOTOUT_FILE(ICASE,TI,IMESH) 
      ENDIF

      IF (GPG%MAXTMP.NE.GPG%MAXTMP .OR. GPG%MAXTMP.EQ.GPG%POSINF .OR. GPG%MAXTMP.EQ.GPG%NEGINF) CALL WRITE_DOTOUT_FILE(ICASE,TI,IMESH) 
      GPG%DT = GPG%DTNEXT
      TI     = TI + GPG%DT
      IF (TI .GE. (GPG%TSTOP(ICASE)-EPSILON_FB)) THEN
      !    DTP= TI-GPG%TSTOP(ICASE)
      !    TI=GPG%TSTOP(ICASE)
      !    GPG%DT = GPG%DT-DTP
         LAST_TIME_STEP=.TRUE.
      ENDIF

   ENDDO ! TIME LOOP

   CALL GET_CPU_TIME(TEND)
   GPG%TUSED(0) = GPG%TUSED(0) + TEND - TSTART
   
   !CALL GET_CPU_TIME(TMID1)
   CALL DEALLOCATE_GPYRO
   CALL CLOSE_FILES
   !CALL GET_CPU_TIME(TMID2) ; GPG%TUSED(53) = GPG%TUSED(53) + TMID2 - TMID1


ENDDO !ICASE


MESSAGE = 'End of simulation reached successfully. Shutting down.' 
!CALL GET_CPU_TIME(TEND)
!GPG%TUSED(0) = GPG%TUSED(0) + TEND - TSTART
CALL SHUTDOWN_GPYRO(MESSAGE)

END PROGRAM GPYRO_STANDALONE