! *****************************************************************************
MODULE GPYRO_MESH
! *****************************************************************************

USE PREC
USE GPYRO_VARS
USE GPYRO_FUNCS, ONLY :  GET_CPU_TIME
USE GPYRO_IO, ONLY: SHUTDOWN_GPYRO

IMPLICIT NONE

CONTAINS


! *****************************************************************************
SUBROUTINE UPDATE_MESH(IMESH)
! *****************************************************************************
INTEGER, INTENT(IN)  :: IMESH
INTEGER              :: IX,IY,IZ
INTEGER              :: NCELLX,NCELLY,NCELLZ
INTEGER              :: IOR, ICOUNT
REAL(EB)             :: dl1,dl2
TYPE(GPYRO_MESH_TYPE), POINTER :: M
TYPE (GPYRO_BOUNDARYS_INFORMATION), POINTER :: BOUNDARYS

! Mesh pointers
G => GPM(IMESH)
M => G%MESH
BOUNDARYS => GP_BOUDARYS(IMESH)

NCELLX= M%NCELLX
NCELLY= M%NCELLY
NCELLZ= M%NCELLZ

IF (M%DIMENSION .EQ. 0) THEN  !OD
   M%DZN (1,1,1) = G%MASS_N(1,1,1) / G%RPN(1,1,1)
   M%DV => M%DZN
   RETURN
ENDIF


IF (GPG%SOLVER_DEFORMATION_MODE .EQ. 1) THEN
   !no deformation
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) &
   !$OMP  COLLAPSE(3)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
      M%DZN (IZ,IX,IY) = M%DZ(IZ,IX,IY) 
      M%DV  (IZ,IX,IY) = G%MASS_N(IZ,IX,IY) / G%RPN(IZ,IX,IY)
   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO
   RETURN
ENDIF

IF (GPG%SOLVER_DEFORMATION_MODE .GE. 2) THEN !
   !2 Deformation But Flux computed as the initial non deformed mesh
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ) &
   !$OMP  COLLAPSE(3)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
      IF (M%IMASK(IZ,IX,IY)) CYCLE
         M%DV (IZ,IX,IY) = G%MASS_N(IZ,IX,IY) / G%RPN(IZ,IX,IY)
         M%DZN(IZ,IX,IY) = M%DV(IZ,IX,IY) / M%DXDY(IZ,IX,IY) 
   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   !Cell size at the boundary
   DO ICOUNT = 1, NGPYRO_FACES_NEEDING_BCS(IMESH)
      IF(.NOT. BOUNDARYS%COMPLETE_CELL_AT_BC(ICOUNT)) CYCLE
   
      IZ  = BOUNDARYS%IZ_GPYRO (ICOUNT)
      IX  = BOUNDARYS%IX_GPYRO (ICOUNT)
      IY  = BOUNDARYS%IY_GPYRO (ICOUNT)
      IOR = BOUNDARYS%IOR_GPYRO(ICOUNT)
      SELECT CASE(IOR)
         CASE( 3) ; M%DZN(IZ-1,IX,IY) = M%DZN(IZ,IX,IY)
         CASE(-3) ; M%DZN(IZ+1,IX,IY) = M%DZN(IZ,IX,IY)
      END SELECT
   ENDDO
   
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IZ,IX,IY) &
   !$OMP  COLLAPSE(2)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
      ! Calculate new z values from new Dz
      !$OMP SIMD
      DO IZ = 3, NCELLZ-1
         dl1=M%DZN(IZ  ,IX,IY)
         dl2=M%DZN(IZ-1,IX,IY)
         M%DZT(IZ,IX,IY) = 0.5D0 * (dl1 + dl2)
      ENDDO

      M%DZT(2     ,IX,IY) = M%DZN(1     ,IX,IY) + 0.5D0 * M%DZN(2       ,IX,IY)
      M%DZT(NCELLZ,IX,IY) = M%DZN(NCELLZ,IX,IY) + 0.5D0 * M%DZN(NCELLZ-1,IX,IY)
      M%DZT(1     ,IX,IY) = M%DZT  (2     ,IX,IY)

      DO IZ = 1, NCELLZ - 1
         M%DZB(IZ,IX,IY) = M%DZT(IZ+1,IX,IY)
      ENDDO
      M%DZB(NCELLZ,IX,IY) = M%DZB(NCELLZ-1,IX,IY)

   ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDIF

IF (GPG%SOLVER_DEFORMATION_MODE .GE. 3) THEN
   !! Deformation And Heat Flux acount for mesh distortient
   !! WARNING: 2D/3D deformation is currently not supported; implementation ongoing…

   CALL COMPUTE_Z_CELL_POSITION(IMESH)

   IF (NCELLX .GT. 1) THEN
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ,dl1,dl2)&
      !$OMP DEFAULT(SHARED) COLLAPSE(3)
      DO IY = 1, NCELLY
      DO IX = 2, NCELLX
      DO IZ = 1, NCELLZ
         dl1=0.5D0*(M%DX(IZ,IX,IY)+M%DX(IZ,IX-1,IY))
         dl2=       M%Z (IZ,IX,IY)-M%Z (IZ,IX-1,IY)

         ! Tangent of the angle alpha formed by the line connecting the two cell centers and the x-axis
         M%TAN_X(IZ,IX,IY) = dl2 / dl1

         ! Using Pythagoras to compute the distance between cell centers
         M%DXW  (IZ,IX,IY) = (dl1**2+dl2**2)**0.5D0

         ! The face surface is assumed to be the average thickness of the surrounding cells times dy
         M%SXW  (IZ,IX,IY) = 0.5D0 * (M%DZN(IZ,IX,IY)+M%DZN(IZ,IX-1,IY))*M%DY(IZ,IX,IY)

      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO


      ! Boundary

      ! Assuming the angle is zero at the boundary
      M%TAN_X(:,1       ,:) = 0D0
      M%TAN_X(:,NCELLX+1,:) = 0D0

      ! The surface at the boundary is assumed to be the cell thickness times dy
      M%SXW  (:,1     ,:) = M%DZN(:,1     ,:)*M%DY(:,1     ,:)
      M%SXE  (:,NCELLX,:) = M%DZN(:,NCELLX,:)*M%DY(:,NCELLX,:)

      ! Boundary cells are treated as half-cells for computing center-to-center distances

      M%DXW  (:,2       ,:) = ( (    M%Z (:,2,:) - M%Z (:,1,:))**2 &
                              +(0.5*M%DX(:,2,:) + M%DX(:,1,:))**2 )**0.5
      M%DXW  (:,1       ,:) = M%DXW(:,2,:)

      M%DXW  (:,NCELLX  ,:) = ( (    M%Z (:,NCELLX-1,:) - M%Z (:,NCELLX,:))**2  &
                              +(0.5*M%DX(:,NCELLX-1,:) + M%DX(:,NCELLX,:))**2 ) **0.5
      
      
      M%DXE  (:,NCELLX,:) = M%DXW(:,NCELLX,:)


      M%SXE   (:,1:NCELLX-1,:) =M%SXW  (:,2:NCELLX,:)
      M%DXE   (:,1:NCELLX-1,:) =M%DXW  (:,2:NCELLX,:)
   ENDIF 

   IF (NCELLY .GT. 1) THEN

   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ,dl1,dl2)&
   !$OMP DEFAULT(SHARED) COLLAPSE(3)
   DO IY = 2, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 1, NCELLZ
      dl1=0.5D0*(M%DY(IZ,IX,IY)+M%DY(IZ,IX,IY-1))
      dl2=       M%Z (IZ,IX,IY)-M%Z (IZ,IX,IY-1)

      ! Tangent of the angle alpha formed by the line connecting the two cell centers and the y-axis
      M%TAN_Y(IZ,IX,IY) = dl2 / dl1

      ! Using Pythagoras to compute the distance between cell centers
      M%DYS  (IZ,IX,IY) = (dl1**2+dl2**2)**0.5D0

      ! The face surface is assumed to be the average thickness of the surrounding cells times dx
      M%SYS  (IZ,IX,IY) = 0.5D0 * (M%DZN(IZ,IX,IY)+M%DZN(IZ,IX,IY-1))*M%DX(IZ,IX,IY)
   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   M%SYN(:,:,1:NCELLY-1)=M%SYS(:,:,2:NCELLY)
   ! Boundary

   ! Assuming the angle is zero at the boundary
   M%TAN_Y(:,:,1       ) = 0D0
   M%TAN_Y(:,:,NCELLY+1) = 0D0

   ! The surface at the boundary is assumed to be the cell thickness times dx
   M%SYS  (:,:,1     ) = M%DZN(:,:,1     )*M%DX(:,:,1     )
   M%SYN  (:,:,NCELLY) = M%DZN(:,:,NCELLY)*M%DX(:,:,NCELLY)

   ! Boundary cells are treated as half-cells for computing center-to-center distances

   M%DYS  (:,:,2       ) = ( (    M%Z (:,:,2) - M%Z (:,:,1))**2 &
                            +(0.5*M%DY(:,:,2) + M%DY(:,:,1))**2 )**0.5
   M%DYS  (:,:,1       ) = M%DYS(:,:,2)

   M%DYS  (:,:,NCELLY  ) = ( (    M%Z (:,:,NCELLY-1) - M%Z (:,:,NCELLY))**2  &
                            +(0.5*M%DY(:,:,NCELLY-1) + M%DY(:,:,NCELLY))**2 ) **0.5
   M%DYN  (:,:,NCELLY) = M%DYS(:,:,NCELLY)

   M%DYN(:,:,1:NCELLY-1)=M%DYS(:,:,2:NCELLY)

   ENDIF
!! WARNING: 2D/3D deformation is currently not supported; implementation ongoing…
ENDIF

CALL CALCULATE_DIFFUSIVE_MESH_FACTORS()
! *****************************************************************************
END SUBROUTINE UPDATE_MESH
! *****************************************************************************



! *****************************************************************************
SUBROUTINE CALCULATE_DIFFUSIVE_MESH_FACTORS()
! *****************************************************************************
! Calculate weighting factors

INTEGER :: NCELLZ,NCELLX,NCELLY
INTEGER :: IZ,IX,IY
TYPE (GPYRO_MESH_TYPE), POINTER :: M
INTEGER :: DEFORMATION_MODEL

DEFORMATION_MODEL=1
M => G%MESH

NCELLZ = M%NCELLZ
NCELLX = M%NCELLX
NCELLY = M%NCELLY

SELECT CASE (GPG%SOLVER_DEFORMATION_MODE)
   CASE(1,2)

      IF (NCELLZ .GT. 1) THEN

         !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ)&
         !$OMP DEFAULT(SHARED) COLLAPSE(3)
         DO IY = 1, NCELLY
         DO IX = 1, NCELLX
         DO IZ = 1, NCELLZ
            M%GZT(IZ,IX,IY) = M%DXDY(IZ,IX,IY)/M%DZT(IZ,IX,IY)
            M%GZB(IZ,IX,IY) = M%DXDY(IZ,IX,IY)/M%DZB(IZ,IX,IY)
         ENDDO
         ENDDO
         ENDDO
         !$OMP END PARALLEL DO
      ENDIF

      IF ((G%NTIMESTEPS .LE. 1) .AND. (NCELLX .GT. 1)) THEN
         ! Only at first iteration, as mesh don't deforme in x direction
         !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ)&
         !$OMP DEFAULT(SHARED) COLLAPSE(3)
         DO IY = 1, NCELLY
         DO IX = 1, NCELLX
         DO IZ = 1, NCELLZ
            M%GXE(IZ,IX,IY) = M%SXE(IZ,IX,IY)/M%DXE(IZ,IX,IY)
            M%GXW(IZ,IX,IY) = M%SXW(IZ,IX,IY)/M%DXW(IZ,IX,IY)
         ENDDO
         ENDDO
         ENDDO
         !$OMP END PARALLEL DO
      ENDIF

      IF ((G%NTIMESTEPS .LE. 1) .AND. (NCELLY .GT. 1)) THEN
         ! Only at first iteration, as mesh don't deforme in y direction
         !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX,IY,IZ)&
         !$OMP DEFAULT(SHARED) COLLAPSE(3)
         DO IY = 1, NCELLY
         DO IX = 1, NCELLX
         DO IZ = 1, NCELLZ
            M%GYN(IZ,IX,IY) = M%SYN(IZ,IX,IY)/M%DYN(IZ,IX,IY)
            M%GYS(IZ,IX,IY) = M%SYS(IZ,IX,IY)/M%DYS(IZ,IX,IY)
         ENDDO
         ENDDO
         ENDDO
         !$OMP END PARALLEL DO
      ENDIF

   CASE(3)
      !BIG PROJECT To be implemented. Good luck!

END SELECT


! *****************************************************************************
END SUBROUTINE CALCULATE_DIFFUSIVE_MESH_FACTORS
! *****************************************************************************

! *****************************************************************************
SUBROUTINE CALCULATE_INTERFACE_MESH_FACTORS()
! *****************************************************************************
! Calculate weighting factors

INTEGER :: NCELLZ,NCELLX,NCELLY,IZ,IX,IY
TYPE (GPYRO_MESH_TYPE), POINTER :: M

M => G%MESH

NCELLZ = M%NCELLZ
NCELLX = M%NCELLX
NCELLY = M%NCELLY

IF (NCELLZ .GT. 1) THEN
    !z-direction:
    G%FB(NCELLZ,:,:) = 0D0
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IZ) SHARED(G,M,NCELLZ)
    DO IZ = 1, NCELLZ-1
    G%FB(IZ,:,:) = M%DZN(IZ+1,:,:) / (M%DZN(IZ+1,:,:) + M%DZN(IZ,:,:))
    ENDDO
    !$OMP END PARALLEL DO 

    G%FT(1,:,:) = 0D0
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IZ) SHARED(G,M,NCELLZ)
    DO IZ = 2, NCELLZ
    G%FT(IZ,:,:) = M%DZN(IZ-1,:,:) / (M%DZN(IZ-1,:,:) + M%DZN(IZ,:,:))
    ENDDO
    !$OMP END PARALLEL DO 
ENDIF

!x-direction:
! As in 2D/3D ther is no deformation in x only calculate at the first time step
IF ((G%NTIMESTEPS .EQ. 1) .AND. (NCELLX .GT. 1)) THEN
   G%FE(:,NCELLX,:) = 0D0
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX) SHARED(G,M,NCELLX)
   DO IX = 1, NCELLX-1
      G%FE(:,IX,:) = M%DX (:,IX+1,:) / (M%DX (:,IX+1,:) + M%DX (:,IX,:))
   ENDDO
   !$OMP END PARALLEL DO 


   G%FW(:,1,:) = 0D0
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX) SHARED(G,M,NCELLX)
   DO IX = 2, NCELLX
      G%FW(:,IX,:) = M%DX (:,IX-1,:) / (M%DX (:,IX-1,:) + M%DX (:,IX,:))
   ENDDO
   !$OMP END PARALLEL DO 

ENDIF

!y-direction:
! As in 2D/3D ther is no deformation only calculate at the first time step
IF ((G%NTIMESTEPS .EQ. 1) .AND. (NCELLY .GT. 1 )) THEN
   G%FN(:,:,NCELLY) = 0D0
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IY) SHARED(G,M,NCELLY)
   DO IY = 1, NCELLY-1
      G%FN(:,:,IY) = M%DY (:,:,IY+1) / (M%DY (:,:,IY+1) + M%DY (:,:,IY))
   ENDDO
   !$OMP END PARALLEL DO 

   G%FS(:,:,1) = 0D0
   !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IX) SHARED(G,M,NCELLX)
   DO IY = 2, NCELLY
      G%FS(:,:,IY) = M%DY (:,:,IY-1) / (M%DY (:,:,IY-1) + M%DY (:,:,IY))
   ENDDO
   !$OMP END PARALLEL DO 
ENDIF

! *****************************************************************************
END SUBROUTINE CALCULATE_INTERFACE_MESH_FACTORS
! *****************************************************************************

!*******************************************************************************
SUBROUTINE INIT_MESH(IMESH)
!*******************************************************************************
INTEGER, INTENT(IN) :: IMESH
TYPE(GPYRO_MESH_TYPE), POINTER :: M

G=>GPM(IMESH)
M => G%MESH

M%DIMENSION = 0
IF (M%NCELLZ .GT.1 ) M%DIMENSION = M%DIMENSION +1
IF (M%NCELLX .GT.1 ) M%DIMENSION = M%DIMENSION +1
IF (M%NCELLY .GT.1 ) M%DIMENSION = M%DIMENSION +1


G%THICKNESS = M%ZDIM


IF (M%DIMENSION .EQ. 0) THEN 
   CALL INIT_MESH_0D(IMESH)
   RETURN
ENDIF

CALL GENERATE_REGULAR_CARTESIAN_GRID(IMESH)
CALL  APPLY_OBST_TO_MESH(IMESH)
CALL CALCULATE_DIFFUSIVE_MESH_FACTORS()

!*******************************************************************************
END SUBROUTINE INIT_MESH
!*******************************************************************************

!*******************************************************************************
SUBROUTINE INIT_MESH_0D(IMESH)
!*******************************************************************************
INTEGER, INTENT(IN) :: IMESH
TYPE(GPYRO_MESH_TYPE), POINTER :: M
REAL(EB) :: dv

G=>GPM(IMESH)
M => G%MESH


!dv=1D0/G%RP(1,1,1)
dv=1D0
M%IMASK(:,:,:) = .FALSE.
M%ID_OBST(:,:,:) = 0
M%DV   (:,:,:)  = dv
M%DX   (:,:,:) = 1D0
M%DY   (:,:,:) = 1D0
M%DXDY (:,:,:) = 1D0

M%DZN  (:,:,:) = dv
M%DZ   (:,:,:) = dv
M%DZT  (:,:,:) = dv
M%DZB  (:,:,:) = dv
M%Z    (:,:,:) = 1D0

!*******************************************************************************
END SUBROUTINE INIT_MESH_0D
!*******************************************************************************

!*******************************************************************************
SUBROUTINE GENERATE_REGULAR_CARTESIAN_GRID(IMESH)
!*******************************************************************************
! Initializes raw grid geometric quantities 
!-------------------------------------------------------------------------------

INTEGER, INTENT(IN) :: IMESH
TYPE(GPYRO_MESH_TYPE), POINTER :: M
INTEGER :: NCELLZ, NCELLX, NCELLY
REAL(EB) :: ZDIM, XDIM, YDIM
INTEGER :: IZ, IX, IY
REAL(EB) :: dx,dy,dz

G=>GPM(IMESH)
M => G%MESH

NCELLZ = M%NCELLZ
NCELLX = M%NCELLX
NCELLY = M%NCELLY

ZDIM = M%ZDIM
XDIM = M%XDIM
YDIM = M%YDIM

G%THICKNESS = ZDIM

IF (M%HALF_CELLS_AT_BC) THEN
   !--------------------------------------------------------------------------
   ! Case 1: Half cells at the boundaries
   !--------------------------------------------------------------------------
   dz=ZDIM / REAL(NCELLZ - 1, EB)
   
   ! Uniform spacing for internal z-thicknesses
   M%DZN(:,:,:) = dz
   M%DZT(:,:,:) = dz
   M%DZB(:,:,:) = dz
   ! Adjust half-thickness on the top and bottom boundaries
   M%DZN(1     ,:,:)= 0.5*dz
   M%DZN(NCELLZ,:,:)= 0.5*dz

   M%DZ (:,:,:) = M%DZN(:,:,:)

ELSE 
   !--------------------------------------------------------------------------
   ! Case 2: Full cells at edges + ghost cells (experimental mode)
   ! NZ=1 and NZ=NCELLZ are ghost cells; internal grid uses NCELLZ-2 real cells
   !--------------------------------------------------------------------------
   dz=ZDIM / REAL(NCELLZ - 2, EB)
   M%DZN(:,:,:) = dz
   M%DZ (:,:,:) = dz
   M%DZT(:,:,:) = dz
   M%DZB(:,:,:) = dz
ENDIF

CALL COMPUTE_Z_CELL_POSITION(IMESH)

!Set x-direction spacing:
IF ((NCELLX .GT. 1) .AND. M%HALF_CELLS_AT_BC) THEN
   dx=XDIM/REAL(NCELLX-1,EB)
   M%DX (:,:,:) = dx

   ! Adjust half-thickness on boundaries
   M%DX (:,1     ,:) = 0.5D0 * dx
   M%DX (:,NCELLX,:) = 0.5D0 * dx

   !Delta-x east and west:
   M%DXE(:,:,:) = dx
   M%DXW(:,:,:) = dx

   ! Set x of each cell:
   M%X(1) = 0D0
   M%X(2) = M%X(1) + M%DX (1,1,1) + 0.5D0*M%DX (1,2,1)
   DO IX = 3, NCELLX-1
      M%X(IX) = M%X(IX-1) + dx
   ENDDO
   M%X(NCELLX) = M%X(NCELLX-1) + 0.5D0 * M%DX (1,NCELLX-1,1) + M%DX (1,NCELLX,1)

ELSEIF ((NCELLX .GT. 1) .AND. (.NOT. M%HALF_CELLS_AT_BC)) THEN

   dx=XDIM/REAL(NCELLX-2,EB) 
   M%DX (:,:,:) = dx

   !Delta-x east and west:
   M%DXE(:,:,:) = dx
   M%DXW(:,:,:) = dx

   ! x-position of the center of the cells
   M%X(1) = -0.5D0 *M%DX (1,1,1) ! Center of the gost cell (not used)
   DO IX = 2, NCELLX
      M%X(IX) = M%X(IX-1) + dx
   ENDDO

ELSE ! NCELLX=1
   M%DX (:,:,:) = 1D0
ENDIF

!Set y-direction spacing:
IF ((NCELLY .GT. 1) .AND. M%HALF_CELLS_AT_BC)  THEN
   dy=YDIM/REAL(NCELLY-1,EB)
   M%DY (:,:,:) = dy
   M%DY (:,:,1    ) = 0.5D0 * dy
   M%DY(:,:,NCELLY) = 0.5D0 * dy

   !Delta-y north and south:
   M%DYN(:,:,:) = dy
   M%DYS(:,:,:) = dy 

   ! Set y of each cell:
   M%Y(1) = 0D0
   M%Y(2) = M%Y(1) + M%DY(1,1,1) + 0.5D0*M%DY(1,1,2)
   DO IY = 3, NCELLY-1
      M%Y(IY) = M%Y(IY-1) + dy
   ENDDO
   M%Y(NCELLY) = M%Y(NCELLY-1) + 0.5D0 * M%DY(1,1,NCELLY-1) + M%DY(1,1,NCELLY)

ELSEIF ((NCELLX .GT. 1) .AND. (.NOT. M%HALF_CELLS_AT_BC)) THEN
   dy = YDIM/REAL(NCELLY-2,EB)
   M%DY  (:,:,:) = dy

   !Delta-y north and south:
   M%DYN(:,:,:) = dy
   M%DYS(:,:,:) = dy

   ! y-position of the center of the cells
   M%Y(1) = -0.5D0 *M%DY (1,1,1) ! Center of the gost cell (not used)
   DO IY = 2, NCELLY
      M%Y(IY) = M%Y(IY-1) + dy
   ENDDO

ELSE !NCELLY=1
   M%DY (:,:,:) = 1D0
ENDIF

! Set Volume and surfaces
!The mesh is orthogonal, the surface and volume are straightforward to compute.

!$OMP PARALLEL DO COLLAPSE(2)
DO IX = 1, NCELLX
DO IY = 1, NCELLY
   !$OMP SIMD
   DO IZ = 1, NCELLZ
      dx = M%DX(IZ,IX,IY)
      dy = M%DY(IZ,IX,IY)
      dz = M%DZ(IZ,IX,IY)

      M%DXDY(IZ,IX,IY) = dx * dy
      M%SYN (IZ,IX,IY) = dx * dz
      M%SYS (IZ,IX,IY) = dx * dz
      M%SXE (IZ,IX,IY) = dy * dz
      M%SXW (IZ,IX,IY) = dy * dz
      M%DV  (IZ,IX,IY) = dx * dy * dz
   END DO
END DO
END DO
!$OMP END PARALLEL DO

!******************************************************************************
END SUBROUTINE GENERATE_REGULAR_CARTESIAN_GRID
!******************************************************************************


!*******************************************************************************
SUBROUTINE COMPUTE_Z_CELL_POSITION(IMESH)
!*******************************************************************************
! Computes the z-coordinate of each cell center based on the cell thicknesses.
!-------------------------------------------------------------------------------
INTEGER, INTENT(IN)::IMESH
INTEGER :: NCELLZ, NCELLX, NCELLY
INTEGER :: IZ, IX, IY
TYPE (GPYRO_MESH_TYPE), POINTER :: M

G => GPM(IMESH)
M => G%MESH

NCELLZ = M%NCELLZ
NCELLX = M%NCELLX
NCELLY = M%NCELLY


IF (M%HALF_CELLS_AT_BC) THEN
   !--------------------------------------------------------------------------
   ! Case 1: Half cells at boundary conditions
   !--------------------------------------------------------------------------

   ! The first z-position is set to zero (reference plane)

   M%Z(1,:,:) = 0.0D0

   ! Compute z of the second layer
   !$OMP PARALLEL DO PRIVATE(IX,IY) DEFAULT(SHARED) COLLAPSE(2)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
      M%Z(2,IX,IY) = M%Z(1,IX,IY) + M%DZN(1,IX,IY) + 0.5D0 * M%DZN(2,IX,IY)
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Compute intermediate z-positions (from 3 to NCELLZ-1)
   !$OMP PARALLEL DO PRIVATE(IX,IY,IZ) DEFAULT(SHARED) COLLAPSE(2)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 3, NCELLZ-1
      M%Z(IZ,IX,IY) = M%Z(IZ-1,IX,IY) + 0.5D0*( M%DZN(IZ-1,IX,IY) + M%DZN(IZ,IX,IY) )
   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Last top cells
   !$OMP PARALLEL DO PRIVATE(IX,IY) DEFAULT(SHARED)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
      M%Z(NCELLZ,IX,IY) = M%Z(NCELLZ-1,IX,IY) + 0.5D0*M%DZN(NCELLZ-1,IX,IY) + M%DZN(NCELLZ,IX,IY)
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

ELSE
   !--------------------------------------------------------------------------
   ! Case 2: No half-cells and ghost cell (beta feature)
   !--------------------------------------------------------------------------
   !$OMP PARALLEL DO PRIVATE(IX,IY,IZ) DEFAULT(SHARED) COLLAPSE(2)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
      M%Z(1,IX,IY) = -0.5D0 * M%DZN(1,IX,IY)   ! ghost cell center
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   !$OMP PARALLEL DO PRIVATE(IX,IY,IZ) DEFAULT(SHARED) COLLAPSE(2)
   DO IY = 1, NCELLY
   DO IX = 1, NCELLX
   DO IZ = 2, NCELLZ
       M%Z(IZ,IX,IY) = M%Z(IZ-1,IX,IY) + 0.5D0*(M%DZN(IZ-1,IX,IY) + M%DZN(IZ,IX,IY))
   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDIF

END SUBROUTINE COMPUTE_Z_CELL_POSITION
!*******************************************************************************


!******************************************************************************
SUBROUTINE APPLY_OBST_TO_MESH(IMESH)
!******************************************************************************
INTEGER, INTENT(IN) :: IMESH
TYPE(GPYRO_MESH_TYPE), POINTER :: M
INTEGER :: IOBST
LOGICAL :: VALID_OBSTS

G=>GPM(IMESH)
M=>G%MESH

!========================================================!
!=========== SET DEFAULTS FOR ENTIRE DOMAIN =============!
!========================================================!

IOBST = 0 ! Correspond to default mesh configuration
M%ID_OBST(:,:,:) = IOBST

!Begin code for masking:
M%IMASK(:,:,:) = .TRUE.
!============ SET GEOMETRY FROM OBST ====================!

! Read info from &GPYRO_GEOM. Note that the OBST's specified in the input deck are read in
! first, and then can be "overwritten" with geometry from the geometry file. 

VALID_OBSTS = .FALSE.
DO IOBST = 1, GPG%NOBST 
   IF (GPG%GEOM(IOBST)%IMESH .NE. IMESH) CYCLE
   VALID_OBSTS = .TRUE.
   CALL SETUP_MASKING(IMESH,IOBST)
ENDDO

! READ GEOMTRY FILE IF PRESENT (Don't work)
CALL READ_GEOMETTRY_FILE(IMESH)

! If no geometry info present, assume the entire domain extents are unmasked:
IF ((.NOT. G%ORIENTATION_FILE_EXISTS) .AND. (.NOT. VALID_OBSTS)) THEN
   M%IMASK(:,:,:) = .FALSE.
ENDIF

!******************************************************************************
END SUBROUTINE APPLY_OBST_TO_MESH 
!******************************************************************************


!******************************************************************************
SUBROUTINE READ_GEOMETTRY_FILE(IMESH)
!******************************************************************************
! (Beta Feature (don't work))
INTEGER, INTENT(IN) :: IMESH
TYPE(GPYRO_MESH_TYPE), POINTER :: M
REAL(EB) :: XB(6)
INTEGER :: ICNUM
INTEGER :: SURF_IDX(1:6)
INTEGER :: NCELLZ, NCELLX, NCELLY
LOGICAL :: GEOMETRY_FILE_EXISTS
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
M=>G%MESH
NCELLZ = M%NCELLZ !Number of cells in z direction
NCELLX = M%NCELLX !Number of cells in x direction
NCELLY = M%NCELLY !Number of cells in y direction


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
! *****************************************************************************
END SUBROUTINE READ_GEOMETTRY_FILE
! *****************************************************************************


! *****************************************************************************
SUBROUTINE SETUP_MASKING(IMESH,ID_OBST)
! *****************************************************************************
INTEGER , INTENT(IN ) :: ID_OBST 
INTEGER , INTENT(IN ) :: IMESH 
TYPE(GPYRO_MESH_TYPE), POINTER :: M
REAL(EB), ALLOCATABLE :: COORDX1(:), COORDX2(:)
REAL(EB), ALLOCATABLE :: COORDY1(:), COORDY2(:)
REAL(EB), ALLOCATABLE :: COORDZ1(:), COORDZ2(:)
INTEGER :: NCELLX,NCELLY,NCELLZ
REAL(EB) :: XB1, XB2, XB3, XB4, XB5, XB6
REAL(EB) :: TARGX1, TARGX2, TARGY1, TARGY2, TARGZ1, TARGZ2
INTEGER :: IX1,IX2,IY1,IY2,IZ1,IZ2

G => GPM(IMESH)
M => G%MESH

NCELLZ = M%NCELLZ
NCELLX = M%NCELLX
NCELLY = M%NCELLY

ALLOCATE(COORDX1(NCELLX), COORDX2(NCELLX))
ALLOCATE(COORDY1(NCELLY), COORDY2(NCELLY))
ALLOCATE(COORDZ1(NCELLZ), COORDZ2(NCELLZ))

XB1 = GPG%GEOM(ID_OBST)%X1
XB2 = GPG%GEOM(ID_OBST)%X2
XB3 = GPG%GEOM(ID_OBST)%Y1
XB4 = GPG%GEOM(ID_OBST)%Y2
XB5 = GPG%GEOM(ID_OBST)%Z1
XB6 = GPG%GEOM(ID_OBST)%Z2


IZ1 = 1
IZ2 = 1
IX1 = 1
IX2 = 1
IY1 = 1
IY2 = 1


IF (NCELLX .GT. 1) THEN 
   COORDX1(:) = M%X(:) - 0.5*M%DX (1,:,1)
   TARGX1     = XB1 
   IX1        = IJK_FROM_XYZ(COORDX1, NCELLX, TARGX1, 1)

   COORDX2(:) = M%X(:) + 0.5*M%DX (1,:,1)
   TARGX2     = XB2 
   IX2        = IJK_FROM_XYZ(COORDX2, NCELLX, TARGX2, 2)
ENDIF

IF (NCELLY .GT. 1) THEN
   COORDY1(:) = M%Y(:) - 0.5*M%DY(1,1,:)
   TARGY1     = XB3 
   IY1        = IJK_FROM_XYZ(COORDY1, NCELLY, TARGY1, 1)

   COORDY2(:) = M%Y(:) + 0.5*M%DY(1,1,:)
   TARGY2     = XB4 
   IY2        = IJK_FROM_XYZ(COORDY2, NCELLY, TARGY2, 2)
ENDIF

IF (GPG%GEOMETRY_IS_UPSIDE_DOWN) THEN
   COORDZ1(:) = M%Z(:,1,1) - 0.5*M%DZN(:,1,1)
   TARGZ1     = M%ZDIM - XB6 
   IZ1        = IJK_FROM_XYZ(COORDZ1, NCELLZ, TARGZ1, 1)

   COORDZ2(:) = M%Z(:,1,1) + 0.5*M%DZN(:,1,1)
   TARGZ2     = M%ZDIM - XB5 
   IZ2        = IJK_FROM_XYZ(COORDZ2, NCELLZ, TARGZ2, 2)
ELSE
   COORDZ1(:) = M%Z(:,1,1) - 0.5*M%DZN(:,1,1)
   TARGZ1     = XB5 
   IZ1        = IJK_FROM_XYZ(COORDZ1, NCELLZ, TARGZ1, 1)

   COORDZ2(:) = M%Z(:,1,1) + 0.5*M%DZN(:,1,1)
   TARGZ2     = XB6
   IZ2        = IJK_FROM_XYZ(COORDZ2, NCELLZ, TARGZ2, 2)
ENDIF

M%IMASK   (IZ1:IZ2,IX1:IX2,IY1:IY2) = .FALSE. ! Cells not masked

M%ID_OBST(IZ1:IZ2,IX1:IX2,IY1:IY2) = ID_OBST

! *****************************************************************************
END SUBROUTINE SETUP_MASKING
! *****************************************************************************

! *****************************************************************************
INTEGER FUNCTION IJK_FROM_XYZ(COORD,NCELL,TARG,ITYPE)
! *****************************************************************************

INTEGER, INTENT(IN) :: NCELL, ITYPE
REAL(EB), INTENT(IN) :: TARG
REAL(EB), DIMENSION(:), INTENT(IN) :: COORD
REAL(EB) :: DIFF, MINDIFF
INTEGER :: I

IF (NCELL .EQ. 1) THEN
   IJK_FROM_XYZ = 1
ELSE
   MINDIFF = 9D9
   IF (ITYPE .EQ. 1) THEN
      DO I = 1, NCELL
         DIFF = ABS(COORD(I) - TARG)
         IF (DIFF .LT. MINDIFF) THEN
            MINDIFF = DIFF
            IJK_FROM_XYZ = I
         ENDIF
      ENDDO
   ELSE
      DO I = 1, NCELL
         DIFF = ABS(COORD(I) - TARG)
         IF (DIFF .LT. MINDIFF) THEN
            MINDIFF = DIFF
            IJK_FROM_XYZ = I
         ENDIF
      ENDDO
   ENDIF
ENDIF

! *****************************************************************************
END FUNCTION IJK_FROM_XYZ
! *****************************************************************************

! *****************************************************************************
END MODULE GPYRO_MESH
! *****************************************************************************
