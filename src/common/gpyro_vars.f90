! *****************************************************************************
MODULE GPYRO_VARS
! *****************************************************************************

USE PREC !Following FDS, single precison is FB and double precision is EB

! For system timing:
INTEGER(8) :: CLOCK_COUNT_RATE, CLOCK_COUNT_MAX 

INTEGER :: IGPYRO_TYPE ! =1 for standalone, =2 for FDS, =3 for GA

! These are the logical units used for i/o. 
INTEGER, PARAMETER :: LUPOINT   = 200, &
                      LUSUM     = 200, & !No longer used
                      LUPROF    = 250, &
                      LUTIME    = 432, &
                      LUINPUT   = 801, &
                      LUOUTPUT  = 802, &
                      LUDOTOUT  = 901, &
                      LUSMV     = 9000, &
                      LUCSV     = 9500
                      
! Misc.
REAL(EB), PARAMETER :: PI = 3.14159265358979323846D0

REAL(EB), PARAMETER :: R = 8.314462 ! Ideal gas constant (J⋅K−1⋅mol−1)
      
REAL(EB), PARAMETER :: SIGMA=5.670373E-8_EB

! "alpha" coefficients for omega function (CWL MS thesis)

REAL(EB), PARAMETER :: AD(6) = [ &
   1.080794_EB, -0.16033_EB, 0.605009_EB, &
  -0.88524_EB,  2.115672_EB, -2.98308_EB ]


! FDS related
REAL(EB), DIMENSION(1:6) :: IGNLOC
REAL(EB) :: IGNPOW,TIGNSTART,TIGNSTOP

! Initial condition TYPE.  
TYPE :: INITIAL_CONDITION_TYPE
   REAL(EB), POINTER, DIMENSION(:) :: YI0
   REAL(EB), POINTER, DIMENSION(:) :: YJ0
   REAL(EB) :: TMP_INITIAL
   REAL(EB) :: TMPG_INITIAL
   REAL(EB) :: P_INITIAL
END TYPE 

! AllBC type: 
TYPE :: ALLBC_TYPE
   INTEGER :: SURF_IDX
   REAL(EB) :: T
   REAL(EB) :: QE
   REAL(EB) :: HC
   REAL(EB) :: NHC   
   REAL(EB) :: TINF
   LOGICAL  :: RERADIATION 
   REAL(EB) :: TFIXED
   REAL(EB) :: MDOTPP
   REAL(EB) :: PRES
   REAL(EB) :: QEG
   REAL(EB) :: HCG
   REAL(EB) :: TINFG
   REAL(EB) :: TFIXEDG
   REAL(EB) :: HM
   REAL(EB), POINTER, DIMENSION(:) :: YJINF
END TYPE

!Geometry TYPE:
TYPE :: GEOM_TYPE
   INTEGER  :: IMESH
   REAL(EB) :: Z1
   REAL(EB) :: Z2
   REAL(EB) :: X1
   REAL(EB) :: X2
   REAL(EB) :: Y1
   REAL(EB) :: Y2
   INTEGER :: ICNUM
   INTEGER :: SURF_IDX(1:6)
END TYPE

! Holds info about solid to solid or solid to solid plus gas reactions. 
! Read in from rxns worksheet
TYPE :: GPYRO_REACTION_TYPE
   CHARACTER(60) :: CFROM
   CHARACTER(60) :: CTO
   INTEGER       :: IFROM !Index from
   INTEGER       :: ITOS  !Index to (solid)
   REAL(EB)      :: Z
   REAL(EB)      :: E
   REAL(EB)      :: DHS
   REAL(EB)      :: DHV
   REAL(EB)      :: CHI
   REAL(EB)      :: ORDER
   REAL(EB)      :: ORDERO2
   INTEGER       :: IKINETICMODEL
   INTEGER       :: IO2TYPE
   REAL(EB)      :: M
   REAL(EB)      :: KCAT
   INTEGER       :: ICAT
   REAL(EB)      :: TCRIT
END TYPE
TYPE (GPYRO_REACTION_TYPE),    POINTER, DIMENSION (:) :: RXN
     
! Homogeneous gaseous reactions TYPE. Holds info about homogeneous gas reactions. 
! Read in from hgrxns worksheet. 
TYPE :: HG_REACTION_TYPE
   CHARACTER(60) :: CREACTANT1
   CHARACTER(60) :: CREACTANT2
   INTEGER       :: IREACTANT1
   INTEGER       :: IREACTANT2
   REAL(EB)      :: P
   REAL(EB)      :: Q
   REAL(EB)      :: B
   REAL(EB)      :: Z
   REAL(EB)      :: E
   REAL(EB)      :: DH
END TYPE
TYPE (HG_REACTION_TYPE), POINTER, DIMENSION (:) :: HGRXN 

! Solid property TYPE.  Holds, of course, solid properties. 
! This is read in from sprops worksheet.
TYPE :: SOLID_PROPERTY_TYPE
   CHARACTER(60), POINTER, DIMENSION (:) :: NAME
   REAL(EB), POINTER, DIMENSION (:)   :: K0Z
   REAL(EB), POINTER, DIMENSION (:)   :: NKZ
   REAL(EB), POINTER, DIMENSION (:)   :: R0
   REAL(EB), POINTER, DIMENSION (:)   :: NR
   REAL(EB), POINTER, DIMENSION (:)   :: C0
   REAL(EB), POINTER, DIMENSION (:)   :: NC
   REAL(EB), POINTER, DIMENSION (:)   :: EMIS
   REAL(EB), POINTER, DIMENSION (:)   :: KAPPA
   REAL(EB), POINTER, DIMENSION (:)   :: TMELT
   REAL(EB), POINTER, DIMENSION (:)   :: DHMELT
   REAL(EB), POINTER, DIMENSION (:)   :: SIGMA2MELT
   REAL(EB), POINTER, DIMENSION (:)   :: GAMMA
   REAL(EB), POINTER, DIMENSION (:)   :: PERMZ
   REAL(EB), POINTER, DIMENSION (:)   :: RS0
   REAL(EB), POINTER, DIMENSION (:)   :: PORE_DIAMETER
   REAL(EB), POINTER, DIMENSION (:)   :: K0X
   REAL(EB), POINTER, DIMENSION (:)   :: NKX
   REAL(EB), POINTER, DIMENSION (:)   :: PERMX
   REAL(EB), POINTER, DIMENSION (:)   :: K0Y
   REAL(EB), POINTER, DIMENSION (:)   :: NKY 
   REAL(EB), POINTER, DIMENSION (:)   :: PERMY
      INTEGER  :: NSSPEC
   INTEGER  :: NRXN

END TYPE
TYPE (SOLID_PROPERTY_TYPE) :: SPROP

! Gas property type.
! This is read in from gprops worksheet.
TYPE :: GAS_PROPERTY_TYPE
   INTEGER :: NGSPEC
   INTEGER :: NHGRXN
   INTEGER :: IBG
   INTEGER :: IO2
   REAL(EB) :: CPG 
   CHARACTER(60), POINTER, DIMENSION (:)   :: NAME
   REAL(EB), POINTER, DIMENSION(:) :: M
   REAL(EB), POINTER, DIMENSION(:) :: SIGMA
   REAL(EB), POINTER, DIMENSION(:) :: EPSOK
   REAL(EB), POINTER, DIMENSION(:) :: C0
   REAL(EB), POINTER, DIMENSION(:) :: NC
   REAL(EB), POINTER, DIMENSION(:,:) :: YIELDS !Heterogeneous
   REAL(EB), POINTER, DIMENSION(:,:) :: HGYIELDS !Homogeneous
END TYPE
TYPE (GAS_PROPERTY_TYPE) :: GPROP


TYPE :: WORK_PROFILE_DUMP
   ! For dump profils 
   REAL(EB), ALLOCATABLE :: BUFFER(:)
   CHARACTER(:), ALLOCATABLE :: LINE
END TYPE 


TYPE :: GPYRO_STORAGE_TYPE

  !DIM (NCELLZ, NCELLX, NCELLY)
  REAL(EB), ALLOCATABLE :: TP_W1(:,:,:), TP_W2(:,:,:)           ! Solid Temperature work array
  REAL(EB), ALLOCATABLE :: HP_W1(:,:,:), HP_W2(:,:,:)           ! Solid Hentalpy work array
  REAL(EB), ALLOCATABLE :: RP_W1(:,:,:), RP_W2(:,:,:)           ! Bulk Density work array
  REAL(EB), ALLOCATABLE :: DLTZ_W1(:,:,:), DLTZ_W2(:,:,:)       ! Z- mesh cell size work array
  REAL(EB), ALLOCATABLE :: RDLTZ_W1(:,:,:), RDLTZ_W2(:,:,:)     ! Mass (rho*dz) work array
  REAL(EB), ALLOCATABLE :: RSP_W1(:,:,:), RSP_W2(:,:,:)         ! Solid Density work array

  REAL(EB), ALLOCATABLE :: M_W1(:,:,:), M_W2(:,:,:)             ! Mean Moleculare weigth work array
  REAL(EB), ALLOCATABLE :: RG_W1(:,:,:), RG_W2(:,:,:)           ! Gas Density work array
  REAL(EB), ALLOCATABLE :: P_W1(:,:,:), P_W2(:,:,:)             ! Pressure work array
  REAL(EB), ALLOCATABLE :: TG_W1(:,:,:), TG_W2(:,:,:)           ! Gas Temperature work array
  REAL(EB), ALLOCATABLE :: HG_W1(:,:,:), HG_W2(:,:,:)           ! Gas Henthalpy work array
  REAL(EB), ALLOCATABLE :: POROSS_W1(:,:,:), POROSS_W2(:,:,:)   ! Porrosity work array

  !DIM (NSSPEC, NCELLZ, NCELLX, NCELLY)
  REAL(EB), ALLOCATABLE :: YI_W1(:,:,:,:), YI_W2(:,:,:,:)       ! Solid mass fractions work array
  REAL(EB), ALLOCATABLE :: XI_W1(:,:,:,:), XI_W2(:,:,:,:)       ! Volume fraction work array
  REAL(EB), ALLOCATABLE :: RYIDZP_W1(:,:,:,:), RYIDZP_W2(:,:,:,:)           ! rho*Yi*Dz rho*Yi*Dz work array
  REAL(EB), ALLOCATABLE :: RYIDZSIGMA_W1(:,:,:,:), RYIDZSIGMA_W2(:,:,:,:)   ! rho*Yi*Dz summation work array

!DIM (NGSPEC, NCELLZ, NCELLX, NCELLY)

  REAL(EB), ALLOCATABLE :: YJG_W1(:,:,:,:), YJG_W2(:,:,:,:)     ! Gas-phase mass fractions work array

END TYPE GPYRO_STORAGE_TYPE


TYPE :: SOLVER_CONVERGENCE_INFO

   LOGICAL, ALLOCATABLE, DIMENSION(:) :: CONVERGED_YIS  ! Convergence of solid species conservation
   LOGICAL, ALLOCATABLE, DIMENSION(:) :: CONVERGED_YJG  ! Convergence of gas species conservation
   LOGICAL :: CONVERGED_TMP   ! Convergence of solid enthalpy solver
   LOGICAL :: CONVERGED_HG    ! Convergence of gas enthalpy solver
   LOGICAL :: CONVERGED_P     ! Convergence of pressure solver
   LOGICAL :: CONVERGED_ALL   ! All solvers converged

   INTEGER , ALLOCATABLE, DIMENSION(:) :: ITER_YIS  ! Number of iterations to converge solid species
   INTEGER , ALLOCATABLE, DIMENSION(:) :: ITER_YJG  ! Number of iterations to converge gas species
   INTEGER :: ITER_TMP  ! Number of iterations to converge solid enthalpy and temperature
   INTEGER :: ITER_HG   ! Number of iterations to converge gas enthalpy and temperature
   INTEGER :: ITER_P    ! Number of iterations to converge pressure
   INTEGER :: ITER_ALL  ! Number of iterations to converge all solvers
   INTEGER :: ITER      ! Number of the actual subiteration

END TYPE


! Boundary condition type - holds all "calculated" quantities for each BC.
! In standalone and genetic algorithm implementations, there is only a single
! one-dimensional boundary condition, but when coupled to FDS there may be 
! hundreds or thousands of separate boundary conditions, and using this 
! "TYPE" is a convenient way to track each one. 
!
TYPE  :: GPYRO_MESH_TYPE

   REAL(EB) :: INITIAL_MASS

   REAL(EB) :: GX   !x-component of gravity vector
   REAL(EB) :: GY   !y-component of gravity vector
   REAL(EB) :: GZ   !z-component of gravity vector   

   REAL(EB) :: X0=0D0   !x absolute coordinate for FDS
   REAL(EB) :: Y0=0D0   !y absolute coordinate for FDS
   REAL(EB) :: Z0=0D0   !z absolute coordinate for FDS
   
   REAL(EB) :: ZDIM=0D0
   REAL(EB) :: XDIM=0D0
   REAL(EB) :: YDIM=0D0

   INTEGER :: NCELLZ
   INTEGER :: NCELLX
   INTEGER :: NCELLY
   INTEGER :: DIMENSION ! 0if 0D, 1 if 1D, 2 if 2D, 3 if 3D
   LOGICAL :: HALF_CELLS_AT_BC = .TRUE.

   REAL(EB), DIMENSION(1:3) :: MESH_CENTROIDS=0D0, MESH_EXTENTS=0D0, MESH_LL=0D0 !, MESH_G

   LOGICAL, POINTER, DIMENSION (:,:,:) :: NEEDSBCT, NEEDSBCB, NEEDSBCW, NEEDSBCE, NEEDSBCN, NEEDSBCS
   INTEGER, POINTER, DIMENSION (:,:,:) :: SURF_IDX_BCT
   INTEGER, POINTER, DIMENSION (:,:,:) :: SURF_IDX_BCB
   INTEGER, POINTER, DIMENSION (:,:,:) :: SURF_IDX_BCW
   INTEGER, POINTER, DIMENSION (:,:,:) :: SURF_IDX_BCE
   INTEGER, POINTER, DIMENSION (:,:,:) :: SURF_IDX_BCN
   INTEGER, POINTER, DIMENSION (:,:,:) :: SURF_IDX_BCS

   !REAL(EB), POINTER, DIMENSION (:,:,:) :: HCRT !contact resistance (top)
   !REAL(EB), POINTER, DIMENSION (:,:,:) :: HCRB !contact resistance (bottom)
   !REAL(EB), POINTER, DIMENSION (:,:,:) :: HCRW !contact resistance (west)
   !REAL(EB), POINTER, DIMENSION (:,:,:) :: HCRE !contact resistance (east)
   !REAL(EB), POINTER, DIMENSION (:,:,:) :: HCRS !contact resistance (south)
   !REAL(EB), POINTER, DIMENSION (:,:,:) :: HCRN !contact resistance (north)

   TYPE(GPYRO_STORAGE_TYPE) :: STORAGE

   REAL(EB), POINTER, DIMENSION (:,:,:) :: TP  !Solid T at point P
   REAL(EB), POINTER, DIMENSION (:,:,:) :: TPN !Solid T at point P (new)
   REAL(EB), POINTER, DIMENSION (:,:,:) :: HP  !Solid enthalpy at point P
   REAL(EB), POINTER, DIMENSION (:,:,:) :: HPN !Solid enthalpy at point P (new)
   REAL(EB), POINTER, DIMENSION (:,:,:) :: TG  !Gas T at point P
   REAL(EB), POINTER, DIMENSION (:,:,:) :: TGN !Gas T at point P (new)
   REAL(EB), POINTER, DIMENSION (:,:,:) :: HG  !Gas enthalpy at point P
   REAL(EB), POINTER, DIMENSION (:,:,:) :: HGN !Gas enthalpy at point P (new)
   REAL(EB), POINTER, DIMENSION (:,:,:) :: RG  !Gas density at cell P
   REAL(EB), POINTER, DIMENSION (:,:,:) :: RGN !Gas density at cell P (new)
   REAL(EB), POINTER, DIMENSION (:,:,:) :: P   !Pressure at point P
   REAL(EB), POINTER, DIMENSION (:,:,:) :: PN  !Pressure at point P (new)
   REAL(EB), POINTER, DIMENSION (:,:,:) :: M   !Mean Moleculare weigth at point P
   REAL(EB), POINTER, DIMENSION (:,:,:) :: MN  !Mean Moleculare weigth at point P(new)
   REAL(EB), POINTER, DIMENSION (:,:,:) :: RP  !Density at point P 
   REAL(EB), POINTER, DIMENSION (:,:,:) :: RPN !Density at point P (new)
   REAL(EB), POINTER, DIMENSION (:,:,:) :: RSP  !Solid density at point P 
   REAL(EB), POINTER, DIMENSION (:,:,:) :: RSPN !Solid density at point P (new)
   REAL(EB), POINTER, DIMENSION (:,:,:) :: RDLTZ  !Old rho*DZ
   REAL(EB), POINTER, DIMENSION (:,:,:) :: RDLTZN !New rho*DZ
   REAL(EB), POINTER, DIMENSION (:,:,:) :: POROSS   !Porosity at point P
   REAL(EB), POINTER, DIMENSION (:,:,:) :: POROSSN  !Porosity at point P (new)
   REAL(EB), POINTER, DIMENSION (:,:,:) :: HCV  !Volumetric heat transfer coefficient
   REAL(EB), POINTER, DIMENSION (:,:,:) :: RE  !Reynolds number
   REAL(EB), POINTER, DIMENSION (:,:,:) :: NU  !Nusselt number

   REAL(EB), POINTER, DIMENSION (:,:,:) :: CPS !Specific heat of the solid
            
   REAL(EB), POINTER, DIMENSION (:,:,:) :: SHP !Positive part of solid h source
   REAL(EB), POINTER, DIMENSION (:,:,:) :: SHM !Negative part of solid h source
   REAL(EB), POINTER, DIMENSION (:,:,:) :: SHNET !Net solid h source

   REAL(EB), POINTER, DIMENSION (:,:,:) :: SHGP !Positive part of gas h source
   REAL(EB), POINTER, DIMENSION (:,:,:) :: SHGM !Negative part of gas h source
            
   !REAL(EB), POINTER, DIMENSION (:,:,:)   :: RGPSIDZ  !gas density * porosity * grid size
   !REAL(EB), POINTER, DIMENSION (:,:,:)   :: RGPSIDZN !Same, but new

   REAL(EB), POINTER, DIMENSION (:,:,:)   :: D12 !Diffusion coefficient

   LOGICAL, POINTER, DIMENSION (:,:,:)   :: IMASK !Geometry mask (True for blocked off, False for not blocked off)
   

   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: YI  !Condensed-species mass fractions
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: YIN !Same, but "new"
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: HI  !Enthalpy of solid species i (hi)
         
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: RRT !Temperature part of reaction rate
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: RRY !Mass fraction part of reaction rate
         
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: YJG  !Gas-phase mass fractions
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: YJGN !Same, but "new"

   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: RYIDZ0  !rho*Yi*Dz (at t = 0)
         
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: RYIDZP  !rho*Yi*Dz at point P
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: RYIDZPN !same, but "new"
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: RYIDZSIGMA  !rho*Yi*Dz summation
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: RYIDZSIGMAN !rho*Yi*Dz0 summation 

   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: UNREACTEDNESS ! 1 minus conversion
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: XI  !Volume fraction of solid species i
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: XIN !Same, but next 
         
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: OMSOLIDFRAC !One minus solid fraction
   !REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: SOLIDFRAC   !Solid fraction
         
   ! OMEGASFBK is Omega (solid) - formation   of condensed species B from reaction k
   ! OMEGASDAK is Omega (solid) - destruction of condensed species A from reaction k
   ! OMEGASFGK is Omega (solid) - formation   of all gases           from reaction k
   ! OMEGASFG  is Omega (solid) - total formation rate of all gases (summed over all k)
   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: OMEGASFBK
   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: OMEGASDAK
   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: OMEGASFGK
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: OMEGASFG

   ! OMEGASFJK is Omega (solid) - formation   of gaseous species J from reaction k
   ! OMEGASDJK is Omega (solid) - destruction of gaseous species J from reaction k
   REAL(EB), POINTER, DIMENSION (:,:,:,:,:) :: OMEGASFJK
   REAL(EB), POINTER, DIMENSION (:,:,:,:,:) :: OMEGASDJK
         
   ! SOMEGA stores the total formaton/destruction rates of solid (condensed) 
   ! species from condensed-phase reactions
   ! GOMEGA stores the total formation/destruction rates of gas-phase
   ! species from condensed-phase reaction
   REAL(EB), POINTER, DIMENSION (:,:,:,:,:) :: SOMEGA
   REAL(EB), POINTER, DIMENSION (:,:,:,:,:) :: GOMEGA

   LOGICAL, POINTER, DIMENSION(:,:,:,:) :: IS_REACTING ! Is a reaction currently happening in this cell?
   LOGICAL :: DEFORMATION ! allow deformation ?
   ! HGRR is homogeneous gaseous reaction rate
   ! OMEGAGDJL is Omega (gaseous) - destruction of gaseous species J from reaction l (L)	   
   ! HGOMEGA stores the formation/destruction of gaseous species from 
   ! homogeneous gaseous reactions
   REAL(EB), POINTER, DIMENSION (:,:,:,:)   :: HGRR 
   REAL(EB), POINTER, DIMENSION (:,:,:,:,:) :: OMEGAGDJL
   REAL(EB), POINTER, DIMENSION (:,:,:,:,:) :: HGOMEGA

   REAL(EB), POINTER, DIMENSION (:,:,:)     :: QSG !Heat transfer from solid to gas
   LOGICAL , POINTER, DIMENSION (:,:,:)     :: CONSUMED !Is a cell completely consumed?

   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: MDOTPPZ !Mass flux of each species in the z-direction
                 
   REAL(EB), POINTER, DIMENSION(:,:) :: THICKNESS    !Thickness

   !New section for 3D solvers:
   REAL(EB), POINTER, DIMENSION (:) :: Z ! Z coordinate
   REAL(EB), POINTER, DIMENSION (:) :: X ! X coordinate
   REAL(EB), POINTER, DIMENSION (:) :: Y ! Y coordinate

   REAL(EB), POINTER, DIMENSION (:,:,:) :: FT 
   REAL(EB), POINTER, DIMENSION (:,:,:) :: FB 
   REAL(EB), POINTER, DIMENSION (:,:,:) :: FE 
   REAL(EB), POINTER, DIMENSION (:,:,:) :: FW 
   REAL(EB), POINTER, DIMENSION (:,:,:) :: FN 
   REAL(EB), POINTER, DIMENSION (:,:,:) :: FS 

   REAL(EB), POINTER, DIMENSION (:,:,:) :: DLTZ  ! DLTZ  - Delta z 
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DLTX  ! DLTX  - Delta x 
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DLTY  ! DLTY  - Delta y 
   
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DLTZN ! DLTZN - Delta z (new), i.e. at t + Dt
   !REAL(EB), POINTER, DIMENSION (:,:,:) :: DLTXN ! DLTXN - Delta x (new), i.e. at t + Dt
   !REAL(EB), POINTER, DIMENSION (:,:,:) :: DLTYN ! DLTYN - Delta y (new), i.e. at t + Dt

   REAL(EB), POINTER, DIMENSION (:,:,:) :: DZT  !Dz top
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DZB  !Dz bottom
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DXW  !Dx west
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DXE  !Dx east
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DYN  !Dy north
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DYS  !Dy south

   REAL(EB), POINTER, DIMENSION (:,:,:) :: DV     !Dx*Dy*Dz (1D->m, 2D->m², 3D -> m³)
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DXDY   !Dx*Dy
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DXDZ   !Dx*Dz
   REAL(EB), POINTER, DIMENSION (:,:,:) :: DYDZ   !Dy*Dz

   REAL(EB), POINTER, DIMENSION (:,:,:) :: KZ   !Effective Thermal conductivity in z direction
   REAL(EB), POINTER, DIMENSION (:,:,:) :: KX   !Effective Thermal conductivity in x direction
   REAL(EB), POINTER, DIMENSION (:,:,:) :: KY   !Effective Thermal conductivity in y direction

   REAL(EB), POINTER, DIMENSION (:,:,:) :: KOCT !(k/c) at top interface
   REAL(EB), POINTER, DIMENSION (:,:,:) :: KOCB !(k/c) at bottom interface
   REAL(EB), POINTER, DIMENSION (:,:,:) :: KOCE !(k/c) at east interface 
   REAL(EB), POINTER, DIMENSION (:,:,:) :: KOCW !(k/c) at west interface
   REAL(EB), POINTER, DIMENSION (:,:,:) :: KOCN !(k/c) at north interface 
   REAL(EB), POINTER, DIMENSION (:,:,:) :: KOCS !(k/c) at south interface

   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: KOCHIT !(k/c)*hi at top interface
   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: KOCHIB !(k/c)*hi at bottom interface
   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: KOCHIE !(k/c)*hi at east interface
   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: KOCHIW !(k/c)*hi at west interface
   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: KOCHIN !(k/c)*hi at north interface
   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: KOCHIS !(k/c)*hi at south interface

   REAL(EB), POINTER, DIMENSION (:,:,:) :: KAPPA ! Absorption coefficient

   REAL(EB), POINTER, DIMENSION (:,:,:) :: SHYI

   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PERMZ    !Permeability in z-direction at point P
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PERMX    !Permeability in x-direction at point P
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PERMY    !Permeability in y-direction at point P

   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PERMONUT !Permeability/viscosity at top interface
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PERMONUB !Permeability/viscosity at bottom interface
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PERMONUE !Permeability/viscosity at east interface
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PERMONUW !Permeability/viscosity at west interface
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PERMONUN !Permeability/viscosity at north interface
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PERMONUS !Permeability/viscosity at south interface

   INTEGER , ALLOCATABLE, DIMENSION (:,:,:)   :: INDICE_OF_BC_GAS_OUT ! Indice of the Bondary cell where the gas get out


   REAL(EB), POINTER, DIMENSION(:,:,:) :: MDOTPPDARCYT !z-direction mass flux (top)
   REAL(EB), POINTER, DIMENSION(:,:,:) :: MDOTPPDARCYB !z-direction mass flux (bottom)
   REAL(EB), POINTER, DIMENSION(:,:,:) :: MDOTPPDARCYE !x-direction mass flux (east)
   REAL(EB), POINTER, DIMENSION(:,:,:) :: MDOTPPDARCYW !x-direction mass flux (west)
   REAL(EB), POINTER, DIMENSION(:,:,:) :: MDOTPPDARCYN !y-direction mass flux (north)
   REAL(EB), POINTER, DIMENSION(:,:,:) :: MDOTPPDARCYS !y-direction mass flux (south)

   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PSIRGDT !porosity * gas density * D (top)
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PSIRGDB !porosity * gas density * D (bottom)
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PSIRGDE !porosity * gas density * D (east)
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PSIRGDW !porosity * gas density * D (west)
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PSIRGDN !porosity * gas density * D (north)
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: PSIRGDS !porosity * gas density * D (south)

   REAL(EB), POINTER, DIMENSION (:,:,:)   :: BUOYANCYZ
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: BUOYANCYX
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: BUOYANCYY

   REAL(EB), POINTER, DIMENSION (:,:,:)   :: RGNT
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: RGNB
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: RGNE
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: RGNW
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: RGNN
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: RGNS

   !REAL(EB), POINTER, DIMENSION(:,:,:) :: GPYRO_CONV_FLUX,GPYRO_DIFF_FLUX,GPYRO_MASS_FLUX
   
   REAL(EB), POINTER, DIMENSION (:,:,:)   :: RESIDUAL_TMP, RESIDUAL_P, RESIDUAL_HG 
   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: RESIDUAL_YIS, RESIDUAL_YJG
   
  ! REAL(EB) :: TLAST !Last time at which GPYRO_PYROLYSIS was called

   INTEGER  :: NTIMESTEPS  ! Number of Actual time step iteration
   
   !Last Time Caled
   REAL(EB) :: TLAST_TMP
   REAL(EB) :: TLAST_YIS
   REAL(EB) :: TLAST_YJG
   REAL(EB) :: TLAST_P
   REAL(EB) :: TLAST_HG
   !Times steps
   REAL(EB) :: DTIME_TMP
   REAL(EB) :: DTIME_YIS
   REAL(EB) :: DTIME_YJG
   REAL(EB) :: DTIME_P
   REAL(EB) :: DTIME_HG
   
   TYPE (GPYRO_BC_TYPE), POINTER, DIMENSION(:,:,:,:) :: GPYRO_BOUNDARY_CONDITION

   ! Variables for orientation: 
   LOGICAL :: ORIENTATION_FILE_EXISTS
   CHARACTER(200) :: ORIENTATION_FILE
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: ORI, K_TENSOR

   !WORK ARRAYS:
   REAL(EB), POINTER, DIMENSION(:,:,:) :: RWORK01, RWORK02, RWORK03, RWORK04, RWORK05, RWORK06, RWORK07, RWORK08, RWORK09, &
                                          RWORK10, RWORK11, RWORK12, RWORK13, RWORK14, RWORK15, RWORK16, RWORK17, RWORK18, &
                                          RWORK19, RWORK20, RWORK21, RWORK22, RWORK23, RWORK24, RWORK25, RWORK26, RWORK27, &
                                          RWORK28, RWORK29, RWORK30, RWORK31, RWORK32, RWORK33, RWORK34, RWORK35, RWORK36, &
                                          RWORK37, RWORK38, RWORK39, RWORK40, RWORK41, RWORK42, RWORK43, RWORK44, RWORK45, & 
                                          RWORK46, RWORK47, RWORK48, RWORK49, RWORK50, RWORK51, RWORK52

   REAL(EB), POINTER, DIMENSION (:,:,:,:) :: RWORK100, RWORK101, RWORK102, RWORK103, RWORK104, RWORK105, RWORK106, RWORK110, &
                                             RWORK111, RWORK112, RWORK120, RWORK121

   LOGICAL, POINTER, DIMENSION(:,:,:) :: LWORK01

   TYPE(WORK_PROFILE_DUMP), ALLOCATABLE :: WORK_DUMP(:)

   TYPE (SOLVER_CONVERGENCE_INFO) :: CONV_INFO


   LOGICAL :: SMOKEVIEW_FILE_OPENED_ALREADY = .FALSE. 

END TYPE



TYPE (GPYRO_MESH_TYPE), ALLOCATABLE, TARGET, DIMENSION(:) :: GPM
TYPE (GPYRO_MESH_TYPE), POINTER :: G

TYPE :: GPYRO_BC_TYPE
   REAL(EB) :: QE
   REAL(EB) :: QENET  
   REAL(EB) :: HC0
   REAL(EB) :: NHC   
   REAL(EB) :: TINF
   REAL(EB) :: TFIXED
   REAL(EB) :: HFIXED
   REAL(EB) :: PRES
   REAL(EB) :: MFLUX
   REAL(EB) :: HM0
   REAL(EB) :: QEG
   REAL(EB) :: HC0G
   REAL(EB) :: TINFG
   REAL(EB) :: TFIXEDG
   REAL(EB) :: HFIXEDG
   REAL(EB), POINTER, DIMENSION(:) :: YJINF
   LOGICAL  :: RERAD
   ! The remaining variables are for fds coupling:
   REAL(EB) :: EMISSIVITY
   REAL(EB) :: T_SURFACE
   REAL(EB) :: T_SURFACE_OLD
   REAL(EB) :: QRADOUT
   REAL(EB) :: QCONF
   REAL(EB), POINTER, DIMENSION(:) :: MASSFLUX  ! Masse flux in Kg.m-2.s-2
   REAL(EB), POINTER, DIMENSION(:) :: MLR       ! Masse Loss Rate in Kg/s
END TYPE

!TYPE (GPYRO_BC_TYPE), ALLOCATABLE, TARGET, DIMENSION(:,:,:,:,:) :: GPYRO_BOUNDARY_CONDITION
TYPE (GPYRO_BC_TYPE), POINTER, DIMENSION(:,:,:,:) :: GPBCP

TYPE :: GPYRO_BOUNDARYS_INFORMATION
   INTEGER , POINTER, DIMENSION(:) :: IMESH_GPYRO
   INTEGER , POINTER, DIMENSION(:) :: IZ_GPYRO
   INTEGER , POINTER, DIMENSION(:) :: IX_GPYRO
   INTEGER , POINTER, DIMENSION(:) :: IY_GPYRO
   INTEGER , POINTER, DIMENSION(:) :: IOR_GPYRO
   INTEGER , POINTER, DIMENSION(:) :: IOR_FDS
   INTEGER , POINTER, DIMENSION(:) :: IMESH_FDS ! IF 0 no coupling ;
                                                ! Positif Coupling with FDS Mesh NB IMESH_FDS;
                                                ! Negatif Coupling with an other Gpyro mesh nb IMESH_FDS
   INTEGER , POINTER, DIMENSION(:) :: I_FDS
   INTEGER , POINTER, DIMENSION(:) :: J_FDS
   INTEGER , POINTER, DIMENSION(:) :: K_FDS
   INTEGER , POINTER, DIMENSION(:) :: IW_FDS
   INTEGER , POINTER, DIMENSION(:) :: RATIO
   REAL(EB), POINTER, DIMENSION(:) :: DIFF
   REAL(EB), POINTER, DIMENSION(:) :: XDIFF
   REAL(EB), POINTER, DIMENSION(:) :: YDIFF
   REAL(EB), POINTER, DIMENSION(:) :: ZDIFF   
   REAL(EB), POINTER, DIMENSION(:) :: Z_GPYRO
   REAL(EB), POINTER, DIMENSION(:) :: X_GPYRO
   REAL(EB), POINTER, DIMENSION(:) :: Y_GPYRO
   REAL(EB), POINTER, DIMENSION(:) :: X_FDS
   REAL(EB), POINTER, DIMENSION(:) :: Y_FDS
   REAL(EB), POINTER, DIMENSION(:) :: Z_FDS
   LOGICAL , POINTER, DIMENSION(:) :: COMPLETE_CELL_AT_BC
   TYPE(GPYRO_BC_TYPE), POINTER, DIMENSION(:) :: GPBC
END TYPE
TYPE (GPYRO_BOUNDARYS_INFORMATION), ALLOCATABLE, DIMENSION (:) :: GP_BOUDARYS

TYPE :: GPYRO_TO_GPYRO_INTERFACE
   INTEGER :: NINTERFACES
   INTEGER , POINTER, DIMENSION(:) :: IMATCH
   INTEGER , POINTER, DIMENSION(:) :: IMESH1
   INTEGER , POINTER, DIMENSION(:) :: IMESH2
   INTEGER , POINTER, DIMENSION(:) :: IZ1
   INTEGER , POINTER, DIMENSION(:) :: IZ2
   INTEGER , POINTER, DIMENSION(:) :: IX1
   INTEGER , POINTER, DIMENSION(:) :: IX2
   INTEGER , POINTER, DIMENSION(:) :: IY1
   INTEGER , POINTER, DIMENSION(:) :: IY2
   INTEGER , POINTER, DIMENSION(:) :: IOR1
   INTEGER , POINTER, DIMENSION(:) :: IOR2
   REAL(EB), POINTER, DIMENSION(:) :: DLT1
   REAL(EB), POINTER, DIMENSION(:) :: DLT2
   REAL(EB), POINTER, DIMENSION(:) :: TMP1
   REAL(EB), POINTER, DIMENSION(:) :: TMP2
   REAL(EB), POINTER, DIMENSION(:) :: K1
   REAL(EB), POINTER, DIMENSION(:) :: K2
   REAL(EB), POINTER, DIMENSION(:) :: Q
END TYPE
TYPE (GPYRO_TO_GPYRO_INTERFACE) :: G2GIFACE,IFACETEMP !GPYRO TO GPYRO INTERFACE

INTEGER, ALLOCATABLE, DIMENSION (:) :: NGPYRO_FACES_NEEDING_BCS

TYPE :: GPYRO_GENERAL_TYPE
   INTEGER :: PROCESS_GPYRO(1:100) = 0, MYID
   LOGICAL :: USE_MPI

   REAL(EB), DIMENSION(1:100,1:3) :: MESH_CENTROIDS, MESH_EXTENTS, MESH_LL 

   INTEGER :: NUM_GPYRO_MESHES = 1
   INTEGER :: NIC !Number of initial conditions
   TYPE(INITIAL_CONDITION_TYPE), POINTER, DIMENSION(:) :: INITIAL_CONDITIONS

   INTEGER :: NSURF_IDX
   TYPE(ALLBC_TYPE), POINTER, DIMENSION(:) :: ALLBC

   INTEGER :: NOBST
   REAL(EB), POINTER, DIMENSION(:) :: ZDIM, XDIM, YDIM
   INTEGER, POINTER, DIMENSION(:) :: NCELLZ, NCELLX, NCELLY
   LOGICAL, POINTER, DIMENSION(:) :: HALF_CELLS_AT_BC
   CHARACTER(100), POINTER, DIMENSION(:) :: GEOMETRY_FILE
   INTEGER, POINTER, DIMENSION(:,:) :: DEFAULT_SURF_IDX
   INTEGER, POINTER, DIMENSION(:) :: DEFAULT_IC
   !REAL(EB), POINTER, DIMENSION(:) :: OFFSETZ, OFFSETX, OFFSETY !Dead option
   TYPE(GEOM_TYPE), POINTER, DIMENSION(:) :: GEOM

   REAL(EB) :: GX   !x-component of gravity vector
   REAL(EB) :: GY   !y-component of gravity vector
   REAL(EB) :: GZ   !z-component of gravity vector
      
   REAL(EB) :: GPYRO_WALL_CLOCK_START

   ! Quantities read in from GPYRO_DUMP namelist group
   CHARACTER(60) :: CASENAME
   LOGICAL  :: DUMP_ENERGY_BALANCE
   INTEGER  :: N_POINT_QUANTITIES
   INTEGER  :: N_PROFILE_QUANTITIES
   INTEGER  :: N_SMOKEVIEW_QUANTITIES
   REAL(EB) :: DTDUMP_GA  
   REAL(EB) :: DTDUMP_POINT
   REAL(EB) :: DTDUMP_PROFILE
   REAL(EB) :: DTDUMP_SMOKEVIEW
   REAL(EB) :: REDUCED_DTDUMP
   REAL(EB) :: TMP_REDUCED_DTDUMP
   
   CHARACTER(60), POINTER, DIMENSION(:) :: POINT_QUANTITY,PROFILE_QUANTITY,SMOKEVIEW_QUANTITY
   REAL(EB)     , POINTER, DIMENSION(:) :: POINT_Z,POINT_X,POINT_Y,PROFILE_COORD1,PROFILE_COORD2,SMOKEVIEW_LOCATION
   CHARACTER(1) , POINTER, DIMENSION(:) :: PROFILE_DIRECTION
   CHARACTER(2) , POINTER, DIMENSION(:) :: SMOKEVIEW_PLANE
   INTEGER      , POINTER, DIMENSION(:) :: POINT_IZ, POINT_IX, POINT_IY, POINT_IMESH, PROFILE_IMESH, SMOKEVIEW_IMESH, SMOKEVIEW_ICELL
   INTEGER      , POINTER, DIMENSION(:) :: PROFILE_ISKIP,PROFILE_IX,PROFILE_IY,PROFILE_IZ,POINT_QUANTITY_INDEX,PROFILE_QUANTITY_INDEX, &
                                           SMOKEVIEW_QUANTITY_INDEX
   LOGICAL  :: DUMP_DETAILED_CONVERGENCE

   
   REAL(EB) :: TAMB !Ambient temperature
   
   REAL(EB) :: RHOINF
         
   INTEGER  :: NCASES  !Number of cases to run
   REAL(EB) :: DT0     !Initial timestep 
   REAL(EB) :: DTNEXT  !Next timestep
   REAL(EB) :: P0      !Background pressure
   
   REAL(EB) :: MAXTMP !Maximum temperature (solid)
   REAL(EB) :: MINTMP !Minimum temperature (solid)
   
   REAL(EB) :: DTMIN_KILL 

   LOGICAL :: NAN
   REAL(EB) :: POSINF
   REAL(EB) :: NEGINF
   
   REAL(EB) :: EPS
   REAL(EB) :: VHLC
         
   REAL(EB) :: HCV
   REAL(EB) :: NU_A
   REAL(EB) :: NU_B
   REAL(EB) :: NU_C
         
   LOGICAL  :: FDSMODE
   REAL(EB) :: FDS_HEAT_OF_COMBUSTION
   INTEGER  :: FDS_MATL_VER
   LOGICAL  :: CONSTANT_DHVOL
   LOGICAL  :: FULL_QSG
   LOGICAL  :: GASES_PRODUCED_AT_TSOLID
         
   LOGICAL  :: NOCONSUMPTION
   REAL(EB) :: TMPTOL
   REAL(EB) :: HTOL
   REAL(EB) :: YITOL
   REAL(EB) :: PTOL
   REAL(EB) :: YJTOL
   REAL(EB) :: HGTOL
   REAL(EB) :: EPS_YIS
   REAL(EB) :: EPS_YJG

   LOGICAL  :: ANISOTROPIC_SPECIES(1:50)

   INTEGER,  POINTER, DIMENSION(:)  :: ICASE
   INTEGER,  POINTER, DIMENSION(:)  :: IMESH
   REAL(EB), POINTER, DIMENSION(:)  :: TSTOP

   LOGICAL,  POINTER, DIMENSION(:) :: ZEROD !Is this a 0D simulation
   REAL(EB), POINTER, DIMENSION(:) :: BETA !Heating rate for TGA
   ! Logicals 
   LOGICAL  :: THERMAL_EQUILIBRIUM
   LOGICAL  :: SOLVE_GAS_ENERGY
   LOGICAL  :: SOLVE_PRESSURE
   LOGICAL  :: SOLVE_GAS_YJ
   LOGICAL  :: NEED_GAS_YJ ! If not SOLVE_GAS_YJ YJ=YJ_INITIAL
   LOGICAL  :: EXPLICIT_T
   LOGICAL  :: USE_TOFH_NEWTON
   LOGICAL  :: SHYI_CORRECTION
   LOGICAL  :: BLOWING
   LOGICAL  :: CONVENTIONAL_RXN_ORDER
   LOGICAL  :: KOZENY_CARMAN
   LOGICAL  :: USE_ANISOTROPIC_SOLID_ENTHALPY_SOLVER
   LOGICAL  :: GEOMETRY_IS_UPSIDE_DOWN
   LOGICAL  :: SOLVE_POROSITY
         
   ! Globals 
   REAL(EB) :: TREF
   REAL(EB) :: TDATUM ! temperature datum for enthalpy
   INTEGER  :: NTDMA_ITERATIONS !For TDMA algorithm
   INTEGER  :: NSSPECIESITERNS,NCONTINUITYITERNS,NCOEFF_UPDATE_SKIP
   REAL(EB) :: ALPHA
   REAL(EB) :: ALPHA_YIS
   REAL(EB) :: ALPHA_YJG
   REAL(EB) :: ALPHA_H
   REAL(EB) :: ALPHA_P
   REAL(EB) :: ALPHA_HG
   REAL(EB) :: DT
   INTEGER  :: CONV_DIFF_SCHEME
   REAL(EB) :: TORTUOSITY_FACTOR
   LOGICAL  :: USE_TORTUOSITY_FACTOR_FOR_FLUX
   LOGICAL  :: USE_CONSTANT_HM0
   REAL(EB) :: CONSTANT_HM0

   ! FDS related parameters
   REAL(EB) :: YY0(1:5)
   REAL(EB) :: TDUMPLAST_GA
   REAL(EB), DIMENSION(1:100) :: TDUMPLAST_SMOKEVIEW, TDUMPLAST_POINT, TDUMPLAST_PROFILE, TDUMPLAST_CSVDUMP !Supports up to 100 meshes
 
   ! Timing
   REAL(EB) :: TUSED(0:60)=0.0_EB
   ! Newton iteration stuff
   REAL(EB), DIMENSION(1:100) :: NEWTON_A,NEWTON_B,NEWTON_C,NEWTON_D,NEWTON_E,NEWTON_F,NEWTON_G

   ! Test for fds coupling
   REAL(EB), DIMENSION(1:100) :: DTNEXT_ARR = 9D9

   ! Used for multi-mesh i/o
   INTEGER, DIMENSION(1:100) :: FIRST_SF_IN_MESH, FIRST_PROF_IN_MESH

   ! Gpyro to Gpyro coupling
   REAL(EB) :: GPYRO_TO_GPYRO_TOLERANCE


   ! Solver control
   CHARACTER(4) :: SOLVER_MAIN_DIRECTIONS
   LOGICAL ::  FIX_DOMAIN_TEMPERATURE

   LOGICAL :: SHUTDOWN_ALREADY_TRIGGERED =.FALSE. !Flag that help shodwn gracefully with OMP

END TYPE

! This is a terrible name for a variable that is used so much, but
! below GPG stands for "GPyro, General". For the longest time, this variable was 
! called FG for "FIST, general" reflecting that this model was originally developed 
! as part of the FIST project...

TYPE (GPYRO_GENERAL_TYPE), TARGET, SAVE :: GPG

! Air properties table
REAL(EB), DIMENSION(1:26,0:8) :: AIR_PROPS_TABLE

! *****************************************************************************
END MODULE GPYRO_VARS
! *****************************************************************************

! *****************************************************************************
MODULE OMP_UTILS
! *****************************************************************************
USE OMP_LIB
IMPLICIT NONE

CONTAINS

FUNCTION MY_OMP_GET_NUM_THREADS() RESULT(N)
INTEGER :: N

N = 1
!$OMP PARALLEL
!$OMP MASTER
!$ IF (OMP_GET_NUM_THREADS() /= 0) THEN
!$    N = OMP_GET_NUM_THREADS()
!$ ENDIF
!$OMP END MASTER
!$OMP BARRIER
!$OMP END PARALLEL

END FUNCTION MY_OMP_GET_NUM_THREADS


FUNCTION OMP_STATUS_MESSAGE() RESULT(MESSAGE)
CHARACTER(LEN=100) :: MESSAGE
INTEGER :: N
N = MY_OMP_GET_NUM_THREADS()
IF (N <= 1) THEN
   MESSAGE = ' OpenMP Disabled'
ELSE
   WRITE(MESSAGE, '(A,I3)') ' OpenMP Enabled; Number of Threads: ', N
END IF
END FUNCTION OMP_STATUS_MESSAGE
! *****************************************************************************
END MODULE OMP_UTILS
! *****************************************************************************
