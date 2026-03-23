! *****************************************************************************
MODULE PREC 
! *****************************************************************************
IMPLICIT NONE
! Just like in FDS...
INTEGER, PARAMETER :: FB = SELECTED_REAL_KIND(6)
INTEGER, PARAMETER :: EB = SELECTED_REAL_KIND(12)

REAL(EB), PARAMETER :: EPSILON_EB=EPSILON(1._EB)                
REAL(FB), PARAMETER :: EPSILON_FB=EPSILON(1._FB)
END MODULE PREC