subroutine AO_to_MO(nBas,c,ERI_AO,ERI_MO)

! Expression of bi-electronic integrals in the MO basis set

  implicit none
!  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: c(nBas,nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: mu,nu,la,si
  integer                       :: p,q,r,s
  double precision,allocatable  :: scr(:,:,:,:)

! Output variables

  double precision,intent(out)  :: ERI_MO(nBas,nBas,nBas,nBas)

! Memory allocation
  
!--------------------------------------
! AO to MO transformation starts here !
!--------------------------------------

end subroutine AO_to_MO
