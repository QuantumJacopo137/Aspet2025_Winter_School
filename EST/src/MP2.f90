subroutine MP2(nBas,nO,nV,e,ERI,ENuc,EHF)

! Compute the HF energy in the MO basis

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EHF
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: i,j
  integer                       :: a,b
  double precision              :: Dijab
  double precision              :: E2d
  double precision              :: E2x
  double precision              :: EcMP2

! Output variables

! Hello world

  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'|          MP2 calculation           |'
  write(*,*)'**************************************'
  write(*,*)

!--------------------------------!
! Compute MP2 correlation energy !
!--------------------------------!

!------------!
! MP2 energy !
!------------!

  write(*,*)
  write(*,'(A32)')           '--------------------------'
  write(*,'(A32)')           ' MP2 calculation          '
  write(*,'(A32)')           '--------------------------'
  write(*,'(A32,1X,F16.10)') ' MP2 correlation energy = ',EcMP2
  write(*,'(A32,1X,F16.10)') ' Direct part            = ',2d0*E2d
  write(*,'(A32,1X,F16.10)') ' Exchange part          = ',-E2x
  write(*,'(A32)')           '--------------------------'
  write(*,'(A32,1X,F16.10)') ' MP2 electronic  energy = ',EHF + EcMP2
  write(*,'(A32,1X,F16.10)') ' MP2 total       energy = ',ENuc + EHF + EcMP2
  write(*,'(A32)')           '--------------------------'
  write(*,*)

end subroutine MP2 
