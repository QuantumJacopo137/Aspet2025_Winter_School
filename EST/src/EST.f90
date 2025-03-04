program EST

! Electronic Structure Theory program

  include 'parameters.h'

  integer                       :: nAt,nBas,nEl,nO,nV
  double precision              :: ENuc,EHF

  double precision,allocatable  :: ZNuc(:)
  double precision,allocatable  :: rAt(:,:)

  integer                       :: nShell
  integer,allocatable           :: TotAngMomShell(:)
  integer,allocatable           :: KShell(:)
  double precision,allocatable  :: CenterShell(:,:)
  double precision,allocatable  :: DShell(:,:)
  double precision,allocatable  :: ExpShell(:,:)

  double precision,allocatable  :: S(:,:)
  double precision,allocatable  :: T(:,:)
  double precision,allocatable  :: V(:,:)
  double precision,allocatable  :: Hc(:,:)
  double precision,allocatable  :: X(:,:)
  double precision,allocatable  :: ERI(:,:,:,:)
  double precision,allocatable  :: ERI_MO(:,:,:,:)
  double precision,allocatable  :: e(:)
  double precision,allocatable  :: c(:,:)

  double precision              :: start_HF,end_HF,t_HF
  double precision              :: start_AOtoMO,end_AOtoMO,t_AOtoMO
  double precision              :: start_MP2,end_MP2,t_MP2
  double precision              :: start_CCD,end_CCD,t_CCD

! Hello World

  write(*,*)
  write(*,*) '***************************'
  write(*,*) '* TCCM winter school 2023 *'
  write(*,*) '***************************'
  write(*,*)

!------------------------------------------------------------------------
! Read input information
!------------------------------------------------------------------------

! Read number of atoms, number of electrons of the system
! nO   = number of occupied orbitals
! nV   = number of virtual orbitals (see below)
! nBas = number of basis functions (see below)
!      = nO + nV

  call read_molecule(nAt,nEl,nO)
  allocate(ZNuc(nAt),rAt(nAt,3))

! Read geometry

  call read_geometry(nAt,ZNuc,rAt,ENuc)

  allocate(CenterShell(maxShell,3),TotAngMomShell(maxShell),KShell(maxShell), &
           DShell(maxShell,maxK),ExpShell(maxShell,maxK))

!------------------------------------------------------------------------
! Read basis set information
!------------------------------------------------------------------------

  call read_basis(nAt,rAt,nBas,nO,nV,nShell,TotAngMomShell,CenterShell,KShell,DShell,ExpShell)

!------------------------------------------------------------------------
! Read one- and two-electron integrals
!------------------------------------------------------------------------

! Memory allocation for one- and two-electron integrals

  allocate(S(nBas,nBas),T(nBas,nBas),V(nBas,nBas),Hc(nBas,nBas),X(nBas,nBas), &
           ERI(nBas,nBas,nBas,nBas),e(nBas),c(nBas,nBas))

! Read integrals

  call read_integrals(nBas,S,T,V,Hc,ERI)  

!------------------------------------------------------------------------
! Orthogonalization X = S^(-1/2)
!------------------------------------------------------------------------

 call orthogonalization_matrix(nBas,S,X) ! FATTO ! **********************************

!------------------------------------------------------------------------
! Compute restricted HF energy
!------------------------------------------------------------------------

    call cpu_time(start_HF)
    call RHF(nBas,nO,S,T,V,Hc,ERI,X,ENuc,EHF,e,c) ! QUESTO
    call cpu_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for HF = ',t_HF,' seconds'
    write(*,*)

!------------------------------------------------------------------------
! AO to MO transformation
!------------------------------------------------------------------------

    allocate(ERI_MO(nBas,nBas,nBas,nBas))

    call cpu_time(start_AOtoMO)
!   call AO_to_MO(nBas,c,ERI,ERI_MO)
    call cpu_time(end_AOtoMO)

    t_AOtoMO = end_AOtoMO - start_AOtoMO
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for AO to MO transformation = ',t_AOtoMO,' seconds'
    write(*,*)
!------------------------------------------------------------------------
! Compute MP2 energy
!------------------------------------------------------------------------

    call cpu_time(start_MP2)
!   call MP2(nBas,nO,nV,e,ERI_MO,ENuc,EHF)
    call cpu_time(end_MP2)

    t_MP2 = end_MP2 - start_MP2
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for MP2 = ',t_MP2,' seconds'
    write(*,*)

!------------------------------------------------------------------------
! Compute CCD energy
!------------------------------------------------------------------------

    call cpu_time(start_CCD)
!   call CCD(nBas,nO,nV,ENuc,EHF,e,ERI_MO)
    call cpu_time(end_CCD)

    t_CCD = end_CCD - start_CCD
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for CCD = ',t_CCD,' seconds'
    write(*,*)

!------------------------------------------------------------------------
! End of EST
!------------------------------------------------------------------------
end program EST
