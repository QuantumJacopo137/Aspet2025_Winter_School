subroutine RHF(nBas,nO,S,T,V,Hc,ERI,X,ENuc,EHF,e,c)

! Perform a restricted Hartree-Fock calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer :: i

  integer,intent(in)            :: nO ! number of occupied 
  double precision,intent(in)   :: S(nBas,nBas) ! overlap
  double precision,intent(in)   :: T(nBas,nBas) ! Kinetic
  double precision,intent(in)   :: V(nBas,nBas) ! Potential
  double precision,intent(in)   :: Hc(nBas,nBas)! Nuclear Hamiltonian
  double precision,intent(in)   :: X(nBas,nBas) ! X matrix
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas) 
  double precision,intent(in)   :: ENuc

! Read number of atoms, number of electrons of the system
! nO   = number of occupied orbitals
! nV   = number of virtual orbitals (see below)
! nBas = number of basis functions (see below)
!      = nO + nV


  double precision :: tmp_matrix(nBas,nBas)

! Local variables

  integer,parameter             :: maxSCF = 64
  double precision,parameter    :: thresh = 1d-5
  integer                       :: nSCF
  double precision              :: Conv
  double precision              :: Gap
  double precision              :: ET,EV,EJ
  double precision              :: EK
  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: P(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: F(:,:),Fp(:,:)
  double precision,allocatable  :: error(:,:)
  
  double precision,external     :: trace_matrix
  double precision     :: trace_JK
  double precision     :: trace_H

! Output variables

  double precision,intent(out)  :: EHF
  double precision,intent(out)  :: e(nBas)
  double precision,intent(out)  :: c(nBas,nBas)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|      Restricted Hartree-Fock calculation     |'
  write(*,*)'************************************************'
  write(*,*)

! Memory allocation

  allocate(cp(nBas,nBas),P(nBas,nBas),      &
           J(nBas,nBas),K(nBas,nBas),F(nBas,nBas),Fp(nBas,nBas), &
           error(nBas,nBas))

! MATRIX DENSITY

  P = 0.d0 !  



! Guess coefficients and eigenvalues

  F(:,:) = Hc(:,:)

! Initialization

  nSCF = 0
  Conv = 1d0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------

  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'| RHF calculation                                  |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','HF energy','|','Conv','|','HL Gap','|'
  write(*,*)'----------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

!   Increment 

    nSCF = nSCF + 1                                             ! ITERATION NUMBER
    Fp = matmul(transpose(X),matmul(F,X))                       ! FIND THE F PRIME
    
    call diagonalize_matrix(nBas, Fp, e)                        ! DIAGONALIZE F PRIME AND FIND C PRIME
    
    cp = Fp                                                     ! THE ARRAY WITH THE EIGENVECTORS IS SET TO CP
    C = matmul(X,cp)                                            ! I CALCULATE C
    
    P = 2.d0 * matmul(C(:,1:nO),transpose(C(:,1:nO)))                              ! UPDATE P
    
    call compute_F(P,Hc,ERI,nBas,F,J,K)                         ! CALCULATE F IN THE MEAN TIME 

    
    !conv = maxval(abs(matmul(F,matmul(P,S))- matmul(S,matmul(P,F)))) ! USE THIS NEW F TO CHECK FOR CONVERGENCE
    conv = maxval(abs(matmul(F,matmul(P,S))-transpose(matmul(F,matmul(P,S))) )) ! USE THIS NEW F TO CHECK FOR CONVERGENCE
    
    tmp_matrix = matmul(P,Hc)
    trace_H = trace_matrix(nBas,tmp_matrix)                     ! COMPUTE THE P*H TRACE

    tmp_matrix = matmul(P,F)
    trace_JK = trace_matrix(nBas, tmp_matrix)                   ! COMPUTE 0.5 * P*(J+K) TRACE

    EHF =0.5d0 * (trace_H  + trace_JK)                          ! COMPUTE THE EHF ENERGY WITH THE NEW F
    Gap = e(nO+1)- e(nO)

!   Dump results
!    write(*,*) 'F matrix at iteration ', nSCF
!    do i = 1,nBas
!        write(*,*) F(i,:)
!    enddo

!    write(*,*) 'P at iteration ', nSCF
!    do i = 1,nBas
!        write(*,*) P(i,:)
!    enddo

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
      '|',nSCF,'|',EHF+ENuc,'|',Conv,'|',Gap,'|'
 
  enddo
  write(*,*)'----------------------------------------------------'
!------------------------------------------------------------------------
! End of SCF loop
!------------------------------------------------------------------------

! Did it actually converge?

  if(nSCF == maxSCF) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    stop

  endif

! Compute final HF energy

  call print_RHF(nBas,nO,e,C,ENuc,ET,EV,EJ,EK,EHF)

end subroutine RHF





subroutine compute_F(P,Hc,ERI,nBas,F,J,K)
implicit none
  integer, intent(in) :: nBas
  double precision,intent(in)   :: Hc(nBas,nBas)! Nuclear Hamiltonian
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: P(nBas,nBas) ! Density matrix
  double precision   :: F(nBas,nBas) ! Density matrix



  double precision   :: J(nBas,nBas) ! Coulomb matrix
  double precision   :: K(nBas,nBas) ! Exchange matrix

  integer :: mu, nu, lambda,sigma

  do mu = 1, nBas
      do nu = 1, nBas
          J(mu,nu) = 0.d0
          K(mu,nu) = 0.d0
          do lambda = 1,nBas
              do sigma = 1, nBas
                  J(mu,nu) = J(mu,nu) + P(lambda,sigma) * ERI(mu,lambda,nu,sigma)
                  K(mu,nu) = K(mu,nu) + 0.5d0 * P(lambda,sigma) * ERI(mu,lambda,sigma,nu)
              enddo
          enddo
          F(mu,nu) = Hc(mu,nu) + J(mu,nu) - K(mu,nu)

      enddo
  enddo
                



end subroutine

