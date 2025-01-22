!-----------------!
 module parameters
!-----------------!

 real(8), parameter :: length=5.12d0   ! length of the box (in Bohr)
 real(8), parameter :: mass = 1822.88839d0 != 1 atomic mass unit (=1 g/mol)
 real(8), parameter :: pi=3.141592653589793d0
 real(8), parameter :: au2kcalmol=627.509d0
 real(8), parameter :: fs2au=41.341373336561d0
 real(8) :: angfreq, barrier
 character(10) :: potentialtype

 end module parameters
!---------------------!

!-----------------!
 program propagate
!-----------------!
 use parameters
 implicit none

 integer :: npoints,ntime,snapshot,i
 real(8) :: alpha,dt,t,dx,x0
 real(8), allocatable :: pot(:),kin(:),psisquare(:)
 complex(8), allocatable :: psi(:),psi0(:),exppot(:),expkin(:)

 open(unit=10,file='wavepacket')
   read(10,*) npoints               !Number of lattice points
   read(10,*) x0                    !Initial position
   read(10,*) alpha                 !Governs the initial width of the wave packet
   read(10,*) dt                    !Propagation time step
   read(10,*) ntime                 !Number of propagation steps
   read(10,*) snapshot              !snapshot frequency 
 close(10)

 open(unit=11,file='potential')
   read(11,*) potentialtype         !harmonic or double well potential
   read(11,*) angfreq               !Angular frequency for harmonic potential
   read(11,*) barrier               !Height of barrier in double well potential (in kcal/mol)
 close(11)

 dt=dt*fs2au                        !convert femtoseconds to atomic units
 angfreq=angfreq/fs2au              !convert femtoseconds to atomic units

 allocate(psi(npoints),psi0(npoints))
 allocate(pot(npoints),exppot(npoints))
 allocate(kin(npoints),expkin(npoints))
 allocate(psisquare(npoints))

 dx=length/dble(npoints)

 call initpsi(npoints,dx,alpha,x0,psi0)              !Obtain initial wavepacket psi0
 call fourier(0,npoints,psi0)                        !Initialize the FFT
 call operators(npoints,dx,dt,pot,kin,exppot,expkin) !Calculate the kinetic and potential operators

 psi=psi0                                            !Set the wavepacket psi at t=0 equal psi0
! open(13,file='Psi_sol.txt')
! do i = 1,npoints
!     write(13,*) (i-1)*dx-length/2., real(psi0(i)), aimag(psi0(i))
! enddo
 do i=0,ntime                                        !Start propagation
    t=i*dt
    if (i>0) then
       psi=psi*exppot                                !Multiply psi with exp(-i*dt*potential)
       call fourier(1,npoints,psi)                   !Forward FFT to momentum space
       psi=psi*expkin                                !Multiply psi with the exp(-i*dt*kinetic operator)
       call fourier(-1,npoints,psi)                  !Backward FFT to position space
    endif
    if (mod(i,snapshot)==0) then                     !Take a snapshot if the remainder of i/snapshot equals 0
       call initgraph(i/snapshot,t)                  !Initialize graph
       psisquare=(abs(psi))**2
       call graphpot(dx,npoints)                     !Plot the potential
       call graphpsi(dx,npoints,psisquare)           !Plot |psi|^2
    endif
 end do                                              !End propagation

 deallocate(psi,psi0)
 deallocate(pot,exppot)
 deallocate(kin,expkin)
 deallocate(psisquare)

 end program propagate
!---------------------!

!------------------------------------------------!
 subroutine initpsi(npoints,dx,alpha,x0,psi0)
!------------------------------------------------!
 implicit none

 integer :: i,j,npoints
 real(8) :: alpha,x,x0,dx
 complex(8) :: psi0(npoints)
 open(13,file='Psi_sol.txt')
! do i = 1,npoints
!     write(13,*) (i-1)*dx-length/2., real(psi0(i)), aimag(psi0(i))
! enddo

 do i=-npoints/2+1,npoints/2
    x=dble(i)*dx
    if (i>0) then
       j=i
    else     
       j=i+npoints
    endif
    psi0(j)=exp(-alpha*(x-x0)**2)

    write(13,*)  x, real(psi0(j)), aimag(psi0(j))
 end do

 
 
 end subroutine initpsi
!----------------------!

!------------------------------------------------------!
 subroutine operators(npoints,dx,dt,pot,kin,exppot,expkin)
!------------------------------------------------------!
 use parameters
 implicit none

 integer :: i,j,npoints
 real(8) :: x,p,b,dt,dx,dp,pot(npoints),kin(npoints)
 complex(8) :: exppot(npoints),expkin(npoints)

 dp=2.d0*pi/length
open(13,file='Energy_his.txt')
do i=-npoints/2+1,npoints/2
    x=dble(i)*dx
    p=dble(i-1)*dp
    if (i>0) then
       j=i
    else
       j=i+npoints
    endif
    if (potentialtype=='harmonic') then
       pot(j)=0.5d0*mass*angfreq**2*x**2
    elseif (potentialtype=='doublewell') then
       pot(j)=barrier*(16.d0*x**4 - 8.d0*x**2 + 1.d0)/au2kcalmol
    endif
    kin(j)=0.5d0*p**2/mass
    exppot(j)=exp(-dt*(0,1)*pot(j))
    expkin(j)=exp(-dt*(0,1)*kin(j))
    write(13,*) x, pot(j), p, kin(j)
 end do
close(13)
 end subroutine operators
!------------------------!

!-----------------------------------!
 subroutine fourier(dir,npoints,psi)
!-----------------------------------!
 implicit none

 integer :: i,npoints,dir
 real(8) :: nr
 complex(8) :: psi(npoints)
 real(8), allocatable, save :: wsave(:)       

 if (dir==1) then
    call dcfftf(npoints,psi,wsave)
    nr=1.d0/dble(npoints)
    do i=1,npoints
       psi(i)=psi(i)*nr
    end do
 elseif (dir==-1) then
    call dcfftb(npoints,psi,wsave)
 elseif (dir==0) then
    if (allocated(wsave)) deallocate(wsave)
    allocate(wsave(4*npoints+20))
    call dcffti(npoints,wsave)
 endif

 end subroutine fourier
!----------------------!
