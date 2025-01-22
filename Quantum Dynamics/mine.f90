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

 integer :: npoints,ntime,snapshot,i,j,tmp, frames
 real(8) :: alpha,dt,t,dx,x0,x
 real(8), allocatable :: pot(:),kin(:),psisquare(:), evolution(:,:), tenergy(:),venergy(:)
 complex(8),allocatable :: evolution_tot(:,:)
 complex(8), allocatable :: psi(:),psi0(:),exppot(:),expkin(:), tmp_psi(:)
 complex(8), allocatable :: coefficients_t(:), coefficients_e(:)
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

 dt = dt*fs2au                        !convert femtoseconds to atomic units
 dx = length/dble(npoints)
 angfreq=angfreq/fs2au              !convert femtoseconds to atomic units

 allocate(psi(npoints),psi0(npoints))
 allocate(pot(npoints),exppot(npoints))
 allocate(kin(npoints),expkin(npoints))
 allocate(psisquare(npoints))

call eval_psi(psi,npoints, length, x0, alpha)
psi0 = psi  ! salvo una copia di psi prima della propagazione per dopo


call operators(npoints,dx,dt,pot,kin,exppot,expkin)


open(13, file= 'Psi_debug.txt')
do i = 1,npoints
    write(13,*) (i-1)*dx-length, real(psi(i)), aimag(psi(i))
enddo
close(13)


! EVOLUZIONE


! evolution è un array che ha come colonne la psi^2 di ogni i-esimo step
! la colonna -1 è il valore di x del punto in question



frames = int(ntime / snapshot)

allocate(evolution(npoints, -1:frames))
allocate(evolution_tot(npoints, -1:ntime))

write(*,*) 'Dimensioni di evolution:', shape(evolution)
do i = -npoints/2+1,npoints/2
    x=dble(i)*dx
    if (i>0) then
       j=i
    else
       j=i+npoints
    endif
    evolution(j,-1) = x
    evolution(j,0) = psi0(j)*conjg(psi0(j))
enddo


evolution_tot(:,-1) = evolution(:,-1)
evolution_tot(:,0) = psi0(:)


call fourier(0,npoints,psi) ! inizializzo per allocare la corretta memoria
tmp = 1
do i = 1,ntime
    psi = psi*exppot
    call fourier(1,npoints,psi)
    psi = psi*expkin
    call fourier(-1,npoints,psi) 
    do j = 1,npoints
        psisquare(j) = dble(psi(j)*conjg(psi(j)))
    enddo
    evolution_tot(:,i) = psi(:)
    if (mod(i,snapshot).eq.0) then
        evolution(:,tmp) = psisquare(:)
        tmp = tmp + 1
    endif
enddo


open(13, file = 'Animation_points.txt')
do i = 1,npoints
    do j = -1,frames
        write(13,'(F10.5)',advance='no') evolution(i,j)
    enddo
    write(13,*)
enddo
close(13)



! energy shit

allocate(tenergy(0:ntime))
allocate(venergy(0:ntime))
allocate(tmp_psi(npoints))


venergy = 0.d0
tenergy = 0.d0

open(14,file='Mechanic_energy.txt')
do i= 0, frames
    tmp_psi(:)=evolution(:,i) 
    call fourier(1,npoints,tmp_psi)
    do j = 1,npoints
        tmp_psi(j) = tmp_psi(j) * kin(j)
    enddo

    call fourier(-1,npoints,tmp_psi)
    tenergy(i) = dot_product(evolution(:,i), tmp_psi(:))

    tmp_psi(:)=evolution(:,i) 
    do j = 1,npoints
        tmp_psi(j) = tmp_psi(j) * pot(j)
    enddo
    venergy(i) =  dot_product(evolution(:,i), tmp_psi(:))
    write(14,*) i, tenergy(i)+venergy(i)
enddo
close(14)


allocate(coefficients_t(0:ntime))

coefficients_t = 0.d0
open(17,file='C_t.txt')
do i = 0, ntime
    do j = 1, npoints
        coefficients_t(i) = coefficients_t(i) + dx * conjg(evolution_tot(j,0))*evolution_tot(j,i)
    enddo
    write(17,*) i*dt, real(coefficients_t(i)), aimag(coefficients_t(i))
enddo
close(17)


allocate(coefficients_e(0:ntime))

coefficients_e = coefficients_t

call fourier(1,ntime+1, coefficients_e)
open(18, file = 'C_e.txt')
do i = 0,ntime
    write(18,*) i, real(coefficients_e(i)), aimag(coefficients_e(i))
enddo
close(18)




end program



subroutine eval_psi(psi,npoints,leng, centre,width)
implicit none
complex(8) :: psi(npoints)
integer :: i, npoints, j
real(8) :: leng,centre, width,x, step
real(8) :: norm


step = leng/npoints

!do i = 1, npoints
!    x = -(leng/2.)+(i-1) * step
!    psi(i) = exp(-width* (x-centre)*(x-centre))
!enddo

do i=-npoints/2+1,npoints/2
    x=dble(i)*step
    if (i>0) then
       j=i
    else
       j=i+npoints
    endif
    psi(j)=exp(-width*(x-centre)**2)

    write(13,*)  x, real(psi(j)), aimag(psi(j))
 end do




norm = 0.d0


do i = 1,npoints
    norm = norm + step * psi(i)*conjg(psi(i))
enddo


psi = psi/sqrt(norm)


end subroutine




!------------------------------------------------------!
 subroutine operators(npoints,dx,dt,pot,kin,exppot,expkin)
!------------------------------------------------------!
 use parameters
 implicit none

 integer :: i,j,npoints
 real(8) :: x,p,b,dt,dx,dp
 real*8 :: pot(npoints),kin(npoints)
 complex(8) :: exppot(npoints),expkin(npoints)

 dp=2.d0*pi/length
! dp = dp / dble(npoints)


 !do i=-npoints/2+1,npoints/2
 !   x=dble(i)*dx
 !   p=dble(i-1)*dp

open(13,file='Energy_debug.txt')
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

! do i=1,npoints
!    x=dble(i-1)*dx-(length/2.)
!    p=-(dble(pi/length))+dp*(i-1)/npoints
 !   p=(i-(npoints/2)-1) * dp
   ! write(*,*) x,p
  !  if (potentialtype=='harmonic') then
  !     pot(i)=0.5d0*mass*angfreq**2*x**2
  !  elseif (potentialtype=='doublewell') then
  !     pot(i)=barrier*(16.d0*x**4 - 8.d0*x**2 + 1.d0)/au2kcalmol
  !  endif
  !  kin(i)=(0.5d0*p*p) /dble(mass)
  !  exppot(i)=exp(-dt*(0,1)*pot(i))
  !  expkin(i)=exp(-dt*(0,1)*kin(i))
  !  write(13,*) x, pot(i), p, kin(i)
 end do
close(13)
 end subroutine operators
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

