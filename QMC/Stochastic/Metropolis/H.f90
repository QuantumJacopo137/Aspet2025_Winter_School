program main
implicit none

real*8 :: x,y,z
real*8 :: ux, uy, uz ! random direction
real*8 :: r(3), u(3), r_tmp(3)
real*8 :: a
real*8,parameter :: dx = 0.1

real*8, external :: e_loc
real*8 :: psi

real*8 :: E_m = 0.d0
real*8,allocatable :: E_L(:,:)
real*8 :: W=0.d0
real*8 :: sigma, sigma_t
real*8,allocatable :: E_statisitcal(:) ! each vector element has the average energy of the n_max points of one MC run 

real*8 :: amplitude
real*8 :: real_random


integer :: i,j,k,tmp, ii
integer :: boundery

integer, parameter :: n_max=100000 ! config per step
integer, parameter :: m = 30
real*8, parameter  :: dL = 0.1 ! walk step 



allocate(E_L(n_max,4))
allocate(E_statisitcal(m))


a = 1.2

do j = 1,m
! initial conditions
E_L = 0.d0
tmp = 1
W = 0.d0
E_m = 0.d0



call random_seed()


call random_number(x)
call random_number(y)
call random_number(z)


r = [x-0.5d0,y-0.5d0,z-0.5d0] ! first step
r = r * 10.d0
!write(*,'(A20,I3)') 'Started step number ', j
!write(*,*) r




do i = 1,n_max
!    do 
    call random_number(ux)
    call random_number(uy)
    call random_number(uz)
    
    ux = (ux-0.5d0) *2.d0
    uy = (uy-0.5d0) *2.d0
    uz = (uz-0.5d0) *2.d0
    u = [ux,uy,uz] 

    r_tmp= r + u * dL
    
    amplitude = min(1., (psi(a,r_tmp)*psi(a,r_tmp))/(psi(a,r)*psi(a,r)))
    call random_number(real_random)

    if (real_random.le.amplitude) then
        r = r_tmp
 !       exit
    endif
    

  !  enddo

    E_L(tmp,1) = r(1)
    E_L(tmp,2) = r(2)
    E_L(tmp,3) = r(3)
    E_L(tmp,4) = e_loc(a,r)
    
    E_m = E_m + E_L(tmp,4)!  * psi(a,r) * psi(a,r) 
    W = W + psi(a,r) * psi(a,r) 
    
    tmp = tmp + 1
enddo


!E_m = E_m / W
E_m = E_m / n_max
E_statisitcal(j) = E_m

enddo





E_m = sum(E_statisitcal) / m

sigma = 0.d0
do i = 1,m
    sigma = sigma + (E_m-E_statisitcal(i))*(E_m-E_statisitcal(i))
enddo


sigma = sigma / (m-1)
    



write(*,*)
write(*,*) '***************************************************************************'
write(*,'(A24, F6.3)') 'Alpha value used: ', a
write(*,'(A24, I8)') 'Points per simulation: ', n_max
write(*,'(A24, I8)') 'Simulations: ', m

write(*,'(A24, F12.6, 2X, A3, 1X, F8.6)') 'True total <E> = ', E_m, '+/-', sqrt(sigma/M)
write(*,'(A24, F6.4)') '\sigma^2 = ', sigma
write(*,*) '***************************************************************************'



open(unit=14, file='Segment_walk')
do i = 1, n_max-1
    write(14,*) E_L(i,1), E_L(i,2), E_L(i,3)
    write(14,*) E_L(i+1,1), E_L(i+1,2), E_L(i+1,3)
    write(14,*) ''
enddo
close(14)










end program




real*8 function psi(a,r)
implicit none
real*8 :: a, dist
real*8 :: r(3)


dist = r(1)*r(1)+r(2)*r(2)+r(3)*r(3)
if (dist.lt.0) then
    write(*,*) 'Error: less than zero argument'
    stop
endif
psi = exp(-a* sqrt(dist))
endfunction


real*8 function kinetic(a,r)
implicit none
real*8 :: a, dist
real*8 :: r(3)



dist = r(1)*r(1)+r(2)*r(2)+r(3)*r(3)
dist = sqrt(dist) 
if (dist.lt.0) then
    write(*,*) 'Error: less than zero argument'
    stop
endif




kinetic= -0.5d0*(a*a-((2*a)/dist))

endfunction


real*8 function potential(r)
implicit none
real*8 :: r(3)
real*8 :: dist


integer :: i,j,k


dist = r(1)*r(1)+r(2)*r(2)+r(3)*r(3)

if (dist.eq.0) then
    write(*,*)  'Error: 0 distance in the denominator'
    stop
endif

if (dist.lt.0) then
    write(*,*) 'Error: less than zero argument'
    stop
endif


potential =- 1/sqrt(dist)

endfunction

real*8 function e_loc(a,r)
  implicit none
  double precision, intent(in) :: a, r(3)

  double precision, external :: kinetic
  double precision, external :: potential

  e_loc = kinetic(a,r) + potential(r)

end function e_loc

subroutine test_potential
    implicit none
    double precision :: r(3)
    double precision :: expected_output
    double precision, external :: potential

    expected_output = -1.d0/15.d0

    r(:) = (/ 2.d0, 5.d0, 14.d0 /)
    if (potential(r) /= expected_output) stop 'Failed'

    r(:) = (/ 5.d0, 14.d0, 2.d0 /)
    if (potential(r) /= expected_output) stop 'Failed'

    r(:) = (/ -2.d0, 5.d0, -14.d0 /)
    if (potential(r) /= expected_output) stop 'Failed'

    r(:) = (/ 5.d0, -14.d0, -2.d0 /)
    if (potential(r) /= expected_output) stop 'Failed'

    r(:) = (/ 0.d0, 9.d0, 12.d0 /)
    if (potential(r) /= expected_output) stop 'Failed'

    r(:) = (/ 9.d0, -12.d0, 0.d0 /)
    if (potential(r) /= expected_output) stop 'Failed'

    r(:) = 0.d0
    expected_output = -huge(1.d0)
    if (potential(r) /= expected_output) stop 'Failed r=0'
    print *, 'potential ok'

end subroutine test_potential






