program main
implicit none

real*8 :: x,y,z
real*8 :: r(3)
real*8 :: a
real*8,parameter :: dx = 0.1

real*8, external :: e_loc
real*8 :: psi

real*8 :: E_m = 0.d0
real*8,allocatable :: E_L(:,:)
real*8 :: W=0.d0
real*8 :: sigma, sigma_t

integer :: i,j,k,tmp, ii
integer :: boundery




boundery = int(5/dx)


allocate(E_L((1+2*boundery)**3,4))




open(unit=13,file='sigma_a.dat')
write(13,*) '# a;   E_m;    sigma^2'
do ii = 1,10
a = real(ii)*0.2
E_L = 0.d0
tmp = 1
W = 0.d0
E_m = 0.d0
do i = -boundery , boundery
    do j =-boundery , boundery
        do k = -boundery , boundery
            if (i.eq.0.and.j.eq.0.and.k.eq.0) cycle
            x = real(i)*dx
            y = real(j)*dx
            z = real(k)*dx
            r = [x,y,z]
            E_L(tmp,1) = x
            E_L(tmp,2) = y
            E_L(tmp,3) = z
            E_L(tmp,4) = e_loc(a,r)
            E_m = E_m + E_L(tmp,4)  * psi(a,r) * psi(a,r) * dx*dx*dx 
            W = W + psi(a,r) * psi(a,r) * dx *dx *dx
            tmp = tmp + 1
        enddo
    enddo
enddo


E_m = E_m / W

            
!write(*,*) '<E> = ',  E_m 

sigma = 0.d0
open(unit=12,file='E_Ls.debug')
do i = 1,tmp-1
!    write(12,*) E_L(i,:)
    r = [E_L(i,1),E_L(i,2),E_L(i,3)]
    sigma = sigma + abs(psi(a,r)*psi(a,r)) * (E_L(i,4)-E_m)*(E_L(i,4)-E_m) * dx * dx * dx

enddo
close(12)

sigma = sigma / W

!write(*,*) 'sigma^2 =', sigma 


write(13,*) a, E_m, sigma

enddo
close(13)

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






