program main
implicit none

real*8 :: x,y,z
real*8 :: r(3)
real*8 :: a
real*8 :: dx = 10.d0/(50.-1.) 
real*8 :: E_a(20)
real*8, external :: e_loc

integer :: i,j,k


a = 0.d0
r = 0.d0

open(unit=11,file='E_x_a.dat')
do i = 1, 50
    r(1) = -5.d0 + (i-1)*dx
do j = 1,20
    a = j* 0.1d0
    E_a(j) = e_loc(a,r) 
enddo
    write(11,*) r(1), E_a(:)
enddo
close(11)


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






