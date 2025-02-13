program calculate_pi

implicit none

real*8 :: x,y,somme
integer :: nbr_pts,n_inside,n_total
integer :: clock, n
integer, allocatable :: seed(:)
integer :: i, pid
integer :: getpid

write(*,*), "combien de test pour converger"
read(*,*), nbr_pts

!!!! Generation Seed
! We do a lot of things here instead of simply call random_seed()
! If you qant to know why, ask J. Cuny

call random_seed(size=n)

allocate(seed(n))

call SYSTEM_CLOCK(COUNT=clock)
pid = getpid()
seed = clock + pid*57 + 37 * (/ (i - 1, i = 1, n) /)
call random_seed(put=seed)  

!!!! End of Generation Seed


n_inside=0
n_total=0

do i=1,nbr_pts

  call random_number(x)
  call random_number(y)

  somme=sqrt(x**2+y**2)
  
  if (somme.le.1) n_inside=n_inside+1

  n_total=n_total+1

  if (modulo(i,10)==0) then
    write(*,*) 4.0d0* dble(n_inside)/dble(n_total)
  endif
enddo

write(*,*) "Final pi"
write(*,*) 4.0d0* dble(n_inside)/dble(n_total)

end program calculate_pi
