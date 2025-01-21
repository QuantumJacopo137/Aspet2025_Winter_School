program main
implicit none

! INPUT VARIABLES
logical :: debug
integer :: n_sites, n_couplings
integer, allocatable :: sites(:,:)
real*8, allocatable :: coupling(:)
real*8, allocatable :: HRM(:,:)


integer :: i,j,k

open(unit=13, file='input.mag')
read(13,*) debug
read(13,*) n_sites
read(13,*) n_couplings

allocate(sites(n_couplings,2), coupling(n_couplings))

do i = 1, n_couplings
    read(13,*) sites(i,1), sites(i,2), coupling(i)
enddo


allocate(HRM(n_sites,n_sites))
HRM = 0.d0

do i = 1,n_sites
    HRM(sites(i,1),sites(i,2)) = coupling(i)
enddo


















end program
