module utilities

public::lecture_xyz             ! Read and save positions in an .xyz file
public::get_distances           ! Calculate square of distances between atoms
public::get_forces              ! Calculate forces from the positions
public::get_pot_ener            ! Calculate energy from the positions
public::output_structure        ! Write output structures 
public::output_energy           ! Write output energies
public::get_kin_ener            ! Calculate kinetic energy from velocities
public::get_T                   ! Calculate T from kinetic energy

contains

subroutine lecture_xyz (filename,tableau,nbr_atomes)
implicit none


integer :: nbr_atomes, i
character*32 :: filename
real*8, allocatable :: tableau(:,:)
character*2 :: tmp
logical :: debug=.false.


open(13, file=filename)
read(13,*) nbr_atomes
read(13,*)

allocate(tableau(nbr_atomes,3))

do i = 1,nbr_atomes
    read(13,*) tmp, tableau(i,:)
enddo
close(13)

if (debug) then
    open(111, file='Debug_input.txt')
    do i = 1,nbr_atomes
        write(111,*) tableau(i,:)
    enddo
endif




endsubroutine lecture_xyz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate square of distances between atoms !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_distances (positions,distances,nbr_atomes)
implicit none

integer :: nbr_atomes, i ,j, k
real*8 :: positions(nbr_atomes,3)
real*8 :: distances(nbr_atomes,nbr_atomes)
real*8 :: diff
logical :: debug = .false.


if (debug) then
open(111, file='Debug_input.txt')
    do i = 1,nbr_atomes
        write(111,*) positions(i,:)
    enddo
close(111)
endif

distances(:,:) = 0.d0
do i = 1,nbr_atomes
    do j = 1,nbr_atomes
        do k = 1,3
            diff = positions(i,k) - positions(j,k)
            diff = diff*diff
            distances(i,j) = distances(i,j) + diff
        enddo
!        distances(j,i) = distances(i,j)
    enddo
enddo

if (debug) then
open(111,file='Distance_debug.txt')
do i= 1,nbr_atomes
!    do j = 1,nbr_atomes
!        write(111,'(F8.5,2X)', advance = 'no') distances(i,j)
!    enddo
    write(111,*) distances(i,:)
    write(111,*)
enddo
close(111)
endif


return


endsubroutine get_distances

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate LJ forces               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_forces (forces,distances,positions,nbr_atomes,sigma,eps)
implicit none
integer :: nbr_atomes, i ,j, k
real*8 :: sigma, eps
real*8 :: positions(nbr_atomes,3)
real*8 :: distances(nbr_atomes,nbr_atomes)
real*8 :: r_ij(nbr_atomes,nbr_atomes)
real*8 :: U(nbr_atomes,nbr_atomes)
real*8 :: forces(nbr_atomes, 3)
real*8 :: frac
real*8 :: tmp
logical :: debug = .false.


r_ij(:,:) = sqrt(distances(:,:))
U = 0.d0
do i = 1,nbr_atomes-1
    do j = i,nbr_atomes
        if (i.eq.j) cycle
        frac = sigma / r_ij(i,j)
        U(i,j) = 24 * (eps/r_ij(i,j)) * (frac**6-2*frac**12)
    enddo
enddo

forces = 0.d0
do i = 1,nbr_atomes
    do j = 1,nbr_atomes
        if (i.eq.j) cycle
        do k = 1,3
            tmp = positions(i,k)-positions(j,k)
            forces(i,k) =forces(i,k) - U(i,j) * tmp / r_ij(i,j)
        enddo
    enddo
enddo



if (debug) then
    open(112,file='Forces_debug.txt', position='append')
    do i =1,nbr_atomes
        write(112,*) forces(i,:), '   ',r_ij(i,:)
    enddo
    write(112,*)
    close(112)
endif
return




endsubroutine get_forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate LJ Energy               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_pot_ener (distances,nbr_atomes,sigma,eps,energy )
implicit none
integer :: nbr_atomes, i,j
real*8 :: sigma,eps
real*8 :: energy, kost
real*8 :: distances(nbr_atomes,nbr_atomes)


energy = 0.d0
kost =4 * eps * sigma**6

do i = 1,nbr_atomes-1
    do j = i+1, nbr_atomes
        energy = energy + kost * ((sigma**6/ distances(i,j)**6)- (1/distances(i,j)**3))
    enddo
enddo


return

endsubroutine get_pot_ener

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Output gemetry                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine output_structure (index,positions,nbr_atomes,energy_pot)
implicit none
integer :: index, nbr_atomes, i
real*8 :: positions(nbr_atomes,3)
real*8 :: energy_pot



if (index.eq.1) then
open(index, file = 'md.xyz', status='replace')
else 
    open(index, file = 'md.xyz', status='old', position='append')
endif

write(index, *) nbr_atomes
write(index,*) '# Potential energy =', energy_pot
do i = 1,nbr_atomes
    write(index,*) 'Ar', positions(i,:)
enddo
write(index, *)
close(index)

endsubroutine output_structure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Ouput Energies                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine output_energy (index,index_loop,timestep,energy_pot,energy_kin,temp)
implicit none
integer :: index, index_loop,i
real*8 :: timestep
real*8 :: energy_kin
real*8 :: energy_pot
real*8 :: temp


open(index, file = 'Output_energies.xyz')
!write(index,'(A8,2X,A14,2X,A14,2X,A14,2X,A14)') '# step', 'time', 'energy pot', 'energy kin', 'temperature'
    write(index,'(I8,2X,F14.6,2X,F14.6,2X,F14.6,2X,F14.6)') index_loop, real(index_loop)*timestep, energy_pot, energy_kin, temp
close(index)
endsubroutine output_energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate Kinetic Energy          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_kin_ener (velo,nbr_atomes,kin_ener,mass)
implicit none
real*8 :: velo(nbr_atomes, 3)
real*8 :: kin_ener, v_mod
real*8 :: mass
integer :: i, k, nbr_atomes




kin_ener=0.d0

do i = 1,nbr_atomes
    v_mod = 0.d0
    do k = 1,3
        v_mod = velo(i,k)*velo(i,k) + v_mod
    enddo
    kin_ener = kin_ener + 0.5d0 * mass * v_mod
enddo



endsubroutine get_kin_ener

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate Instantaneous T         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_T (kin_ener,temp,nbr_atomes)
implicit none
integer :: nbr_atomes
real*8 :: kin_ener,temp
real*8, parameter :: k_b = 3.1668114e-6
temp = 2*kin_ener

temp = temp /(k_b*3*nbr_atomes)

endsubroutine get_T


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Apply Thermostat                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine thermostat(velo,thermo,temp,temp_0,val,nbr_atomes,timestep)
implicit none
integer :: i,j, nbr_atomes, thermo ! se thermo = 1 allora è quello di riscala altrimenti con 2 barren
real*8 :: velo(nbr_atomes,3)
real*8 :: temp,temp_0, val, timestep  ! val è il tau di accoppiamento in caso di barren
! temp0 = temp simulazione, l'altra è quella istantanea
real*8 :: lambda

if (thermo.eq.1) then
    lambda = temp_0 / temp
    lambda = sqrt(lambda)

    velo = lambda * velo
!    write(*,*) lambda
else if (thermo.eq.2) then
    lambda = temp_0 / temp
    lambda = lambda - 1
    lambda = lambda * timestep / val
    lambda = lambda +1
    lambda = sqrt(lambda)
!    write(*,*) lambda
    
    velo = lambda * velo
else
    write(*,*) 'Error: Invalid thermostat option!'
    stop
endif


return
endsubroutine thermostat

end module utilities
