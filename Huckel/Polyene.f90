program main
implicit none

real*8,allocatable :: H(:,:), D(:,:), autov(:)
real*8, allocatable :: work(:)
real*8 :: alpha, beta

integer :: i,j,k
integer :: n
integer :: info, lwork

logical :: cyclic

character*64 :: fomt

alpha = 0.d0
lwork = 3*n - 1
allocate(work(lwork))


write(*,*) 'How many carbon atoms?'
read(*,*) n


write(*,*) 'Set value of beta:'
read(*,*) beta

allocate(H(n,n), D(n,n), autov(n))

write(*,*) 'Is it cyclic (t/f)?'
read(*,*) cyclic


! Inizializzo H

H = 0

do i=1,n-1
    H(i,i+1) = beta
    H(i+1,i) = beta
    H(i,i) = alpha
enddo

H(n,n) = alpha
if (cyclic) then
    H(1,n) = beta
    H(n,1) = beta
endif



! printing hamiltonian
write(fomt, '(A,I0,A)') '(', n, 'F10.3)'

open(unit=12,file='Hamilotonian.txt')
do i=1,n
    write(12,fomt) H(i,:)
enddo

write(12,*) ''
write(12,*) '==========================================='
write(12,*) ''

! diagonalizzazione
D=H
call dsyev('V', 'U', n, D, n, autov, work, lwork, info)

! Controllo errori
    if (info /= 0) then
        print *, 'Errore nella diagonalizzazione: info =', info
        stop
    end if

write(12,*) 'Autovalori:'
write(12,fomt) autov 

write(12,*) ''
write(12,*) '==========================================='
write(12,*) ''

write(12,*) 'Autovettori:'
do i = 1, n
    write(12,*) '-------------------------------------------'
    write(12,*) 'Numero ',i,':'
    write(12,fomt) D(i,:)
enddo



close(12)


deallocate(H)
deallocate(D)
deallocate(autov)
deallocate(work)



end program
