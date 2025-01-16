program main
implicit none

real*8,allocatable :: H(:,:), D(:,:), autov(:)
real*8, allocatable :: work(:)
real*8 :: alpha, beta_p, beta_m

integer :: i,j,k, starting,ending
integer :: ii, tmp
integer :: n
integer :: info, lwork
integer,allocatable :: indices(:) 

logical :: cyclic

character*64 :: fomt
character*4 :: string

write(*,*) 'How many carbon atoms?'
read(*,*) n
alpha = 0.d0
lwork = 3*n - 1
allocate(work(lwork))



allocate(H(n,n), D(n,n), autov(n))
allocate(indices(n))



write(*,*) 'Set value of beta(+) (remember, IT IS NEGATIVE):'
read(*,*) beta_p
write(*,*) ''

write(*,*) 'Set value of beta(-) (remember, IT IS NEGATIVE):'
read(*,*) beta_m

write(*,*) 'Is it cyclic (t/f)?'
read(*,*) cyclic


! Inizializzo H

H = 0

do i=2,n-1
    if (mod(i,2).eq.1) then
        H(i+1,i) = beta_p
        H(i-1,i) = beta_p
    else
        H(i+1,i) = beta_m
        H(i-1,i) = beta_m
    endif
    H(i,i) = alpha
enddo
H(2,1) = beta_p

if (mod(n,1).eq.1) then
    H(n-1,n) = beta_p
else
    H(n-1,n) = beta_m
endif



H(1,1) = alpha
H(n,n) = alpha



if (cyclic) then
    H(1,n) = H(n-1,n)
    H(n,1) = H(2,1)
endif





write(*,*) '========================================================='
write(*,*) 'Which eigenvectors and eigenvalues to you want to print?'
write(*,*) ''
write(*,*) 'Write one number/word at time, when you are done type "."'

indices = 0
ii = 1
do 
    read(*,*) string
    
    if (string.eq.'.'.or. ii.eq.n+1) exit
    
    if (string.eq.'all') then
        do i = 1,n
            indices(i) = i
        enddo
        exit
    endif

    if (string.eq.'from') then
        read(*,*) starting
        write(*,*) 'To?'
        read(*,*) ending
        tmp = starting
        do i = ii, (ending-starting)+1
            indices(i) = tmp
            tmp= tmp+ 1 
        enddo
        ii = ii + 1
        cycle
    endif
    read(string, '(I3)') tmp
    indices(ii) = tmp
    ii = ii + 1
enddo



! printing hamiltonian
write(fomt, '(A,I0,A)') '(', n, 'F10.3)'
! DIAGONALIZZAZIONE
D=H
call dsyev('V', 'U', n, D, n, autov, work, lwork, info)

! Controllo errori
    if (info /= 0) then
        print *, 'Errore nella diagonalizzazione: info =', info
        stop
    end if


open(unit=12,file='Results.out')
open(unit=13,file='plot.txt')
write(12,*) 'Undiagonalized Huckel Hamiltonian:'
do i=1,n
    write(12,fomt) H(i,:)
enddo

write(12,*) ''
write(12,*) '==========================================='
write(12,*) ''



write(12,*) 'Autovettori e autovalori:'
do i = 1, n
    if (indices(i).eq.0) exit
    write(12,*) '-------------------------------------------'
    write(12,'(A18,1X,I3,A1)') 'Autovettore numero ',indices(i),':'
    write(12,fomt) D(indices(i),:)
    write(12,'(A15, 1X, F10.3)') 'Con autovalore:', autov(indices(i))
enddo

write(13,*) '# there are ', n,' atoms'
write(13,*) '# index, eigenvectors, eigenvalues' 
do i=1,n
    write(13,*) i-1, D(i,:), autov(i)
enddo

close(13)
close(12)


deallocate(H)
deallocate(D)
deallocate(autov)
deallocate(work)



end program
