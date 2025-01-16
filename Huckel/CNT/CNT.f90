program main
implicit none

real*8,allocatable :: H(:,:), D(:,:), autov(:)
real*8, allocatable :: work(:)
real*8 :: alpha, beta

integer :: i,j,k, starting,ending
integer :: ii, tmp
integer :: n, len
integer :: info, lwork
integer,allocatable :: indices(:) 

logical :: cyclic

character*64 :: fomt
character*4 :: string

write(*,*) 'How many C atoms in the circumference?'
read(*,*) n
write(*,*)
write(*,*) 'How many units is the CNT long?'
read(*,*) len



alpha = 0.d0
lwork = 3*len*n - 1
allocate(work(lwork))



allocate(H(len*n,len*n), D(len*n,len*n), autov(len*n))
allocate(indices(len*n))



write(*,*) 'Set value of beta (remember, IT IS NEGATIVE):'
read(*,*) beta


cyclic = .true.

! Inizializzo H

H = 0

! easy, just a bunch of closed polyenes
do j=0,len-1
do k=1,n-1
    i = k+j*n
    H(i,i+1) = beta
    H(i+1,i) = beta
    H(i,i) = alpha
enddo

H(n,n) = alpha

H(j*n+1,j*n+n) = beta
H(j*n+n,j*n+1) = beta

enddo


! now the fun part, the interaction between two adjacent belts
if (len.ge.2) then ! at least 2 layers to be a nantube, duh!
    do j = 0,len-2
        do i = 1, n, 2 
            H(j*n+i, (j+1)*n+1+mod(i,n)) = +beta !CHANGE ME BACK
            H((j+1)*n+1+mod(i,n),j*n+i) = +beta
        enddo
    enddo


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
write(fomt, '(A,I0,A)') '(', len*n, 'F10.3)'
! DIAGONALIZZAZIONE
D=H
call dsyev('V', 'U', len*n, D, len*n, autov, work, lwork, info)

! Controllo errori
    if (info /= 0) then
        print *, 'Errore nella diagonalizzazione: info =', info
        stop
    end if


open(unit=12,file='Results.out')
open(unit=13,file='plot.txt')
write(12,*) 'Undiagonalized Huckel Hamiltonian:'
do i=1,n*len
    write(12,fomt) H(i,:)
enddo


open(unit= 15, file='Hamilotonian.txt')
do i=1,n*len
    write(15,fomt) H(i,:)
enddo
close(15)

write(12,*) ''
write(12,*) '==========================================='
write(12,*) ''



write(12,*) 'Autovettori e autovalori:'
do i = 1, n*len
    if (indices(i).eq.0) exit
    write(12,*) '-------------------------------------------'
    write(12,'(A18,1X,I3,A1)') 'Autovettore numero ',indices(i),':'
    write(12,fomt) D(indices(i),:)
    write(12,'(A15, 1X, F10.3)') 'Con autovalore:', autov(indices(i))
enddo


open(unit=14, file = 'E_n.txt')
write(13,*) '# there are ', n,' atoms'
write(13,*) '# index, eigenvectors, eigenvalues' 
do i=1,n*len
    write(13,*) i-1, D(i,:)
    write(14,*) i, autov(i)
enddo

close(14)
close(13)
close(12)


deallocate(H)
deallocate(D)
deallocate(autov)
deallocate(work)



end program
