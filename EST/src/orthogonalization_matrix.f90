subroutine orthogonalization_matrix(nBas,S,X)
implicit none

integer :: nBas
integer :: i,j 
double precision :: S(nBas,nBas), X(nBas,nBas), eigen(nBas), S_D(nBas,nBas), U(nBas,nBas)
double precision :: XT(nBas,nBas)
logical :: debug

debug = .false.


U = S


call diagonalize_matrix(nBas,U,eigen)


eigen(:) = 1.d0/dsqrt(eigen(:))
S_D=0.d0

call ADAt(nBas,U,eigen,X)




!XT = matmul(X,matmul(S,transpose(X)))
!do i = 1,nBas
!    write(*,*) XT(i,:)
!enddo



if (debug) then
write(*,*) 'S'
do i = 1,nBas
    write(*,*) S(i,:)
enddo

write(*,*) 'S_D'
do i = 1,nBas
    write(*,*) S_D(i,:)
enddo



write(*,*) 'X'
do i = 1,nBas
    write(*,*) X(i,:)
enddo

XT= transpose(X)

endif


end subroutine


!diagonalize_matrix(N,A,e)

