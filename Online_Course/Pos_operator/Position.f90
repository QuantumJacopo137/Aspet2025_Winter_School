program main
    implicit none

    integer :: dimention, uni, int_rndm
    integer :: i, j

    real(8) :: starting_point
    real(8) :: step
    real(8), allocatable :: X(:,:),P(:,:), comm(:,:), Psi(:), Psi2(:), Pos(:)
    real(8) :: rndm, pi =3.14159265358979

    character(1) :: ssave
    character(50) :: fomt

    ! Input dall'utente
    write(*,*) 'How many points do you want?'
    read(*,*) dimention

    write(*,*) 'Type the step: '
    read(*,*) step

    write(*,*) 'Type the starting point: '
    read(*,*) starting_point

    do
        write(*,*) ''
        write(*,*) 'You want the output on the terminal or in a file? Type ''T'' for terminal or ''F'' for file'
        read(*,*) ssave

        select case (ssave)
            case ('T', 't')
                uni = 6 ! Terminale
                exit
            case ('F', 'f')
                uni = 11 ! File
                exit
            case default
                write(*,*) 'Not valid entry, try again!'
        end select
    end do

    ! Alloco memoria
    allocate(Psi(dimention))
    allocate(X(dimention, dimention))
    ALLOCATE(P(dimention, dimention))
    ALLOCATE(Comm(dimention, dimention))
    ALLOCATE(Pos(dimention))
    allocate(Psi2(dimention))
    

! ==================================
    open(unit=12, file='To_plot.txt')
    ! Inizializza la matrice
    X = 0.d0
    P = 0.d0
    do i = 1, dimention
        X(i, i) = starting_point + (i - 1) * step
    end do

    do i = 2,dimention-1
        P(i-1,i) = -1
        P(i+1,i) = +1
    enddo
    P(2,1) = 1
    P(dimention-1,dimention) = -1
    P = P/(2*step)

    ! Costruisco il formato dinamico
    write(fomt, '(A,I0,A)') '(', dimention, 'F10.4)'

    comm = matmul(X,P)-matmul(P,X)
    
!    call random_seed()
!    call random_number(rndm)
    Psi = 0.d0
    int_rndm = int((dimention-1)*rndm)
!    Psi(int_rndm) = 1


    ! Apro file se necessario
    if (uni == 11) open(unit=11, file='Output.txt')


    ! Scrivo la matrice
    write(uni,*) '=====$ POSITION OPERATOR $======'
    do i = 1, dimention
        write(uni, fomt) X(i, :)
    end do
    write(uni,*) '=====$ MOMENTUM OPERATOR $======'
    do i = 1, dimention
        write(uni, fomt) P(i, :)
    end do
    write(uni,*) '=====$     COMMUTATOR    $======'
     do i = 1, dimention
        write(uni, fomt) comm(i, :)
    end do
    ! Chiudo il file se aperto

    Pos=0.d0
    do i = 1,dimention
        Pos(i) = pi*X(i,i) / X(dimention,dimention)
        Psi(i) = sin(Pos(i))
    enddo
    
    Psi2 = matmul(comm, Psi)
    
    do i = 1, dimention
        write(12,*) Pos(i), Psi(i), Psi2(i) 
    enddo

write(uni,*) '=====$        PSI        $======'
    write(uni,fomt) Psi
write(uni,*) '=====$        PSI2        $======'
    write(uni,fomt) Psi2

    if (uni == 11) close(11)

    ! Dealloco la memoria
    deallocate(Psi)
    deallocate(X)
    


 close(12)
end program

