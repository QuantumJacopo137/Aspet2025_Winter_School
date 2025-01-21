      PROGRAM Heisenberg
      implicit none
      logical*1 :: debug
      logical*1, allocatable :: tmpdeter(:)
      logical*1, allocatable :: det(:,:), det1(:), det2(:)
      integer :: i, j, k, l, dim, nsites, nal, ndeter, factorial
      integer :: nms_spaces, numspace
      integer, allocatable :: tab(:), size_ms_subspace(:)
      integer :: accu, compt, na_min, sizemax
      real*8 ValCoupl, ms
      real*8, allocatable :: coupling(:,:), mat(:,:), eigval0(:)
      real*8, allocatable :: energies(:,:)
 

Ccccccccccccccccccccccccccccccc
C    Reading
      open(11, file="input", form="formatted")
C    Debug logical, true for verbose output
      read (11,*) debug
C    Number of sites
      read (11,*) nsites
C    Allocating coupling matrix
      allocate (coupling(nsites,nsites))
      if(debug) then
       write(6,*) 'Calling the reading of the coupling matrix'
      endif
      call readinput(debug,nsites,coupling)
      close(11)
Ccccccccccccccccccccccccccccccc
c write the informations read in the external file. Check that it is OK
      write(6,*) 'Coupling matrix between sites i and j'
      do i=1,nsites ! questo
       write(6,101,advance='no') i
      enddo
      write(6,*)
      do i=1,nsites
       write(6,102,advance='no') i
       do j=1,nsites
        write (6,103,advance='no') coupling(i,j) ! anche lui
       enddo
       write(6,*)
      enddo
Ccccccccccccccccccccccccccccccc
c Size of the largest Ms subspace
!      na_min= int(mod(nsites,2)) + nsites/2     !!! minimum number of alpha spins, is a function of nsites
      na_min = (nsites+1)/2
      sizemax= factorial(nsites)
      sizemax= sizemax / (factorial(sizemax-na_min))
      sizemax = sizemax / (factorial(na_min))
      !!! size of the largest matrix, is a function of nsites and na_min
      !nms_spaces= (nsites+1)/2 + mod(nsites,2) !!! number of different ms spaces, is a function of nsites
      nms_spaces = nsites /2  +1 
      if (debug) then
          write(*,*) 'na_min =', na_min
          write(*,*) 'size max =', sizemax
          write(*,*) 'nms spaces =', nms_spaces
          stop
      endif
ccccccccccccccccccccccccccccccc
      write(6,*) '*********************'
      write(6,84) 'Number of sites        :', nsites
      write(6,84) 'Number of ms spaces    :', nms_spaces
      write(6,85) '#alpha in Ms min, size :', na_min, sizemax
      write(6,*) '*********************'
ccccccccccccccccccccccccccccccc
c  !!! Understand what contains tmpdeter, energies and size_ms_subspace arrays
      allocate (tmpdeter(nsites),det1(nsites),det2(nsites))
      allocate (energies(nms_spaces+1,0:sizemax))
      allocate (size_ms_subspace(nms_spaces))
c    !!! You have to initialize the arrays  size_ms_subspace and energies
Ccccccccccccccccccccccccccccccc
c     Loop over the possible alpha/beta=nsites-alpha repartition of the spins
      write(6,*) '******************************************'
      write(6,*) 'Loop over the spaces of determinants of decreasing
     &Ms value'
      numspace=0
      do nal=nsites,nsites-nms_spaces+1,-1
       numspace=numspace+1
       ms=(2*nal-nsites)/2.0
       energies(numspace,0)=ms
       ndeter=factorial(nsites)/(factorial(nsites-nal)*factorial(nal))
       allocate (tab(ndeter),eigval0(ndeter))
       allocate (det(ndeter,nsites))
       write(6,*) '******************************************'
       write(6,*) 
       write(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       write(6,83) '!!! Space number:',numspace  , ' !!! ' !!! Print the correct variable
       write(6,80)'Nb alph/bet spin, Ms Value:', nal, nsites-nal, ms !!!  Print the correct variables
       write(6,81)'Number of determinants          :', ndeter !!! Print the correct variable
       do i=1,ndeter
        tab=0
       enddo
       accu=0
       compt=0
ccc  Generation of all the determinants of a given Ms value
       call deter(nsites,nal,tab,ndeter,accu,compt)  !!! This subroutine is given
ccc   Convert the above integer in a suite of nsites boolean
       call Integ_to_Bit(nsites,ndeter,tab,det) !!! This subroutine is given
       write(6,*) '*******'
       write(6,*) 'Logical value of the determinants'
       do i=1,ndeter !!! Complete the end of the loop
        write(6,*) det(i,:)  !!! Print the correct arrray
       enddo
       write(6,*) '*******'
ccc  Fill the Hamiltonian matrix
       if(debug) then
        write(6,*) 'Generation of the Hamiltonian matrices'
       endif
       allocate (mat(ndeter,ndeter)) !!! What is the dimension of this array?
       do i=1,ndeter !!! Complete the end of the loop
       do j=i,ndeter !!! Complete the end of the loop
    	 do k=1,nsites !!! Complete the end of the loop
    	  det1(k)=det(i,k)
    	  det2(k)=det(j,k)
    	 enddo
C Call the subroutine that fills the matrix
         mat(i,j)=ValCoupl(debug,nsites,i,j,det1,det2,coupling) !!! You have to write this subroutine
	     mat(j,i)=mat(i,j) !!! The matrix is hermitic
	   enddo
       enddo   
C A nice printing of the Heisenberg Matrix
       write(6,*) 'The Hamiltonian matrix in the det basis'
       write(6,104,advance='no') ' Deter   '
       write(6,105,advance='no') ' '
       do i=1,ndeter
C        write(6,111,advance='no') i
        do l=1,nsites
         write(6,107,advance='no') det(i,l)
        enddo
        write(6,105,advance='no') '  '
       enddo
       write(6,*)
       do i=1,ndeter
C        write(6,112,advance='no') i
        do l=1,nsites
         write(6,107,advance='no') det(i,l)
        enddo
        write(6,105,advance='no') '  '
        do j=1,ndeter
         write (6,113,advance='no') mat(j,:) ! NOT SURE
        enddo
        write(6,*)
       enddo
CCCCCCCCCC
ccc  Diagonalize the Hamiltonian matrix
       if(debug) then
        write(6,*) 'Call of the diagonalization'
       endif
       call diasym(mat,eigval0,ndeter)   !!! It is given
cc energies contains the energies of each subspace
cc size_ms_subspace contains the size of each subspace
       do i=1,ndeter
       energies(numspace,i)=eigval0(i)
      enddo
      size_ms_subspace(numspace)=ndeter
ccc
C     Writing of the results
       write(6,*) 'Energies and decomposition of the wave functions in
     &the determinant basis set (in line)'
       write(6,104,advance='no') ' Energy   '
       write(6,105,advance='no') ' '
       do k=1,ndeter
        do l=1,nsites
         write(6,107,advance='no') det(k,l)
        enddo
        write(6,105,advance='no') '  '
       enddo
       write(6,*)
       do k=1,ndeter
        write(6,99,advance='no') eigval0(k), '    '
        write(6,100) (mat(l,k), l=1,ndeter) ! EVEN THIS
       enddo
ccc
       deallocate(tab,mat,det,eigval0)
      enddo
      deallocate(tmpdeter,det1,det2)
ccc  Summary
      if (debug) then
       do i=1,nms_spaces
        write(6,*) 'Space num', i, 'size',size_ms_subspace(i)
        write(6,*) ' Ms and Energies'
        do j=0, size_ms_subspace(i)
         write(6,99) energies(i,j)
        enddo
       enddo
      endif    
ccc  Associate the S value to the energy
      if (debug) then
       write(6,*) 'Entering the subroutine that identify the Ms value'
      endif
      call identify_Ms(debug,nms_spaces,size_ms_subspace,energies)  !!!  You have to write this subroutine. You can follow the way proposed in identify_MS.f
c The procedure will give the S value in the last column of array "energies"
cccccccccccccccccccccccccc
      write(6,*) '******************************************'
      write(6,*) 'Final results'
      Write(6,*) 'Summary: Ms value and energies in this subspace'
      do i=1,nms_spaces
       write(6,61,advance='no') energies(i,0)
       do j=1,size_ms_subspace(i)
        write(6,62,advance='no') energies(i,j)
       enddo
       write(6,*)
      enddo
      write(6,64,advance='no') ' '

      Write(6,*) 'S values and below energy'
      do j=1,size_ms_subspace(nms_spaces)
       write(6,63,advance='no') energies(nms_spaces+1,j)
      enddo
      write(6,*)
      do j=1,size_ms_subspace(nms_spaces)
       write(6,62,advance='no') energies(nms_spaces,j)
      enddo
      write(6,*)

      write(6,*)
      open(13, file="E_and_S_output", form="formatted")
      write(13,*) 'Energie, S'
      do j=1,size_ms_subspace(nms_spaces)
       write(13,65) energies(nms_spaces,j),energies(nms_spaces+1,j)
      enddo
      close(13)
Ccccccccccccccccccccccccccccccc

61    format(F5.1)
62    format(F12.3)
63    format(F12.1)
64    format(a5)
65    format(F12.6,F6.1)
71    format(a15)
72    format(F4.2)
73    format(i4)
74    format(F10.6)
80    format(a34,i3,i3,F5.1)
81    format(a34,i3)
83    format(a17,i3,a5)
84    format(a23,i5,a1)
85    format(a23,i5,i5)
99    format(F8.4,a2)
100   format(50F10.4)
101   format(i10)
102   format(i2)
103   format(50F10.4)
104   format(a10)
105   format(a4,2x)
107   format(L)
111   format(i10)
112   format(i2)
113   format(50F10.4)
10000 end PROGRAM Heisenberg
