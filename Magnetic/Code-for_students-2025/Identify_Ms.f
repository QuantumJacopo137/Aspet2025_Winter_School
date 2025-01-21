      SUBROUTINE identify_Ms(debug,nms_spaces,size_ms_subspace,energies)

cc Find the Ms associated to the different energies
cc debug: verbose mode
cc nms_spaces
cc size_ms_subspace
cc energies

      implicit none
      logical*1 :: debug,qequal,found
      integer, intent(in) :: nms_spaces , size_ms_subspace(nms_spaces)
      real*8, intent(out) ::
      c energies(nms_spaces+1,0:size_ms_subspace(nms_spaces))
      real*8, allocatable :: mstab1(:), mstab2(:), tmpmstab(:)
      real*8 diff
      integer :: i,j,k,l,compt
CCCCCCCCCC
      allocate (tmpmstab(size_ms_subspace(nms_spaces)))
c
c The vector with the S value is initialize to the max Ms value
c If two reals have a difference smaller than diff we consider that they are equal
      diff=1E-8
      do i=1,size_ms_subspace(nms_spaces)
       tmpmstab(i)=energies(1,0)
      enddo
      if(debug) then
       write(6,*) '###############################'
       write(6,*) 'In subroutine Identify_Ms'
       write(6,*) 'Number of space', nms_spaces
       write(6,*) 'Sizes',(size_ms_subspace(j), j=1,nms_spaces)
      endif
ccccc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



C Loop over the Ms subspaces to compare the energies of an Ms space (variable i)
C to the energies of the space Ms-1 (i+1)
C mstab1 contains the S values of the energies of i subspace
C mstab2 contains the S values of the energies of i+1 subspace
      do i=1,nms_spaces-1
       allocate (mstab1(size_ms_subspace(i))
       allocate (mstab2(size_ms_subspace(i+1))
C For Ms max, only one determinant whose S value is Ms max
       if (i.eq.1) then
        mstab1(1)=energies(1,0) ! MAYBE I DON'T KNOW PLEASE E-SCEMAMA HELP ME
        
        else
C mstab1 gets the value obtained in the previous loop over Ms
        do j=1,size_ms_subspace(i)
         mstab1(j)=tmpmstab(j)
        enddo
       endif
       if (debug) then
        write(6,*) '*****************'
        write(6,98) 'comparing space=',i,' ms=', energies(i,0), 
     c             ' size=', size_ms_subspace(i)
        write(6,98) '      and space=',i+1,' ms=', energies(i+1,0),
     c             ' size=', size_ms_subspace(i+1)
        do k=1,size_ms_subspace(i)
         write(6,99,advance='no') mstab1(k)
        enddo
        write(6,*)
       endif
ccccc







C k is the index of the state in the i+1 subspace
       k=1

C j is the index of the state in the i subspace
       do j=1,size_ms_subspace(i)
        if (debug) then
         write(6,*) '*********'
	 write(6,*) 'i,j', i,j
	 write(6,97) 'Energy of the first subspace', energies(i,j)
	 write(6,97) 'S value of the first subspace',mstab1(j)
        endif 
        found=.false.

C looks inside the second space an energy identical to the jth ones in
C the first space, till it finds it



    do while (.not.found)
	 if (debug) then
	  write(6,*) 'i+1,k', i+1,k
	  write(6,97) 'Energy of the second subspace', energies(i+1,k)
	 endif 
	 if (k==size_ms_subspace(nms_spaces)+2) then !E-SCEMAMA PLEASE 
	  write(6,*) 'Strange, the value was not found. I stop'
	  stop
	 endif
C if the kth value of the second space is not equal to the jth value of
C the first one, then
C          the S value of the kth energy  [mstab2(k)]
C          is equal [=] 
C          to the Ms of the
C          second space

C          k is incremented by one
	if (.not.abs(energies(i+1,k)-energies(i,j)).le.diff) then ! E-SCEMAMA PLEASE
	  mstab2(k)=energies(i+1,0) ! E-SCEMAMA PLEASE
	
      if(debug) write(6,*) 'Different energies, k', k  ! USELESS


C if the kth value of the second space is = to the jth value of the
C first one then 
C    the S value of kth value of space 2 is = to the jth value of space 1 
C    k is incremented by one and we go out the do while loop (j->j+1)
     else       
	    mstab2(k) = mstab1(j)
	    found=.true. ! HERE YOU CHECK IF ITS THERE (mstab2(k)==mstab1(j))

        if(debug) then
          write(6,*) 'Same energies, k', k
          write(6,*) 'S value ',mstab1(j)
        endif
	endif


C Incrementation of k
	    k=k+1

        enddo
       enddo




C if the last energy of the first space has been associated to an energy
C of the second state but there are still energies in the 2nd space then
C       all the remaining energies of space 2 have an S value equal to
C       the Ms value of the second space
       
       if (k.le.size_ms_subspace(i+1)+1) then ! SCEMAMAMAMAMAMAMAMAMAMMAAAM
        do l=k,size_ms_subspace(i+1)
         mstab2(l)=energies(i+1,0)
        enddo
        if (debug) then
         write(6,*)'End of the energies of the first set'
         write(6,*)'The S value of the remaining energies of the second'
         write(6,*)'set are then ', energies(i+1,0)
        endif
       endif
       if (debug) then
        write(6,*) 'S values of the first subspace'
        do k=1,size_ms_subspace(i)
         write(6,99,advance='no') mstab1(k)
        enddo
        write(6,*)
        write(6,*) 'S values of the second subspace'
        do k=1,size_ms_subspace(i+1)
         write(6,99,advance='no') mstab2(k)
        enddo
        write(6,*)
       endif

C the table that contains, at the end, the S value of all the states,
C gets the value of the second space and will be used, in the next step
C of the loop as the values of the first space
       do k=1,size_ms_subspace(i+1)
        tmpmstab(k)=mstab2(k)
       enddo
       deallocate (mstab1,mstab2)
      enddo
C the last column of the energies table contains the S value
      do i=1,size_ms_subspace(nms_spaces)
       energies(nms_spaces+1,i)=tmpmstab(i)
      enddo
c
97    format(a30,F8.3)
98    format(a16,i3,a5,F8.1,a6,i3)
99    format(F5.1,a2)
      End SUBROUTINE Identify_Ms
