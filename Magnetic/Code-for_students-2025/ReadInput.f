      SUBROUTINE ReadInput(debug,nsites,coupling)

cc This subroutine reads the input coupling matrix
cc if debug is true the program is verbose
cc nsites is the number of sites
cc coupling is the output coupling matrix

      implicit none
      logical*1 debug
      integer :: i, j, c1, c2, ncoupl, nsites
      real*8 :: coupl
      real*8 coupling(nsites,nsites)

CCCCCCCCCC
      if(debug) then
       write(6,*) 'Entering in ReadInput'
      endif
C      open(11, file="input", form="formatted")
C    Initialization of the coupling matrix
      do i=1,nsites
       do j=1,nsites
        coupling(i,j)=0.0d0
       enddo
      enddo  
CCCCCCCCCC
C    Number of couplings to be read
      read (11,*) ncoupl
C    Reading of the coupling matrix
c    The couplings are given as
c    Site1 Site2 coupling_intensity ... N times  
      do i=1,ncoupl
       read (11,*) c1, c2, coupl
       coupling(c1,c2)=coupl
       coupling(c2,c1)=coupl
      enddo
      if(debug) then
       write(6,*) 'End of ReadInput'
      endif
CCCCCCCCCC
      End SUBROUTINE ReadInput
