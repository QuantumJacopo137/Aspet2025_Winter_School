      FUNCTION ValCoupl(debug,nsites,i,j,det1,det2,coupling) result(val)

cc Gives the coupling between two determinants, either identical either different
cc debut: verbose mode
cc nsites is the number of sites
cc i is the number of the first determinant
cc j is the number of the second determinant
cc det1 is the first determinant
cc det2 is the second determinant
cc coupling is the coupling matrix between sites

      implicit none
      logical*1 :: debug, det3(nsites)
      logical*1, intent(in) :: det1(nsites), det2(nsites)
      integer, intent(in) :: i, j
      integer :: nsites, k, l, diffspin, position_diff(2)
      real*8, intent(in) :: coupling(nsites,nsites)
      real*8 :: val
CCCCCCCCCC
      val=0.0
      diffspin=0
      if (i==j) then
       if(debug) then
        write(6,*) 'Diagonal term'
       endif
       ! devi sommare tutti i termini di coupling per le coppie uguali
       do k = 1,nsites-1
           do l = k,nsites
            if (det1(k).neqv.det1(l)) then
                val = val + coupling(k,l)
            endif
           enddo
       enddo
    

               
       !?????????????????
       !?????????????????
       !?????????????????
      else
       diffspin = 0
       do k = 1,nsites   
        det3(k)=xor(det1(k),det2(k))
           if (det3(k)) then
               diffspin = diffspin + 1
               position_diff(diffspin) = k
           endif 
           if (diffspin.ge.3) then
               val = 0
               return
           endif
       enddo
       
       val = coupling(position_diff(1),position_diff(2))

       if(debug) then
        write(6,*) 'Extradiagonal term'
        write(6,*) 'det1', det1, 'det2', det2, 'xor', det3
       endif

!       ?????????????????
!       ?????????????????
!       ?????????????????
!       You may use the xor function (look what it gives) in the following way
!       ?????????????????
!       ?????????????????
!       ?????????????????
      endif

CCCCCCCCCC
      End FUNCTION ValCoupl

