      subroutine Integ_to_Bit(nsites,ndeter,tab,det)
       implicit none
       integer  , intent(in)  :: nsites, ndeter
       integer  , intent(in)  :: tab(ndeter)
       logical*1, intent(out) :: det(ndeter,nsites)
       integer :: i, k, l
        
        do i=1,ndeter
         do k=1,nsites
          det(i,nsites-k+1)=btest(tab(i),k-1)
         enddo
C         write(6,*) tab(i), (det(i,l), l=1,nsites)
        end do
       end

