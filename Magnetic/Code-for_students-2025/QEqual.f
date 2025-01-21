      FUNCTION QEqual(diff,x,y) result(qeq)

cc Answer if x and y are different by less than diff

      implicit none
      logical*1 :: qeq
      real*8, intent(in) :: diff, x, y
CCCCCCCCCC
      if (abs(???).le.abs(???)) then
       qeq=.true.
      else
       qeq=.false.
      endif
CCCCCCCCCC
      End FUNCTION QEqual

