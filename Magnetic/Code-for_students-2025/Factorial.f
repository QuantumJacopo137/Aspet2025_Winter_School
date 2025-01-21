      recursive function factorial(n)  result(f)

cc Give, recursively, the value of n!=n*(n-1)*(n-2)*...*1

        integer, intent(in) :: n
        integer :: f
        
        if(n == 0) then
          f=1
          return
        else
          f=n*factorial(n-1)
        end if
      end FUNCTION factorial
      
