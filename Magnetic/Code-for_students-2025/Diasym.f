      SUBROUTINE diasym(a,eig,n)

cc This subroutine call the diagonalization subroutine
cc a is the input matrix and also the output diagonalized matrix
cc eig contains the eigenvalures (energies)
cc n is the dimension of the matrix

      implicit none
      integer n,l,inf
      real*8  a(n,n),eig(n),work(n*(3+n/2))
      l=n*(3+n/2)
      call dsyev('V','U',n,a,n,eig,work,l,inf)
      end SUBROUTINE diasym
