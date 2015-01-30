c-------------------------------------------------------------------------------------------------c
      subroutine romberg(func, a, b, integral, eps, nprn, nmin)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Here is my own version of romberg where the algorithm is taken from the wikipedia article on   c
c  romberg integration.  I'm not exactly sure how this compares to qromb since that one seems to  c
c  use a different algorithm but whatevs.                                                         c
c                                                                                                 c
c  func - external user supplied function                                                         c
c  a - lower bound on the integral                                                                c
c  b - upper bound on the integral                                                                c
c  integral - the result of the integral within the given error                                   c
c  eps - the user given error of the integral                                                     c
c  nprn - flag for printing out results as the integral is being done: 0 - nothing printed,       c
c    1 - print results.                                                                           c
c  nmin - the minimum number of iterations before the subroutine can finish.                      c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

c input parameters
      double precision a, b, integral, eps
      integer nprn, nmin
      
c parameters used in this subroutine only
      integer i, j, m, n, nmax
      parameter (nmax=30)
      double precision R(0:nmax,0:nmax), h, summ, func

c user supplied external function
      external func

c resetting the matrix that holds all of the results
      do i=0, nmax
        do j=0, i
          R(i,j) = 0.d0
        enddo
      enddo

c starting the construction of the R matrix      
      R(0,0) = (b-a)*(func(a)+func(b))/2.d0

c output for testing purposes
      if(nprn.gt.0) then
        write(*,fmt='(E16.8,1X)') R(0,0)
      endif
      
      n = 1
      do while(n.le.nmax)
        h = (b-a)/2.d0**n

c generate a sum of intermediate points        
        summ = 0.d0
        do i=1,2**(n-1)
          summ = summ + func(a + h*(2.d0*i-1.d0))
        enddo

c start the next line of the R matrix        
        R(n,0) = R(n-1,0)/2.d0 + h*summ

c fill out the rest of the n-th row
        do m=1, n
          R(n,m) = (4.d0**m*R(n,m-1)-R(n-1,m-1))/(4.d0**m - 1.d0)
        enddo

c output that shows the results of the nth row
        if (nprn.gt.0) then
          write(*,1010) (R(n,i),i=0,n)
 1010     format(10(E16.8,1X))
        endif

c check to see if our integral is good enough
        if ((n.ge.nmin).and.(dabs(R(n,n-1)-R(n,n)).le.eps*R(n,n))) then
          integral = R(n,n)
          return
        endif
        
        n = n+1
      enddo

c if we go through an arbitrary number of iterations we stop      
      write(*,*) 'Too many steps in Romberg'
      n = n-1
      write(*,*) 'Value of integral: ',R(n,n)
      write(*,*) 'Relative error in integral: ',dabs(R(n,n)-R(n,n-1))/R(n,n)
      pause

      integral = R(n,n)
      
      return
      end