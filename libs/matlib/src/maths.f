c-------------------------------------------------------------------------------------------------c
      function dabscmplx(c)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Returns the absolute value of the complex number c.                                            c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

c function definition and input parameters
      double precision dabscmplx
      double complex c

      dabscmplx = dsqrt(dble(c*conjg(c)))

      return
      end
      
      
      
c-------------------------------------------------------------------------------------------------c
      function factorial(n)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This function calculates the factorial for the user supplied integer value n.                  c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      
c function definition and input parameters
      integer factorial, n

c parameters for this function only
      integer i
      
      factorial = 1
      
      if (n.eq.0) then
        return
      else if (n.lt.0) then
        pause 'trying to calculate the factorial of a negative number'
        stop
      else
        do i=2,n
          factorial = factorial*i
        enddo
      endif
      
      return
      end