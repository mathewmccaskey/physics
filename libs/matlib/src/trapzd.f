c-------------------------------------------------------------------------------------------------c
      subroutine trapzd(func,a,b,s,n)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This routine computes the nth stage of refinement of an extended trapezoidal rule.  func is    c
c  input as the name of the function to be integrated between limits a and b, also input.  When   c
c  called with n=1, the routine returns as s the crudent estimate of the integral.  Subsequent    c
c  calls with n=2,3,... will improve the accuracy of s by adding 2^(n-2) additional interior      c
c  points.  s should not be modified between sequential calls.                                    c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

      integer n
      double precision a,b,s,func
      external func
      integer it,j
      double precision del,sum,tnm,x

      if (n.eq.1) then
        s = 0.5d0*(b-a)*(func(a)+func(b))
      else
        it = 2**(n-2)
        tnm = it

c This is the spacing of the points to be added
        del = (b-a)/tnm
        x = a+0.5d0*del
        sum = 0.d0
        do j=1,it
          sum = sum+func(x)
          x = x+del
        enddo

c This replaces s by its refined value
        s = 0.5d0*(s+(b-a)*sum/tnm)
      endif

      return
      end