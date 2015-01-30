c-------------------------------------------------------------------------------------------------c
      subroutine qromb(func,a,b,ss)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Returns as ss the integral of the function func from a to b.  Integration is performed by      c
c  Romberg's method of order 2K, where, e.g. K=2 is Simpson's rule.  Parameters: EPS is the       c
c  fractional accuracy desired as determined by the extrapolation error estimate; JMAX limits     c
c  the total number of steps; K is the number of poitnes used in the extrapolation                c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

      integer jmax, jmaxp, k, km
      double precision a, b, func, ss, eps
      external func
      parameter (eps=1.0d-6, jmax=20, jmaxp=jmax+1, k=5, km=k-1)

c uses polint, trapzd
      integer j

c these store the successive trapezoidal approximations and their
c relative stepsizes
      double precision dss, h(jmaxp), s(jmaxp)
      h(1) = 1.d0
      do j=1,jmax
        call trapzd(func,a,b,s(j),j)
        if (j.ge.k) then
          call polint(h(j-km),s(j-km),k,0.d0,ss,dss)
          if (dabs(dss).le.(eps*dabs(ss))) return
        endif
        s(j+1)=s(j)

c This is a key step: The factor is 0.25 even though the stepsize decrease is only 0.5
c This makes the extrapolation a polynomial in h^2.
        h(j+1)=0.25d0*h(j)
      enddo

      pause 'too many steps in qromb'
      end