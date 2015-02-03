c-------------------------------------------------------------------------c
c 1-D function integration
c-------------------------------------------------------------------------c
      subroutine romberg(func, a, b, integral, eps, iter)
c-------------------------------------------------------------------------c
c                                                                         c
c  Returns the integral of the function func from a to b.  Integration    c
c  is performed by Romberg's method of order 2K, where, e.g. K=2 is       c
c  Simpson's rule.  Based on numerical recipe's qromb.                    c
c                                                                         c
c  Parameters:                                                            c
c  func - the user supplied double precision function                     c
c  a - lower bound of the integral                                        c
c  b - upper bound of the integral                                        c
c  integral - result returned by romberg                                  c
c  eps - the fractional accuracy desired as determined by the             c
c        extrapolation error estimate                                     c
c  iter - the minimum number of iterations (to improve accuracy)          c
c                                                                         c
c  jmax - limits the total number of steps                                c
c  k - the number of points used in the extrapolation.                    c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      
c input parameters
      double precision func, a, b, integral, eps
      integer iter
      external func

c parameters for this subroutine only
      integer jmax, jmaxp, k, km
      parameter (jmax=50, jmaxp=jmax+1, k=5, km=k-1)

c uses polint, trapzd
      integer j

c these store the successive trapezoidal approximations and their
c relative stepsizes
      double precision dss, h(jmaxp), s(jmaxp)

      h(1) = 1.d0
      do j=1,jmax
        call trapzd(func,a,b,s(j),j)
        if ((j.ge.k).and.(j.ge.iter)) then
          call polint(h(j-km),s(j-km),k,0.d0,integral,dss)
          if (dabs(dss).le.eps*dabs(integral)) return
        endif
        s(j+1)=s(j)
        
c This is a key step: The factor is 0.25 even though the stepsize decrease is only 0.5
c This makes the extrapolation a polynomial in h^2.
        h(j+1)=0.25d0*h(j)
      enddo
      
c      pause 'too many steps in romberg'
      end



c-------------------------------------------------------------------------c
      subroutine trapzd(func, a, b, s, n)
c-------------------------------------------------------------------------c
c                                                                         c
c  This routine computes the nth stage of refinement of an extended       c
c  trapezoidal rule.  func is input as the name of the function to be     c
c  integrated between limits a and b, also input.  When called with n=1,  c
c  the routine returns as s the crudent estimate of the integral.         c
c  Subsequent calls with n=2,3,... will improve the accuracy of s by      c
c  adding 2^(n-2) additional interior points.  s should not be modified   c
c  between sequential calls.                                              c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none

      integer n
      double precision a, b, s, func
      external func
      integer it, j
      double precision del, summ, tnm, x

      if (n.eq.1) then
        s = 0.5d0*(b-a)*(func(a)+func(b))
      else
        it = 2**(n-2)
        tnm = it

c This is the spacing of the points to be added
        del = (b-a)/tnm
        x = a+0.5d0*del
        summ = 0.d0
        do j=1,it
          summ = summ+func(x)
          x = x+del
        enddo

c This replaces s by its refined value
        s = 0.5d0*(s+(b-a)*summ/tnm)
      endif
      return
      end



c-------------------------------------------------------------------------c
      subroutine polint(xa, ya, n, x, y, dy)
c-------------------------------------------------------------------------c
c                                                                         c
c  Given arrays xa and ya, each of length n, and given a value x, this    c
c  routine returns a value y, and an error estimate dy.  If P(x) is the   c
c  polynomial of degree N=1 such that P(xa_i) = ya_i, i=1...n then the    c
c  returned value y = P(x).                                               c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none

      integer n, nmax
      double precision dy, x, y, xa(n), ya(n)
      parameter (nmax=10)
      integer i, m, ns
      double precision den, dif, dift, ho, hp, w, c(nmax), d(nmax)

      ns = 1
      dif = dabs(x-xa(1))

c Here we find the index ns of the closest table entry
      do i=1,n
        dift = dabs(x-xa(i))
        if (dift.lt.dif) then
          ns = i
          dif = dift
        endif

c and initialize the tableau of c's and d's
        c(i) = ya(i)
        d(i) = ya(i)
      enddo

c This is the initial approximation to y
      y = ya(ns)
      ns = ns-1

c For each column of the tableau,
c we loop over the current c's and d's and update them
      do m=1,n-1
        do i=1,n-m
          ho = xa(i)-x
          hp = xa(i+m)-x
          w = c(i+1)-d(i)
          den = ho-hp

c This error can occur only if two input xa's are (to within roundoff) identical.
c          if (den .eq. 0.d0) pause 'failure in polint'
          den = w/den

c Here the c's and d's are updated
          d(i) = hp*den
          c(i) = ho*den
        enddo
        if(2*ns .lt. n-m) then
          dy = c(ns+1)
        else
          dy = d(ns)
          ns = ns-1
        endif
        y = y+dy
      enddo
      return
      end