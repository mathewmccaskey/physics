c-------------------------------------------------------------------------------------------------c
      subroutine spline(x,y,n,yp1,ypn,y2)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e. y_i=f(x_i), with          c
c  x1<x2<..<xn, and given values yp1 and ypn for the first derivative of the interpolating        c
c  function at points 1 and n, respectively, this routine returns an array y2(1:n) of length n    c
c  which contains the second derivatives of the interpolating function at the tabulated points    c
c  xi.  If yp1 and/or ypn are equal to 1e30 or larger, the routine is signaled to set the         c
c  corresponding boundary condition for a natural spline, with zero second derivative on that     c
c  boundary.                                                                                      c
c                                                                                                 c
c  Parameter: NMAX is the largest anticipated value of n                                          c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

      integer n,nmax
      double precision yp1,ypn,x(n),y(n),y2(n)
      parameter (nmax=500)
      integer i,k
      double precision p,qn,sig,un,u(nmax)

c The lower boundary condition is set either to be "natural"
      if (yp1.gt.0.99d30) then
        y2(1) = 0.d0
        u(1) = 0.d0

c or else to have a specified first derivative
      else
        y2(1) = -0.5d0
        u(1) = (3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif

c This is the decomposition loop of the tridiagonal algorithm.
c y2 and u are used for temporary storage of the decomposed factors.
      do i=2,n-1
        sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p = sig*y2(i-1)+2.d0
        y2(i) = (sig-1.d0)/p
        u(i) = (6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     .    /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo

c The upper boundary condition is set either to be "natural"
      if (ypn.gt..99d30) then
        qn = 0.d0
        un = 0.d0

c or else to have a specified first derivative
      else
        qn = 0.5d0
        un = (3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.d0)

c This is the backsubstitution loop of the tridiagonal algorithm
      do k=n-1,1,-1
        y2(k) = y2(k)*y2(k+1)+u(k)
      enddo

      return
      end