c-------------------------------------------------------------------------------------------------c
      subroutine rk4(y,dydx,n,x,h,yout,derivs)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Given values for the variables y(1:n) and their derivatives dydx(1:n) known at x, use the      c
c  fourth-order Runge-Kutta method to advance the solution over an interval h and return the      c
c  incremented variables as yout(1:n), which net not be a distinct array from y.  The user        c
c  supplies the subroutine derivs(x,y,dydx) which returns derivatives dydx at x.                  c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

c input parameters
      integer n
      double precision h, x, dydx(n), y(n), yout(n)
      external derivs

c parameters used in this subroutine only
      parameter (nmax = 50)
      integer i
      double precision h6, hh, xh, dym(nmax), dyt(nmax), yt(nmax)

      hh = h*0.5d0
      h6 = h/6.d0
      xh = x + hh
      do i=1,n
        yt(i) = y(i) + hh*dydx(i)
      enddo
      call derivs(xh,yt,dyt)
      do i=1,n
        yt(i) = y(i) + hh*dyt(i)
      enddo
      call derivs(xh,yt,dym)
      do i=1,n
        yt(i) = y(i) + h*dym(i)
        dym(i) = dyt(i) + dym(i)
      enddo
      call derivs(x+h,yt,dyt)
      do i=1,n
        yout(i) = y(i) + h6*(dydx(i)+dyt(i)+2.d0*dym(i))
      enddo

      return
      end