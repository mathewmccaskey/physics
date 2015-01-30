c-------------------------------------------------------------------------------------------------c
      subroutine rkf(derivs,n,y,x,h,yout,yerr)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This subroutine calculates a step in solving an ODE using the Runga-Kutta-Cash-Karp method.    c
c                                                                                                 c
c  derivs = external subroutine for calculating derivatives                                       c
c  n = number of ODEs being solved                                                                c
c  y = inital point y in the step                                                                 c
c  x = initial point x in the step                                                                c
c  h = stepsize                                                                                   c
c  yout = output of the new y-value                                                               c
c  yerr = estimated error of the step                                                             c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

c input parameters
      integer n
      double precision y(n),x(n),h,yout(n),yerr(n)
      external derivs

c parameters for this subroutine only
      integer i
      double precision k1(n),k2(n),k3(n),k4(n),k5(n),k6(n), ytemp(n)
      double precision c2,c3,c4,c5,c6
      double precision a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65
      double precision b1,b2,b3,b4,b5,b6
      double precision b1s,b2s,b3s,b3s,b4s,b5s,b6s
      parameter (c2=1.d0/5.d0,c3=3.d0/10.d0,c4=3.d0/5.d0,c5=1.d0,c6=7.d0/8.d0,a21=1.d0/5.d0,
     . a31=3.d0/40.d0,a32=9.d0/40.d0,a41=3.d0/10.d0,a42=-9.d0/10.d0,a43=6.d0/5.d0,
     . a51=-11.d0/54.d0,a52=5.d0/2.d0,a53=-70.d0/27.d0,a54=-35.d0/27.d0,
     . a61=1631.d0/55296.d0,a62=175.d0/512.d0,a63=575.d0/13824.d0,a64=44275.d0/110592.d0,
     . a65=253.d0/4096.d0,b1=37.d0/378.d0,b2=0.d0,b3=250.d0/621.d0,b4=125.d0/594.d0,
     . b5=0.d0,b6=512.d0/1771.d0,b1s=2825.d0/27648.d0,b2s=0.d0,b3s=18575.d0/48384.d0,
     . b4s=13525.d0/55296.d0,b5s=277.d0/14336.d0,b6s=1.d0/4.d0)

c many setup calculations are required
      call derivs(x,y,k1)
      do i=1,n
        ytemp(i) = y(i) + h*a21*k1(i)
      enddo
      call derivs(x+h*c2,ytemp,k2)
      do i=1,n
        ytemp(i) = y(i) + h*(a31*k1(i) + a32*k2(i))
      enddo
      call derivs(x+h*c3,ytemp,k3)
      do i=1,n
        ytemp(i) = y(i) + h*(a41*k1(i) + a42*k2(i) + a43*k3(i))
      enddo
      call derivs(x+h*c4,ytemp,k4)
      do i=1,n
        ytemp(i) = y(i) + h*(a51*k1(i) + a52*k2(i) + a53*k3(i) + a54*k4(i))
      enddo
      call derivs(x+h*c5,ytemp,k5)
      do i=1,n
        ytemp(i) = y(i) + h*(a61*k1(i) + a62*k2(i) + a63*k3(i) + a64*k4(i) + a65*k5(i))
      enddo
      call derivs(x+h*c6,ytemp,k6)

c time to calculate the result and error
      do i=1,n
        yout(i) = y(i) + h*(b1s*k1 + b2s*k2 + b3s*k3 + b4s*k4 + b5s*k5 + b6s*k6)
        yerr(i) = h*((b1-b1s)*k1 + (b2-b2s)*k2 + (b3-b3s)*k3 + (b4-b4s)*k4 + (b5-b5s)*k5
     .    + (b6-b6s)*k6)
      enddo

      return
      end