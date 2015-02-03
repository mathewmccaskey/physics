c-----------------------------------------------------------c
      program relicdensity
c-----------------------------------------------------------c
      implicit none
      include 'relic.inc'
      include 'matlib.inc'

c parameters for the odeintegrator
      integer nok,nbad,nvar
      double precision Yo, x_1, x_2, eps, h1, hmin

      parameter (eps=1.d-5,nvar=1)
      external derivative, rkqs

c DM parameters
      double precision oh2, n_now, oh2_approx, taacs
      double precision Y_now, x_F, s, msq, fxn

c setting up the boundaries and set sizes
      x_1 = 15.d0
      x_2 = 1000.d0
      h1 = 0.1d0
      hmin = 0.0d0

c setting up the DM paramters
      g = 1.d0
      m = 200.d0

c initialize the relativistic degrees of freedom
      call readinfile('datafiles/g_star.txt',100,tempgstar,gstar)
      call readinfile('datafiles/g_star_S.txt',100,tempgstarS,gstarS)
      g_star_S = getvalue(m/x_1,tempgstarS,gstarS,100)

c Initial value of the number density per entropy density, Y
      Yo = 45.d0/(4.d0*pi**4)*(g/g_star_S)*(pi/2.d0)**(0.5d0)*x_1**(1.5d0)*dexp(-x_1)

      x_pass = 15.d0

c Thermally averaged annihilation cross section check
c      do i=1, 1000
c         write(*,*) i, taacs(dble(i))
c        write(*,*) dble(i)/101.d0, fxn(dble(i)/101.d0)
c      enddo

c      s = 4.d0*m**2
c      msq = 10000.d0/(s-50.d0**2)**2/(8.d0*pi)*(1.d0-4.d0*50.d0**2/s)**(0.5)

c      write(*,*) 'taacs_today = ',1/s*msq

c  So this integral is very difficult to do numerically
      call odeint(Yo,nvar,x_1,x_2,eps,h1,hmin,nok,nbad,derivative,rkqs)

      oh2 = Yo*2889.2d0*m/(1.05d-5)

c Instead we will use the approximation solution outlined in Kolb and TUrner
      call approxsolution(x_F,Y_now)

      n_now = Y_now*2889.2d0
      oh2_approx = n_now*m/(1.05d-5)

      write(*,*) 'Relic density from the chemical rate equation ',oh2
c      write(*,*) 'The freezeout temperature is ',x_F
c      write(*,*) 'Relic density approximated to be ',oh2_approx

      end



c-----------------------------------------------------------c
      subroutine derivative(x_i,Y_i,dYdx)
c-----------------------------------------------------------c
c                                                           c
c  This subroutine defines the derivative of Y with respect c
c  to x at any point x along the ODE integrator.            c
c                                                           c
c-----------------------------------------------------------c
      implicit none
      include 'relic.inc'
      include 'matlib.inc'

c input parameters
      double precision x_i, Y_i, dYdx

c parameters in this subroutine only
      double precision s, H_m, Y_eq, cs, taacs

      g_star = getvalue(m/x_i,tempgstar,gstar,100)
      g_star_S = getvalue(m/x_i,tempgstarS,gstarS,100)

      s = 2.d0*pi**2/45.d0*g_star_S*(m/x_i)**3
      H_m = 0.602d0*mpl/(g_star**(0.5)*m**2)
      
c      Y_eq = 45.d0/(4.d0*pi**4)*(g/g_star_S)*x_i**2*bessk(2,x_i)

      Y_eq = 45.d0/(4.d0*pi**4)*(g/g_star_S)*x_i**2*(pi/(2.d0*x_i))**(0.5)*dexp(-x_i)

      write(*,fmt='(5(E15.8,1x))') pi, x_i, g, g_star_S, Y_eq


      cs = taacs(x_i)

c      dYdx = -x_i*cross*s*H_m*(Y_i-Y_eq)*(Y_i+Y_eq)
      dYdx = -x_i*cs*s*H_m*(Y_i-Y_eq)*(Y_i+Y_eq)
      
      dYdx = -dsqrt(pi/(45*g_star))*g_star_S*mpl*m/x_i**2*cs*(Y_i-Y_eq)*(Y_i+Y_eq)
      
c this is for testing purposes
c      write(*,fmt='(4(E14.8,1x))') x_i,Y_i,Y_eq,dYdx

      return
      end



c-----------------------------------------------------------c
      function taacs(x)
c-----------------------------------------------------------c
c                                                           c
c This function calculates the thermally averaged           c
c  annihilation cross section.  The derivation was based    c
c  off the discussion in section 5.2 of Kolb and Turner.    c
c                                                           c
c-----------------------------------------------------------c
      implicit none
      include 'relic.inc'
      include 'matlib.inc'

c input parameters
      double precision x, taacs

c external function to be integrated over
      double precision fxn
      external fxn

c passing x so that the function can use it
      x_pass = x

c the last argument calls the function over 10,000 times
      call qromb(fxn,0.d0,1.d0-1.d-10,taacs)

      taacs = 1d-10

      return
      end



c-----------------------------------------------------------c
      function fxn(beta)
c-----------------------------------------------------------c
c                                                           c
c external function being integrated over using vegas.      c
c                                                           c
c-----------------------------------------------------------c
      implicit none
      include 'relic.inc'
      include 'matlib.inc'

c input parameters
      double precision beta

c variables used in this function only
      double precision fxn, s, msq
      double precision mx, my, gamy      

c setting dependent variables
      s = 4.d0*m**2/(1.d0-beta**2)

      temp = m/x_pass

c      fxn = x_pass/(m**2*bessk(2,x_pass)**2)*beta**2/(1.d0-beta**2)**(2.5d0)*
c     .   bessk1(s**0.5/temp)*msq

      mx = 200.d0
      my = 500.d0
      g = 1.d0
      gamy = 17.67d0
      msq = g**2*my*gamy*4.d0*mx**2/((s-my**2)**2+my**2*gamy**2)/4.d0

      fxn = x_pass**(1.5d0)/(m**2*pi**(0.5d0))*beta**2/(1.d0-beta**2)**(2.25d0)*msq*
     .  dexp(2.d0*x_pass*(1.d0-1.d0/(1.d0-beta**2)**(0.5d0)))

      return
      end



c-----------------------------------------------------------c
      subroutine approxsolution(xf,yinf)
c-----------------------------------------------------------c
c                                                           c
c  This subroutine uses the approximate solution worked out c
c  in Kolb and Turner.                                      c
c                                                           c
c-----------------------------------------------------------c
      implicit none
      include 'relic.inc'
      include 'matlib.inc'

c input parameters
      double precision xf, yinf

c Parameters used in this subroutine only
      double precision xf_old, cs, s, msq
      double precision mx, gamy, my
      integer n

      n = 0
      xf = 100.d0
      xf_old = 10.d0

      s = 4.d0*m**2
      my = 500.d0
      g = 1.d0
      gamy = 17.67d0
      msq = g**2*my*gamy*4.d0*m**2/((s-my**2)**2+my**2*gamy**2)/4.d0
      cs = 1.d0/s*msq

      do 100 while (dabs(xf-xf_old) .gt. 0.00001d0)
        xf_old = xf

        g_star_S = getvalue(m/xf,tempgstarS,gstarS,100)
        g_star = getvalue(m/xf,tempgstar,gstar,100)

        xf = dlog(0.038d0*dble(n+1)*g_star_S/g_star**(0.5)*mpl*m*cs)
     .   -(dble(n)+0.5d0)*dlog(dlog(0.038d0*dble(n+1)*g_star_S/g_star**(0.5)*mpl*m*cs))

 100  continue

      xf = 18.35d0

      yinf = 3.79d0*dble(n+1)*g_star**(0.5)/g_star_S*xf/m/mpl/cs

      return
      end
