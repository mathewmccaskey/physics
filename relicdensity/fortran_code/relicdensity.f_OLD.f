c-----------------------------------------------------------c
      program relicdensity
c-----------------------------------------------------------c
      implicit none
      include 'relic.inc'
      include 'matlib.inc'
      include 'nrvars.inc'

c parameters for the odeintegrator
      integer nok,nbad,nvar
      double precision Yo, x_1, x_2, eps, h1, hmin

      parameter (eps=1.d-3,nvar=1)
      external derivative, rkqs

c DM parameters
      double precision oh2, n_now, oh2_approx, taacs
      double precision Y_now, x_F, s, msq

c setting up the boundaries and set sizes
      x_1 = 20.d0
      x_2 = 1000.d0
      h1 = 0.1d0
      hmin = 0.0d0

c setting up the DM paramters
      g = 1.d0
      m = 1000.d0

c initialize the relativistic degrees of freedom
      call readinfile('datafiles/g_star.txt',100,tempgstar,gstar)
      call readinfile('datafiles/g_star_S.txt',100,tempgstarS,gstarS)
      g_star_S = getvalue(m/x_1,tempgstarS,gstarS,100)

c Initial value of the number density per entropy density, Y
      Yo = 45.d0/(4.d0*pi**4)*(g/g_star_S)*(pi/2.d0)**(0.5d0)*x_1**(1.5d0)*dexp(-x_1)

c Thermally averaged annihilation cross section check

c      do i=1, 100
c       write(*,*) 10*i, taacs(dble(10.d0*i))
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

      write(*,*) 'Relic density from boltzmann equation ',oh2
      write(*,*) 'The freezeout temperature is ',x_F
      write(*,*) 'Relic density approximated to be ',oh2_approx

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
      include 'nrvars.inc'

c input parameters
      double precision x_i, Y_i, dYdx

c parameters in this subroutine only
      double precision s, H_m, Y_eq, cross, taacs

      g_star = getvalue(m/x_i,tempgstar,gstar,100)
      g_star_S = getvalue(m/x_i,tempgstarS,gstarS,100)

      s = 2.d0*pi**2/45.d0*g_star_S*(m/x_i)**3
      H_m = 0.602d0*mpl/(g_star**(0.5)*m**2)
      Y_eq = 45.d0/(4.d0*pi**4)*(g/g_star_S)*(pi/2.d0)**(0.5d0)*x_i**(1.5d0)*dexp(-x_i)
      cross = taacs(x_i)

      dYdx = -x_i*cross*s*H_m*(Y_i-Y_eq)*(Y_i+Y_eq)

c this is for testing purposes
      write(*,fmt='(4(E14.8,1x))') x_i,Y_i,Y_eq,dYdx

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
      include 'nrvars.inc'

c input parameters
      double precision x, taacs

c parameters for the vegas integration
      integer init, itmax, ncall, ndim, nprn
      double precision avgi, chi2a, sd, xoff, region(6), eps
      double precision fxn
      external fxn

      eps = 1.d-3

c We have a 3 dimensional integral over E1, E2, and cos(theta_12)
      ndim = 3

c Velocities are going from 0 to 1 and cos(theta_12) is going from -1 to 1
      region(1) = 0.d0
      region(1+ndim) = 1.d0
      region(2) = 0.d0
      region(2+ndim) = 1.d0
      region(3) = -1.d0
      region(3+ndim) = 1.d0

c settings to get a fast and accurate integral
      init = 0
      ncall = 10000
      itmax = 5

c This is set so that nothing gets printed
      nprn = -1

c These are to find the value of the integral, error, and uncertainty
      avgi =  0.d0
      sd = 0.d0
      chi2a = 0.d0

      x_pass = x

      call vegas(region,ndim,fxn,init,ncall,itmax,nprn,avgi,sd,chi2a)

      init = 1
      ncall = 100000
      itmax = 1

      call vegas(region,ndim,fxn,init,ncall,itmax,nprn,avgi,sd,chi2a)

      taacs = avgi
      taacs = 1d-8

      return
      end



c-----------------------------------------------------------c
      function fxn(var,wgt)
c-----------------------------------------------------------c
c                                                           c
c external function being integrated over using vegas.      c
c                                                           c
c-----------------------------------------------------------c
      implicit none
      include 'relic.inc'
      include 'matlib.inc'
      include 'nrvars.inc'

c input parameters
      double precision var(3),wgt

c variables used in this function only
      double precision fxn, b1, b2, cos12, n_eq, s, msq, phase
      double precision mx, my, gamy      

c setting dependent variables
      b1 = var(1)
      b2 = var(2)
      cos12 = var(3)

      temp = m/x_pass

      n_eq = g*(m**2/(2.d0*pi*x_pass))**(1.5d0)

      s = 2.d0*m**2*(1.d0+(1.d0-b1**2)**(-0.5)*(1.d0-b2**2)**(-0.5)*(1.d0-b1*b2*cos12))

      mx = 200.d0
      my = 500.d0
      g = 0.1d0
      gamy = 17.67d0
c      gamy = 8.d0/pi*4.d0*pi/g**2*(mx/730.d0)**3*(1.d0+mx/(2.d0*my))*(1.d0-my/(2.d0*mx)**2
      msq = g**2*my*gamy*4.d0*mx**2/((s-my**2)**2+my**2*gamy**2)/4.d0

c      msq = 100000.d0/(s-50.d0**2)**2/(8.d0*pi)*(1.d0-4.d0*50.d0**2/s)**(0.5)

      fxn = 2*pi**2*g**2/((2.d0*pi)**6*n_eq**2)*
     .  m**2*b1**2/(1.d0-b1**2)**2 * m**2*b2**2/(1.d0-b2**2)**2*
     .  dexp(x_pass*(1.d0-1.d0/dsqrt(1.d0-b1**2)))*
     .  dexp(x_pass*(1.d0-1.d0/dsqrt(1.d0-b2**2)))*msq

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
      include 'nrvars.inc'

c input parameters
      double precision xf, yinf

c Parameters used in this subroutine only
      double precision xf_old, cross, s, msq
      integer n

      n = 0
      xf = 100.d0
      xf_old = 10.d0

      s = 4.d0*m**2
      msq = 100000.d0/(s-50.d0**2)**2/(8.d0*pi)*(1.d0-4.d0*50.d0**2/s)**(0.5)
      cross = 1.d0/s*msq

      do 100 while (dabs(xf-xf_old) .gt. 0.00001d0)
        xf_old = xf

        g_star_S = getvalue(m/xf,tempgstarS,gstarS,100)
        g_star = getvalue(m/xf,tempgstar,gstar,100)

        xf = dlog(0.038d0*dble(n+1)*g_star_S/g_star**(0.5)*mpl*m*cross)
     .   -(dble(n)+0.5d0)*dlog(dlog(0.038d0*dble(n+1)*g_star_S/g_star**(0.5)*mpl*m*cross))

 100  continue

      yinf = 3.79d0*dble(n+1)*g_star**(0.5)/g_star_S*xf/m/mpl/cross

      return
      end
