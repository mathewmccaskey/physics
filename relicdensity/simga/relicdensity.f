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
      double precision oh2, n_now, oh2_approx
      double precision Y_now, x_F, s, msq, fxn

c setting up the boundaries and set sizes
      x_1 = 15.d0
      x_2 = 1000.d0
      h1 = 0.1d0
      hmin = 0.0d0

c setting up the DM paramters
      g = 1.d0
      m = 1000.d0

c setting the number of points in the thermally averaged cross setion array
      npoints = 21

c initialize the relativistic degrees of freedom
      call readinfile('datafiles/g_star.txt',100,tempgstar,gstar)
      call readinfile('datafiles/g_star_S.txt',100,tempgstarS,gstarS)
      call readinfile('datafiles/taacs1.txt',npoints,xtaacs,taacs)
      g_star_S = lineint(m/x_1,tempgstarS,gstarS,100)

c Initial value of the number density per entropy density, Y
c      Yo = 45.d0/(4.d0*pi**4)*(g/g_star_S)*(pi/2.d0)**(0.5d0)*x_1**(1.5d0)*dexp(-x_1)
      Yo = 45.d0/(4.d0*pi**4)*(g/g_star_S)*(x_1**2+0.d0)*bessk(2,dsqrt(x_1**2+0.d0))

      x_pass = x_1

c  So this integral is very difficult to do numerically
      call odeint(Yo,nvar,x_1,x_2,eps,h1,hmin,nok,nbad,derivative,rkqs)

      oh2 = Yo*2889.2d0*m/(1.05d-5)

      write(*,*) 'Relic density from the chemical rate equation ',oh2


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
      double precision s, H_m, Y_eq, cs

      g_star = lineint(m/x_i,tempgstar,gstar,100)
      g_star_S = lineint(m/x_i,tempgstarS,gstarS,100)

      s = 2.d0*pi**2/45.d0*g_star_S*(m/x_i)**3
      H_m = 0.602d0*mpl/(g_star**(0.5)*m**2)
      
c      Y_eq = 45.d0/(4.d0*pi**4)*(g/g_star_S)*x_i**2*bessk(2,x_i)
c      Y_eq = 45.d0/(4.d0*pi**4)*(g/g_star_S)*x_i**2*(pi/(2.d0*x_i))**(0.5)*dexp(-x_i)
      Y_eq = 45.d0/(4.d0*pi**4)*(g/g_star_S)*(x_i**2+0.d0)*bessk(2,dsqrt(x_i**2+0.d0))


      cs = 2.d0*lineint(x_i,xtaacs,taacs,npoints)

c      dYdx = -x_i*cross*s*H_m*(Y_i-Y_eq)*(Y_i+Y_eq)
c      dYdx = -x_i*cs*s*H_m*(Y_i-Y_eq)*(Y_i+Y_eq)
      
      dYdx = -dsqrt(pi/(45*g_star))*g_star_S*mpl*m/x_i**2*cs*(Y_i-Y_eq)*(Y_i+Y_eq)
      
c this is for testing purposes
c      write(*,fmt='(4(E14.8,1x))') x_i,Y_i,Y_eq,dYdx

      return
      end