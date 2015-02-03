ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Mathew McCaskey
c University of Wisconsin-Madison
c Spring 2007
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c---------------------------------------------------------------------
      subroutine relicdensity
c  Calculates the relic density of the singlet DM partilces
c
c---------------------------------------------------------------------
      include 'smsing.inc'
      include 'hdecay.inc'
	real*8 oh2last

c  calls the subroutine to calculate the relic density
c	print *,'mdm = ',mdm
	xfactor = 30d0

      call odeintegrator(mdm/xfactor)
	if(oh2.eq.'NAN') then 
		flag_fail = 1
		return
	endif
	xfactor = xfactor - 1d0
1234	oh2last = oh2
      call odeintegrator(mdm/xfactor)
	if(oh2.eq.'NAN') then 
		flag_fail = 1
		return
	endif
c	print *,oh2,oh2last
	if(dabs(oh2-oh2last)/oh2.gt.0.02d0) then
		xfactor = xfactor - 1d0
		goto 1234
	endif
	oh2unc = dabs(oh2-oh2last)/oh2
	

c  verify the relic denstiy with the approximate formula
c      call checkoh2

      return
      end


c-------------------------------------------------------------------
      subroutine odeintegrator(T1)
c Solves the ordinary differential equation to get relic density
c
c-------------------------------------------------------------------
      include 'smsing.inc'

c Using the template from numerical recipes
      integer nbad,nok,nvar
      real*8 Fo,T1,T2,eps,h1,hmin,qo

      parameter (eps=1.0d-3,nvar=1)
      external derivative, rkqsdp

c      T1=mdm/20d0      
c	if(mdm.lt.5d0) T1 = mdm/10d0
c	if(T1.lt.1d0) T1 = 1d0
      T2=1d-3        ! 1 MeV

      h1 = (T2-T1)/(1.0d4)                   !step size
      hmin = 0.0d0

c Initial value of the relic density function
      Fo = dlog(dsqrt((mdm/(2.d0*Pi*T1))**3))-mdm/T1

      call odeintdp(Fo,nvar,T1,T2,eps,h1,hmin,nok,nbad,
     *      derivative,rkqsdp)

      oh2 = dexp(Fo)*mdm*2.742d8

      print *,'Relic density for T1 = ',mdm/T1,' is ',oh2
c TAKE THIS OUT WHEN DONE!!!!!!!!
      stop

      return
      end


c-------------------------------------------------------------------
      subroutine derivative(T,F,dFdT)
c  This defines the derivative of the relic density at the given point T
c  So that we can use it with odeint and get the relic density
c
c------------------------------------------------------------------
      include 'smsing.inc'

      real*8 T,F,dFdT,no
      real*8 Mpl,gT,taacs
      parameter (Mpl=1.2209d19)

c First get the degrees of freedom based on temperature   
      gT = getgT(T)

c This calculates the thermally averaged annihilation cross section
c      call thermalave(T,taacs)
      taacs = thermalave(T)

      qo = dlog(dsqrt((mdm/(2.d0*Pi*T))**3)) - mdm/T

c differenial equation to find the relic density
      dFdT = taacs*Mpl*dsqrt(45.d0/(4.d0*Pi**3*gT))
     &   *dexp(F)*(1d0-dexp(2.d0*(qo-F)))

c This is to check the density we calculate with the equilibrium density
      no = dsqrt((mdm/(2.d0*Pi*T))**3)*dexp(-mdm/T)

      write(*,*) mdm/T,F,dexp(F),no,dfdt
c      write(1,*) mdm/T,F,dexp(F),no,dfdt
      return
      end


C---------------------------------------------------------------------
      double precision function getgT(T)
c
c This function gets the temperature dependent degrees of freedom
c
c---------------------------------------------------------------------
      include 'smsing.inc'

      real*8 T,gT

      if (T .gt. mt) then
        gT = 28.d0+7.d0/8.d0*90.d0
      else if (T .gt. mhsm) then
        gT = 28.d0+7.d0/8.d0*78.d0
      else if (T .gt. mz) then
        gT = 27.d0+7.d0/8.d0*78.d0
      else if (T .gt. mw) then
        gT = 24.d0+7.d0/8.d0*78.d0
      else if (T .gt. mb) then
        gT = 18.d0+7.d0/8.d0*78.d0                              
      else if (T .gt. mtau) then
        gT = 18.d0+7.d0/8.d0*66.d0
      else if (T .gt. mc) then
        gT = 18.d0+7.d0/8.d0*62.d0
      else if (T .gt. ms) then
        gT = 18.d0+7.d0/8.d0*50.d0
      else if (T .gt. mmu) then
        gT = 18.d0+7.d0/8.d0*38.d0
      else if (T .gt. md) then
        gT = 18.d0+7.d0/8.d0*34.d0
      else if (T .gt. mu) then
        gT = 18.d0+7.d0/8.d0*22.d0
      else if (T .gt. me) then
c By the time T gets here we will have already passed freezout temp
c so there is no need to exclude any more degrees of freedom
        gT = 18.d0+7.d0/8.d0*10.d0
      endif

      getgT = gT

      return
      end


c---------------------------------------------------------------------
      double precision function thermalave(T)
c This calculates the thermally averaged annihilation cross section
c
c---------------------------------------------------------------------
      include 'smsing.inc'

      real*8 T,taacs,K2,beta
      integer init,itmax,ncall,ndim,nprn
      real*8 avgi,chi2a,sd,xoff,region(4),eps
      real*8 fxn
      external fxn

      eps = 1.d-3

c We have a 2 dimensional integral over beta and T
      ndim = 2

c We are integrating beta from 0 to 1 in beta and epsilon around T
      region(1) = 0.d0
      region(1+ndim) = 1.d0
      region(2) = T
      region(2+ndim) = T+eps

c Playing around with the vegas driver these values give an accurate and 
c fast integral
      init = 0
      ncall = 5000
      itmax = 1

c This is set so that nothing is printed
      nprn = -1

c These are to find the value of the integral, error, and uncertainty
      avgi = 0.d0
      sd = 0.d0
      chi2a = 0.d0

      call vegasdp(region,ndim,fxn,init,ncall,itmax,nprn,avgi,sd,chi2a)

c      write(*,*) T,sd/avgi

c      taacs = avgi/eps
      thermalave = avgi/eps 

      return
      end


c-----------------------------------------------------------------------
      double precision function fxn(p,wgt)
c external function used for vegasdp
c
c-----------------------------------------------------------------------
      include 'smsing.inc'

      real*8 p(2),beta,wgt
      real*8 coeff,T,bessel1,bessel2
      real*8 ancross,s,besselratio

c p(1) = beta, p(2) = T
      beta = p(1)
      T = p(2)

c Center of mass energy
      s = 4.d0*mdm**2/(1-beta**2)

c      bessel2 = besskdp(2,mdm/T)

c Coefficient in front of the integral
c      coeff = 2.d0/(mdm*T*(bessel2**2.d0))
      coeff = 2.d0/(mdm*T)
      
c Call to calculate the coannihilation cross section given beta
      call annihilcross(ancross,beta)

c      bessel1 = bessk1dp(dsqrt(s)/T)

       besselratio = mdm*dsqrt(2d0/(pi*dsqrt(s)*T))*
     .     dexp((2d0*mdm-dsqrt(s))/T)*
     .          (1d0+3d0/(8d0*dsqrt(s)/T)+
     .               3d0*(-5d0)/(2d0*(8d0*dsqrt(s)/T)**2)+
     .                   105d0/(1024d0*(dsqrt(s)/T)**3))
     .  /(1d0+15d0/(8d0*mdm/T)+105d0/(128d0*(mdm/T)**2)-
     .                   315d0/(1024d0*(mdm/T)**3))**2

c	print *,dsqrt(s),mdm,t
c	print *,besselratio
c	print *,bessel1/bessel2**2
c	print *
c	pause

c      fxn = coeff*ancross*(beta**2.d0)*s/(1.d0-beta**2.d0)**(2.5d0)
c     &       *bessel1
      fxn = coeff*ancross*(beta**2.d0)*s/(1.d0-beta**2.d0)**(2.5d0)
     &       *besselratio

c	write(*,*)coeff,mdm/T,ancross,beta,s,bessel1

      return
      end


c---------------------------------------------------------------------
      subroutine annihilcross(ancross,beta)
c     Calculates the annihilation cross section for the singlet
c
c---------------------------------------------------------------------
      include 'smsing.inc'
      include 'hdecay.inc'
      
      real*8 ancross,s,z,a,b,beta

c Sets the annihilation cross section to be zero
      ancross = 0.d0

c Sets the center of mass energy 
      s = 4.d0*mdm**2/(1.d0-beta**2)

c This is a very common denomenator in all the cross sections
      z = (s-mhsm**2)**2+wdth(1)**2*mhsm**2

c These are needed for SS->HH due to the t and u channels
      a = mhsm**2-s/2.d0
      b = s/2.d0*dsqrt(1-mhsm**2/mdm**2*(1-beta**2))

c Contribution from SS->hh and SS->h->hh (s channel)
c and SS->S->HH (t and u channel)
      if (mdm .gt. mhsm) then
        ancross = ancross + (2.d0*del2)**2/(64.d0*Pi*s)
     &    *dsqrt(1.d0-4.d0*mhsm**2/s)
     &    *(1.d0+(6.d0*lambda)**2*vsm**2/(16.d0*z)
     &          +(6.d0*lambda)*vsm**2*(s-mhsm**2)/(2.d0*z)
     &          +(2.d0*del2)**2*vsm**4/4.d0*(2.d0/(a**2-b**2)+1.d0/(a*b)
     &               *dlog((a-b)/(a+b)))
     &          +((2.d0*del2)*vsm**2+(2.d0*del2)*(6.d0*lambda)*vsm**4
     &               *(s-mhsm**2)/(4.d0*z))*dlog((a+b)/(a-b))/b)

c        ancross = ancross + del2**2/(64.d0*Pi*s)
c     &     *dsqrt(1.d0-4.d0*mhsm**2/s)
c     &     *(1.d0+lambda**2*vsm**2/(16.d0*z)
c     &           +lambda*vsm**2*(s-mhsm**2)/(2.d0*z)
c     &           +del2**2*vsm**4/4.d0*(2.d0/(a**2-b**2)+1.d0/(a*b)
c     &                *dlog((a-b)/(a+b)))
c     &           +(del2*vsm**2+del2*lambda*vsm**4*(s-mhsm**2)/(4.d0*z))
c     &                *dlog((a+b)/(a-b))/b)

c This is the old annihilation cross section that does not include
c the t and u channel singlet exchange
c        ancross = ancross + del2**2/4.d0*(1.d0+lambda*vsm**2/(2*((
c     &    s-mhsm**2)**2+wdth(1)**2*mhsm**2))*(lambda*vsm**2/2.d0+
c     &    s-mhsm**2))*(1.d0-4.d0*mhsm**2/s)**(0.5d0)/(16.d0*Pi*s)
      endif

c Contribution from SS->W+W- 
      if (mdm .gt. mw) then
        ancross = ancross+(2.d0+(1.d0-s/(2.d0*mw**2))**2)*(2.d0*del2)**2
     &    *mw**4/(16.d0*Pi*s*z)*(1.d0-4.d0*mw**2/s)**(0.5d0)
c        ancross = ancross+(2.d0+(1.d0-s/(2.d0*mw**2))**2)*del2**2
c     &    *mw**4/(16.d0*Pi*s*z)*(1.d0-4.d0*mw**2/s)**(0.5d0)
      endif

c Contribution from SS->ZZ
      if (mdm .gt. mz) then
        ancross = ancross+(2.d0+(1.d0-s/(2.d0*mz**2))**2)*(2.d0*del2)**2
     &    *mz**4/(32.d0*Pi*s*z)*(1.d0-4.d0*mz**2/s)**(0.5d0)
c        ancross = ancross + (2.d0+(1.d0-s/(2.d0*mz**2))**2)*del2**2
c     &    *mz**4/(32.d0*Pi*s*z)*(1.d0-4.d0*mz**2/s)**(0.5d0)
      endif

c Contributions from SS->fbar f
c SS->e+e-
      if (mdm .gt. me) then
        ancross = ancross + 
     &    (2.d0*del2)**2*vsm**2*le**2
     &    *(1.d0-4.d0*me**2/s)**(1.5d0)/(32.d0*Pi*z)
c        ancross = ancross + 
c     &  del2**2*vsm**2*le**2*(1.d0-4.d0*me**2/s)**(1.5d0)/(32.d0*Pi*z)
      endif

c SS->mu+mu-
      if (mdm .gt. mmu) then
        ancross = ancross + 
     &    (2.d0*del2)**2*vsm**2*lmu**2
     &    *(1.d0-4.d0*mmu**2/s)**(1.5d0)/(32.d0*Pi*z)
c        ancross = ancross + 
c     &  del2**2*vsm**2*lmu**2*(1.d0-4.d0*mmu**2/s)**(1.5d0)/(32.d0*Pi*z)
      endif

c SS->tau+tau-
      if (mdm .gt. mtau) then
        ancross = ancross + 
     &  (2.d0*del2)**2*vsm**2*ltau**2*(1.d0-4.d0*mtau**2/s)**(1.5d0)
     &    /(32.d0*Pi*z)
c        ancross = ancross + 
c     &  del2**2*vsm**2*ltau**2*(1.d0-4.d0*mtau**2/s)**(1.5d0)
c     &    /(32.d0*Pi*z)
      endif

c SS->u ubar
      if (mdm .gt. mu) then
        ancross = ancross + 
     &    (2.d0*del2)**2*vsm**2*lu**2
     &    *(1.d0-4.d0*mu**2/s)**(1.5d0)/(32.d0*Pi*z)
c        ancross = ancross + 
c     &  del2**2*vsm**2*lu**2*(1.d0-4.d0*mu**2/s)**(1.5d0)/(32.d0*Pi*z)
      endif

c SS->d dbar         
      if (mdm .gt. md) then
        ancross = ancross + 
     &    (2.d0*del2)**2*vsm**2*ld**2
     &    *(1.d0-4.d0*md**2/s)**(1.5d0)/(32.d0*Pi*z)
c        ancross = ancross + 
c     &  del2**2*vsm**2*ld**2*(1.d0-4.d0*md**2/s)**(1.5d0)/(32.d0*Pi*z)
      endif

c SS->s sbar             
      if (mdm .gt. ms) then
        ancross = ancross + 
     &    (2.d0*del2)**2*vsm**2*ls**2
     &    *(1.d0-4.d0*ms**2/s)**(1.5d0)/(32.d0*Pi*z)
c        ancross = ancross + 
c     &  del2**2*vsm**2*ls**2*(1.d0-4.d0*ms**2/s)**(1.5d0)/(32.d0*Pi*z)
      endif

c SS->c cbar
      if (mdm .gt. mc) then
        ancross = ancross + 
     &    (2.d0*del2)**2*vsm**2*lc**2
     &    *(1.d0-4.d0*mc**2/s)**(1.5d0)/(32.d0*Pi*z)
c        ancross = ancross + 
c     &  del2**2*vsm**2*lc**2*(1.d0-4.d0*mc**2/s)**(1.5d0)/(32.d0*Pi*z)
      endif

c SS->b bbar
      if (mdm .gt. mb) then
        ancross = ancross + 
     &    (2.d0*del2)**2*vsm**2*lb**2
     &    *(1.d0-4.d0*mb**2/s)**(1.5d0)/(32.d0*Pi*z)
c        ancross = ancross + 
c     &  del2**2*vsm**2*lb**2*(1.d0-4.d0*mb**2/s)**(1.5d0)/(32.d0*Pi*z)
      endif

c SS->t tbar
      if (mdm .gt. mt) then
        ancross = ancross + 
     &    (2.d0*del2)**2*vsm**2*lt**2
     &    *(1.d0-4.d0*mt**2/s)**(1.5d0)/(32.d0*Pi*z)
c        ancross = ancross + 
c     &  del2**2*vsm**2*lt**2*(1.d0-4.d0*mt**2/s)**(1.5d0)/(32.d0*Pi*z)
      endif

      return
      end

c---------------------------------------------------------------------
      subroutine checkoh2
c Calculates the approximate value for the relic density and 
c compares it with the value calculated value from odeint
c
c---------------------------------------------------------------------
      include 'smsing.inc'
      real*8 x1,x2,tol,nwave,gT,T,cross
      real*8 getgT,zbrentdp,thermalave,func
      integer i
      external func

      x1 = 10.d0
      x2 = 100.d0

      tol = 0.01d0
      xf = 20d0!zbrentdp(func,x1,x2,tol)

      T = mdm/xf
      gT = getgT(T)
      cross = thermalave(T)

c nwave = 0 - s wave scattering
c nwave = 1 - p wave scattering
      nwave = 1.d0

      oh2approx = 1.07d9*xf**(1d0+nwave)*(1d0+nwave)/
     &              (dsqrt(gT)*1.22d19*cross)

c      write(*,*) oh2,oh2approx,oh2/oh2approx,xf,mdm

c      write(*,*) 'Relic density',oh2,oh2approx,xf,mdm

      return
      end


c---------------------------------------------------------------------
      double precision function func(x)
c external function used for zbrentdp
c
c---------------------------------------------------------------------
      include 'smsing.inc'

      real*8 x,T,cross,nwave
      real*8 getgT,thermalave

      T = mdm/x
      gT = getgT(T)

      cross = thermalave(T)
      nwave = 1.d0

c      func = x-dlog(0.038d0*dsqrt(gT*x)*1.2209d19*mdm*cross)

      func = x-dlog(0.038d0*(1.d0+nwave)*dsqrt(gT*x**(1.d0+nwave))
     &        *1.2209d19*mdm*cross)

c      xf = dlog(0.038d0*(1d0+nwave)*dsqrt(gT)*1.22d19*mdm*cross)
c      xf = xf - (.5d0+nwave)*dlog(xf)

      return
      end
