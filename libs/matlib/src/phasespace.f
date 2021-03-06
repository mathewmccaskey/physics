c-------------------------------------------------------------------------------------------------c
      subroutine phasespace2(pinit,mf1,mf2,pf1,pf2,wgt,x,vegasint)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Calculates the two-body decay of an initial state with momentum pinit into two final state     c
c  particles of mass mf1 and mf2 with momenta pf1 and pf2.  The phase space weight is calculated  c
c  for the randomly generated momenta of the final state particles                                c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

c input parameters
      double precision pinit(0:3), mf1, mf2, wgt, x(2)
      integer vegasint

c parameters appearing in this subroutine only
      double precision ctheta, root_s, piboost(0:3), pfabs,
     . stheta

c setting the rest energy of the initial state and the final particle mometum
      root_s = invmass(pinit)

      if (root_s.lt.(mf1+mf2)) then
       write(*,*) 'Not a valid two-body decay'
       stop
      endif

      pfabs = dsqrt(lambda(root_s**2,mf1**2,mf2**2))/(2.d0*root_s)

c if vegas is used to integrate then use the input array x
      if (vegasint.eq.1) then
       ctheta = 2.d0*x(1)-1.d0
       phi = 2.d0*pi*x(2)
      else if (vegasint.ne.1) then
       ctheta = 2.d0*random(idum)-1.d0
       phi = 2.d0*pi*random(idum)
      endif

      stheta = dsqrt(1.d0-ctheta**2)

c The momenta of the final state particles are set by the random variables
      pf1(0) = dsqrt(pfabs**2+mf1**2)
      pf1(1) = pfabs*stheta*dcos(phi)
      pf1(2) = pfabs*stheta*dsin(phi)
      pf1(3) = pfabs*ctheta

      pf2(0) = dsqrt(pfabs**2+mf2**2)
      do x1=1,3
       pf2(x1) = pf1(x1)*gmunu(x1,x1)
      enddo

c Setting up the boost momentum to initial momentum
      do x1=0,3
       piboost(x1) = pinit(x1)*gmunu(x1,x1)
      enddo

      call boost(piboost,pf1)
      call boost(piboost,pf2)

c Phase space weight
      wgt = pfabs/(4.d0*pi*root_s)

      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine phasespace3(pinit,mf1,mf2,mf3,pf1,pf2,pf3,wgt,x,vegasint)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Calculates the three-body decay of an initial state with momentum pinit into three final       c
c  state particles of mass mf1, mf2, and mf3 with momenta pf1, pf2, and pf3.  The phase space     c
c  weight is calculated for the randomly generated momenta of the final state particles.          c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

c input parameters
      double precision pinit(0:3), mf1, mf2, mf3, wgt, x(5)
      integer vegasint

c parameters appearing in this subroutine only
      double precision ctheta1, phi1, ctheta2, phi2, stheta1, stheta2, 
     . mxsq, mxsqmin, mxsqmax, pxyabs, p12abs, px(0:3),
     . root_s, pxboost(0:3), piboost(0:3)

c checking if this is a valid three body decay
      root_s = invmass(pinit)

      if (root_s.lt.(mf1+mf2+mf3)) then
       write(*,*) 'Not a valid three body decay'
       stop
      endif

      mxsqmin = (mf1+mf2)**2
      mxsqmax = (root_s-mf3)**2

c if vegas is used to integrate then use the input array x
      if (vegasint.eq.1) then
       mxsq = (mxsqmax-mxsqmin)*x(1)+mxsqmin
       ctheta1 = 2.d0*x(2)-1.d0
       phi1 = 2.d0*pi*x(3)
       ctheta2 = 2.d0*x(4)-1.d0
       phi2 = 2.d0*pi*x(5)
      else if (vegasint.ne.1) then

c randomly generate values for 5 free parameters in 3 body phase space
       mxsq = (mxsqmax-mxsqmin)*random(idum)+mxsqmin
       ctheta1 = 2.d0*random(idum)-1.d0
       phi1 = 2.d0*pi*random(idum)
       ctheta2 = 2.d0*random(idum)-1.d0
       phi2 = 2.d0*pi*random(idum)
      endif

      stheta1 = dsqrt(1.d0-ctheta1**2)
      stheta2 = dsqrt(1.d0-ctheta2**2)

      pxyabs = dsqrt(lambda(root_s**2,mxsq,mf3**2))/(2.d0*root_s)
      p12abs = dsqrt(lambda(mxsq,mf1**2,mf2**2))/(2.d0*dsqrt(mxsq))

c the momentum of the X particle is set by the free params
      px(0) = dsqrt(pxyabs**2+mxsq)
      px(1) = pxyabs*stheta1*dcos(phi1)
      px(2) = pxyabs*stheta1*dsin(phi1)
      px(3) = pxyabs*ctheta1

c energy momentum conservation of the pinit -> px + pf3 system
      pf3(0) = dsqrt(pxyabs**2+mf3**2)
      do x1=1,3
       pf3(x1) = px(x1)*gmunu(x1,x1)
      enddo

c the momentum of the 1st daughter particle is set by the mass of the
c X system and the last 2 free parameters
      pf1(0) = dsqrt(p12abs**2+mf1**2)
      pf1(1) = p12abs*stheta2*dcos(phi2)
      pf1(2) = p12abs*stheta2*dsin(phi2)
      pf1(3) = p12abs*ctheta2

c energy momentum conservation of the X -> f2 + f3 system
      pf2(0) = dsqrt(p12abs**2+mf2**2)
      do x1=1,3
       pf2(x1) = pf1(x1)*gmunu(x1,x1)
      enddo

c setting up the boost vectors needed
      do x1=0,3
       pxboost(x1) = px(x1)*gmunu(x1,x1)
       piboost(x1) = pinit(x1)*gmunu(x1,x1)
      enddo

c boost p1,p2 into the rest frame of their mother particle
      call boost(pxboost,pf1)
      call boost(pxboost,pf2)

c now boost p1,p2,p3 into the lab frame using original mother momentum
      call boost(piboost,pf1)
      call boost(piboost,pf2)
      call boost(piboost,pf3)

c phase space weight
      wgt = pxyabs*p12abs*(mxsqmax-mxsqmin)/
     .  (32.d0*pi**3*dsqrt(mxsq)*root_s)

      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine phasespace4(pinit,mf1,mf2,mf3,mf4,pf1,pf2,pf3,pf4,wgt,x,vegasint)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Calculates the three-body decay of an initial state with momentum pinit into three final       c
c  state particles of mass mf1, mf2, and mf3 with momenta pf1, pf2, and pf3.  The phase space     c
c  weight is calculated for the randomly generated momenta of the final state particles.          c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

c input parameters
      double precision pinit(0:3), mf1, mf2, mf3, mf4, wgt, x(8)
      integer vegasint

c parameters appearing in this subroutine only
      double precision ctheta1, phi1, ctheta2, phi2, ctheta3, phi3,
     . stheta1, stheta2, stheta3, mxsq, mxsqmax, mxsqmin, 
     . mysq, mysqmax, mysqmin, pxyabs, p12abs, p34abs, px(0:3), py(0:3),
     . root_s, pxboost(0:3), pyboost(0:3), piboost(0:3)

c setting the mass of the mother particle and the max momentum of d1
      root_s = invmass(pinit)

      if (root_s.lt.(mf1+mf2+mf3+mf4)) then
       write(*,*) 'Not a valid four body decay'
       stop
      endif

      mxsqmin = (mf1+mf2)**2
      mysqmin = (mf3+mf4)**2
      mxsqmax = (root_s-dsqrt(mysqmin))**2
      mysqmax = (root_s-dsqrt(mxsqmin))**2

c if vegas is used to integrate then use the input array x
      if (vegasint.eq.1) then
       mxsq = (mxsqmax-mxsqmin)*x(1)+mxsqmin
       mysq = (mysqmax-mysqmin)*x(2)+mysqmin
       ctheta1 = 2.d0*x(3)-1.d0
       phi1 = 2.d0*pi*x(4)
       ctheta2 = 2.d0*x(5)-1.d0
       phi2 = 2.d0*pi*x(6)
       ctheta3 = 2.d0*x(7)-1.d0
       phi3 = 2.d0*pi*x(8)
      else if (vegasint.ne.1) then

c randomly generate values for 5 free parameters in 3 body phase space
       mxsq = (mxsqmax-mxsqmin)*random(idum)+mxsqmin
       mysq = (mysqmax-mysqmin)*random(idum)+mysqmin
       ctheta1 = 2.d0*random(idum)-1.d0
       phi1 = 2.d0*pi*random(idum)
       ctheta2 = 2.d0*random(idum)-1.d0
       phi2 = 2.d0*pi*random(idum)
       ctheta3 = 2.d0*random(idum)-1.d0
       phi3 = 2.d0*pi*random(idum)
      endif

c sets the sine of the thetas
      stheta1 = dsqrt(1.d0-ctheta1**2)
      stheta2 = dsqrt(1.d0-ctheta2**2)
      stheta3 = dsqrt(1.d0-ctheta3**2)

c checks to see if randomly generated mxsq and mysq is valid
      if ((mxsq+mysq).gt.root_s**2) then
       wgt = 0.d0
       return
      endif

c sets the absolute value of the momenta we will be using
      pxyabs = dsqrt(lambda(root_s**2,mxsq,mysq))/(2.d0*root_s)
      p12abs = dsqrt(lambda(mxsq,mf1**2,mf2**2))/(2.d0*dsqrt(mxsq))
      p34abs = dsqrt(lambda(mysq,mf3**2,mf4**2))/(2.d0*dsqrt(mysq))

c setting up the X and Y 4-momenta from the free params
      px(0) = dsqrt(pxyabs**2+mxsq)
      px(1) = pxyabs*stheta1*dcos(phi1)
      px(2) = pxyabs*stheta1*dsin(phi1)
      px(3) = pxyabs*ctheta1

      py(0) = dsqrt(pxyabs**2+mysq)
      do x1=1,3
       py(x1) = px(x1)*gmunu(x1,x1)
      enddo

c the momentum of the daughter particles are set by the other free params
      pf1(0) = dsqrt(p12abs**2+mf1**2)
      pf1(1) = p12abs*stheta2*dcos(phi2)
      pf1(2) = p12abs*stheta2*dsin(phi2)
      pf1(3) = p12abs*ctheta2

      pf3(0) = dsqrt(p34abs**2+mf3**2)
      pf3(1) = p34abs*stheta3*dcos(phi3)
      pf3(2) = p34abs*stheta3*dsin(phi3)
      pf3(3) = p34abs*ctheta3

c momentum 4-vectors of the other daughter particles
      pf2(0) = dsqrt(p12abs**2+mf2**2)
      pf4(0) = dsqrt(p34abs**2+mf4**2)
      do x1=1,3
       pf2(x1) = pf1(x1)*gmunu(x1,x1)
       pf4(x1) = pf3(x1)*gmunu(x1,x1)
      enddo

c setting up the boost vectors needed
      do x1=0,3
       pxboost(x1) = px(x1)*gmunu(x1,x1)
       pyboost(x1) = py(x1)*gmunu(x1,x1)
       piboost(x1) = pinit(x1)*gmunu(x1,x1)
      enddo

c boost p1,p2,p3,p4 into the rest frame of their respective mother particle
      call boost(pxboost,pf1)
      call boost(pxboost,pf2)

      call boost(pyboost,pf3)
      call boost(pyboost,pf4)

c now boost p1,p2,p3,p4 into the lab frame using original mother momentum
      call boost(piboost,pf1)
      call boost(piboost,pf2)
      call boost(piboost,pf3)
      call boost(piboost,pf4)

c phase space weight
      wgt = p12abs*p34abs*pxyabs*(mxsqmax-mxsqmin)*(mysqmax-mysqmin)/
     .  (256.d0*pi**5*dsqrt(mxsq*mysq)*root_s)   

      return
      end