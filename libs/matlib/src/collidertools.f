c-------------------------------------------------------------------------------------------------c
      subroutine boost(p,q)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Boosts the 4-momentum vector q into the rest frame of the 4-momentum vector p.                 c
c                                                                                                 c
c  Note: inputs p and q cannot come from the same variable. i.e. call boost(p1,p1) does not       c
c  work. Instead use a dummy momentum: call boost(p_dum,p1) where p_dum = p1.                     c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

c Parameters appearing in the this subroutine only
      double precision mp, gam, qnew(0:3)

c Check if the rest frame can be boosted into
      mp = invmass(p)

      if (mp.le.1d-20) then
       write(*,*) 'Cannot boost into rest frame of a massless particle'
       stop
      endif

      gam = p(0)/mp

      qnew(0) = p4dot(p,q)/mp
      do x1=1,3
       qnew(x1) = q(x1) + (p3dot(p,q)/(mp*(p(0)+mp))-q(0)/mp)*p(x1)
      enddo

      do x1=0,3
       q(x1) = qnew(x1)
      enddo

      return
      end



c-------------------------------------------------------------------------------------------------c
      function deltaphi(p,q)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Returns the delta phi between two momentum vectors p and q.                                    c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

      if ((pt(p).lt.1d-20).or.(pt(q).lt.1d-20)) then
       pause 'Cannot calculate delta phi with zero pt'
       stop
      endif

      deltaphi = dacos((p(1)*q(1)+p(2)*q(2))/(pt(p)*pt(q)))

      return
      end



c-------------------------------------------------------------------------------------------------c
      function deltaR(p,q)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Returns the delta R between two 4-vectors p and q                                              c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

      deltaR = dsqrt((rapidity(p)-rapidity(q))**2+deltaphi(p,q)**2)

      return
      end



c-------------------------------------------------------------------------------------------------c
      function Et(p)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Returns the transverse energy of a 4-vector p                                                  c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

      Et = dsqrt(invmass(p)**2 + pt(p)**2)

      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine getangles(p,theta,phi)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Subroutine that gives the angles theta and phi in spherical coordinates given an incoming      c
c  4-momentum p                                                                                   c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

      if (p3dot(p,p).gt.0.d0) then
       theta = dacos(p(3)/dsqrt(p3dot(p,p)))
      else 
       theta = 0.d0
      endif

      if (pt(p).gt.0.d0) then
       phi = dacos(p(1)/pt(p))
      else
       phi = 0.d0
      endif

      if (p(2).lt.0.d0) phi = 2.d0*pi-phi

      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine getpolarization(p,spin,pol)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Given the four momentum p and the polarization desired this subroutine gives the polarization  c
c  vector pol                                                                                     c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

c input parameters
      integer spin
      double complex pol(0:3)

c initilize the polarization vector
      do x1=0,3
       pol(x1) = dcmplx(0.d0,0.d0)
      enddo

c zero momentum check
      if (p3dot(p,p).eq.0.d0) then
       if (spin.eq.1) pol(1) = dcmplx(1.d0,0.d0)
       if (spin.eq.2) pol(2) = dcmplx(1.d0,0.d0)
       if (spin.eq.3) pol(3) = dcmplx(1.d0,0.d0)
       return
      endif

c helicity convention, spin=1 for minus, spin=2 for plus, 
c spin=3 for longitudinal (massive only)
      if (spin.eq.1) then
       if (pt(p).eq.0.d0) then
        pol(1) = dcmplx(1.d0/dsqrt(2.d0),0.d0)
        pol(2) = dcmplx(0.d0,-1.d0/dsqrt(2.d0))
       else
        pol(1) = dcmplx(p(1)*p(3),p(2)*dsqrt(p3dot(p,p)))/
     .    (dsqrt(2.d0*p3dot(p,p))*pt(p))
        pol(2) = dcmplx(p(2)*p(3),-1.d0*p(1)*dsqrt(p3dot(p,p)))/
     .    (dsqrt(2.d0*p3dot(p,p))*pt(p))
        pol(3) = dcmplx(-1.d0*pt(p),0.d0)/dsqrt(2.d0*p3dot(p,p))
       endif

      else if (spin.eq.2) then
       if (pt(p).eq.0.d0) then
        pol(1) = dcmplx(-1.d0/dsqrt(2.d0),0.d0)
        pol(2) = dcmplx(0.d0,-1.d0/dsqrt(2.d0))
       else
        pol(1) = dcmplx(-1.d0*p(1)*p(3),p(2)*dsqrt(p3dot(p,p)))/
     .    (dsqrt(2.d0*p3dot(p,p))*pt(p))
        pol(2) = -1.d0*dcmplx(p(2)*p(3),p(1)*dsqrt(p3dot(p,p)))/
     .    (dsqrt(2.d0*p3dot(p,p))*pt(p))
        pol(3) = dcmplx(pt(p),0.d0)/dsqrt(2.d0*p3dot(p,p))
       endif

      else if (spin.eq.3) then
       if (invmass(p).eq.0.d0) then
        write(*,*) 'Cannot have long. pol. with massless particle'
        stop
       endif
       pol(0) = dcmplx(dsqrt(p3dot(p,p))/invmass(p),0.d0)
       pol(1) = dcmplx(p(1)*p(0)/(invmass(p)*dsqrt(p3dot(p,p))),0.d0)
       pol(2) = dcmplx(p(2)*p(0)/(invmass(p)*dsqrt(p3dot(p,p))),0.d0)
       pol(3) = dcmplx(p(3)*p(0)/(invmass(p)*dsqrt(p3dot(p,p))),0.d0)
      endif

      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine getu(p,spin,uout)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Given the four momentum p and the helicity spin this subroutine gives the particle wave        c
c  function u                                                                                     c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

c input parameters
      integer spin

c parameters used in this subroutine only
      double precision absp

c getting the magnitude of momentum and polar angles
      absp = dsqrt(p3dot(p,p))
      call getangles(p,theta,phi)

      do x1=1,4
       uout(x1) = dcmplx(0.d0,0.d0)
      enddo

c helicity convention, spin=1 for minus, spin=2 for plus
      if (spin.eq.1) then
       uout(1) = -1.d0*dsqrt(p(0)+absp)*
     .  dcmplx(dcos(phi)*dsin(theta/2.d0),
     .    -1.d0*dsin(phi)*dsin(theta/2.d0))
       uout(2) = dsqrt(p(0)+absp)*dcmplx(dcos(theta/2.d0),0.d0)

c if absp > p(0) then we will assume the particle is massless
       if (absp.lt.p(0)) then
        uout(3) = -1.d0*dsqrt(p(0)-absp)*
     .   dcmplx(dcos(phi)*dsin(theta/2.d0),
     .    -1.d0*dsin(phi)*dsin(theta/2.d0))
        uout(4) = dsqrt(p(0)-absp)*dcmplx(dcos(theta/2.d0),0.d0)
       endif
      else if (spin.eq.2) then
       if (absp.lt.p(0)) then
        uout(1) = dsqrt(p(0)-absp)*dcmplx(dcos(theta/2.d0),0.d0)
        uout(2) = dsqrt(p(0)-absp)*
     .   dcmplx(dcos(phi)*dsin(theta/2.d0),dsin(phi)*dsin(theta/2.d0))
       endif
       uout(3) = dsqrt(p(0)+absp)*dcmplx(dcos(theta/2.d0),0.d0)
       uout(4) = dsqrt(p(0)+absp)*
     .  dcmplx(dcos(phi)*dsin(theta/2.d0),dsin(phi)*dsin(theta/2.d0))
      else
       write(*,*) 'Not a valid fermion spin'
       stop
      endif

      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine getv(p,spin,vout)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Given the four momentum p and the helicity spin this subroutine gives the anti-particle wave   c
c  function v                                                                                     c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

c input parameters
      integer spin

c parameters used in this subroutine only
      double precision absp

c getting the magnitude of momentum and polar angles
      absp = dsqrt(p3dot(p,p))
      call getangles(p,theta,phi)

      do x1=1,4
       vout(x1)=dcmplx(0.d0,0.d0)
      enddo

c helicity convention, spin=1 for minus, spin=2 for plus
      if (spin.eq.1) then

c if absp > p(0) then we will assume that the particle is massless
       if (absp.lt.p(0)) then
        vout(1) = -1.d0*dsqrt(p(0)-absp)*dcmplx(dcos(theta/2.d0),0.d0)
        vout(2) = -1.d0*dsqrt(p(0)-absp)*
     .   dcmplx(dcos(phi)*dsin(theta/2.d0),dsin(phi)*dsin(theta/2.d0))
       endif
       vout(3) = dsqrt(p(0)+absp)*dcmplx(dcos(theta/2.d0),0.d0)
       vout(4) = dsqrt(p(0)+absp)*
     .  dcmplx(dcos(phi)*dsin(theta/2.d0),dsin(phi)*dsin(theta/2.d0))
      else if (spin.eq.2) then
       vout(1) = -1.d0*dsqrt(p(0)+absp)*
     .  dcmplx(dcos(phi)*dsin(theta/2.d0),
     .    -1.d0*dsin(phi)*dsin(theta/2.d0))
       vout(2) = dsqrt(p(0)+absp)*dcmplx(dcos(theta/2.d0),0.d0)
       if (absp.lt.p(0)) then
        vout(3) = dsqrt(p(0)-absp)*
     .   dcmplx(dcos(phi)*dsin(theta/2.d0),
     .    -1.d0*dsin(phi)*dsin(theta/2.d0))
        vout(4) = -1.d0*dsqrt(p(0)-absp)*dcmplx(dcos(theta/2.d0),0.d0)
       endif
      else
       write(*,*) 'Not a valid fermion spin'
       stop
      endif

      return
      end



c-------------------------------------------------------------------------------------------------c
      function gmunu(mu,nu)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Returns the minkowski metric                                                                   c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

c input parameters 
      integer mu, nu

      if ((mu.gt.3).or.(mu.lt.0).or.(nu.gt.3).or.(nu.lt.0)) then
       pause 'Not a valid index for gmunu'
       stop
      endif

      if (mu.eq.nu) then
       if (mu.eq.0) then
        gmunu = 1.d0
       else
        gmunu = -1.d0
       endif
      else
       gmunu = 0.d0
      endif

      return
      end

       

c-------------------------------------------------------------------------------------------------c
      function invmass(p)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Calculates the invariant mass from a momentum 4-vector                                         c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

      invmass = p4dot(p,p)

      if (invmass.le.-1.d0) then
       pause 'Negative invariant mass squared:'
c       pause invmass
       stop
      else if (invmass.le.1.d-10) then
       invmass = 0.d0
      else
       invmass = dsqrt(invmass)
      endif

      return
      end



c-------------------------------------------------------------------------------------------------c
      function lambda(a,b,c)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Calculate the lambda function which shows up in momentum and phase space calculations.         c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

c input parameters
      double precision a, b, c

      lambda = a**2 + b**2 + c**2 - 2.d0*(a*b + a*c + b*c)

      return
      end




c-------------------------------------------------------------------------------------------------c
      function p3dot(p,q)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Returns the 3-vector dot product                                                               c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

      p3dot = p(1)*q(1)+p(2)*q(2)+p(3)*q(3)

      return
      end



c-------------------------------------------------------------------------------------------------c
      function p4dot(p,q)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Returns the 4-vector dot product                                                               c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

      p4dot = p(0)*q(0)-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)
       
      return
      end



c-------------------------------------------------------------------------------------------------c
       subroutine plus(p,q,pf)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Returns the addition of two 4-vectors (p1,p2)                                                  c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
       implicit none
       include 'matlib.inc'

c input parameters
       double precision pf(0:3)

       do x1=0,3
        pf(x1) = p(x1) + q(x1)
       enddo

       return
       end



c-------------------------------------------------------------------------------------------------c
      subroutine projectleft(uin,uout)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Subroutine to give the left handed projection of wave function uin.                            c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

c initializing the output wave functioin
      do x1=1,4
       uout(x1) = dcmplx(0.d0,0.d0)
      enddo

c set the first two values of output to the wave function uin
      uout(1) = uin(1)
      uout(2) = uin(2)

      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine projectright(uin,uout)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Subroutine to give the right handed projection of wave function uin.                           c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

c initializing the output wave functioin
      do x1=1,4
       uout(x1) = dcmplx(0.d0,0.d0)
      enddo

c set the second two values of output to the wave function uin
      uout(3) = uin(3)
      uout(4) = uin(4)

      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine pslash(vec,matrix)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Subroutine to create a slashed matrix with an input four vector vec.  The final output is in   c
c  matrix. Requires that the subroutine initlib has been called.                                  c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

c input parameters
      double complex vec(0:3),matrix(0:3,0:3)

      call setgamma(gamma0,gamma1,gamma2,gamma3)

      do x1=0,3
       do x2=0,3
        matrix(x1,x2) = gamma0(x1,x2)*vec(0)-gamma1(x1,x2)*vec(1)-
     .   gamma2(x1,x2)*vec(2)-gamma3(x1,x2)*vec(3)
       enddo
      enddo

      return
      end



c-------------------------------------------------------------------------------------------------c
      function pt(p)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Returns the transverse momentum of 4-vector p                                                  c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

      pt = dsqrt(p(1)*p(1)+p(2)*p(2))

      return
      end



c-------------------------------------------------------------------------------------------------c
      function rapidity(p)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Returns the rapidity of the momentum 4-vector p                                                c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

      rapidity = 0.5d0*dlog(max(1d-20,(p(0)+p(3)))/
     .                      max(1d-20,(p(0)-p(3))))

      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine smearp(p,a,b)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Energy smearing of the form delta E/E = a/sqrt(E) + b                                          c
c  p is replaced with the smeared p.                                                              c
c                                                                                                 c
c  for jets: a = 0.5, b = 0.03                                                                    c
c  for leptons: a = 0.1, b = 0.007                                                                c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'
      
c parameters used in this subroutine only
      double precision a,b,del
      
      del = randgauss(idum,0.d0,dsqrt(a**2/p(0)+b**2))
      
      do x1=0,3
        p(x1) = p(x1)*(1.d0+del)
      enddo
      
      return
      end
     
      
      
c-------------------------------------------------------------------------------------------------c
      function theta3(p,q)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Returns the angle between the moentum 3-vectors of the 4-vectors p1,p2                         c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

      if ((p3dot(p,p).lt.1.d-30).or.(p3dot(q,q).lt.1.d-30)) then
       pause 'Cannot calculate angle with zero momentum'
       stop
      endif

      theta3 = dacos(p3dot(p,q)/dsqrt(p3dot(p,p)*p3dot(q,q)))

      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine ubar(uin,ubarout)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Given an input wave function u this subroutine yields the u bar wave function                  c
c  (u^dagger*gamma^0).                                                                            c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

c parameters used in this subroutine only
      double complex dum(1,4)

      call setgamma(gamma0,gamma1,gamma2,gamma3)

      do x1=1,4
       dum(1,x1) = dcmplx(0.d0,0.d0)
       ubarout(1,x1) = dcmplx(0.d0,0.d0)
      enddo

c ubar = u^dagger*gamma^0
      call dagger(4,1,uin,dum)
      call matrixmultiplycmplx(1,4,4,4,dum,gamma0,ubarout)

      return
      end