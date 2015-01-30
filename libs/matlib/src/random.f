c-------------------------------------------------------------------------------------------------c
      function random(idum)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Randomly generates a double precision number between 0 and 1.  First set idum to a negative    c
c  number to initialize the seed and afterward just call random(idum)                             c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      double precision random,am,eps,rnmx
      parameter (im1=2147483563,im2=2147483399,am=1.d0/im1,imm1=im1-1,
     . ia1=40014, ia2=40692, iq1=53668, iq2=52774, ir1=12211, ir2=3791,
     . ntab=32, ndiv=1+imm1/ntab, eps=1.2d-7, rnmx=1.d0-eps)
      integer idum2, j, k, iv(ntab), iy
      save iv, iy, idum2
      data idum2/123456789/, iv/ntab*0/, iy/0/

c initialize
      if (idum.le.0) then

c prevents idum = 0
       idum = max(-idum,1)
       idum2 = idum

c Load the shuffle table (after 8 warm-ups)
       do j=ntab+8,1,-1
        k = idum/iq1
        idum = ia1*(idum-k*iq1)-k*ir1
        if (idum.lt.0) idum = idum + im1
        if (j.le.ntab) iv(j) = idum
       enddo
       iy=iv(1)
      endif

c Start here when not initializing
      k = idum/iq1

c Compute idum=mod(ia1*idum,im1) without overflows by Schrage's method
      idum = ia1*(idum-k*iq1)-k*ir1
      if (idum.lt.0) idum = idum + im1
      k = idum/iq2

c Compute idum2=mod(ia2*idum2,im2) likewise
      idum2 = ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2 = idum2 + im2

c Will be in the range 1:ntab
      j = 1+iy/ndiv

c Here idum is shuffled, idum and idum2 are combined to generate output
      iy = iv(j)-idum2
      iv(j) = idum
      if(iy.lt.1) iy = iy + imm1

c The endpoints are not expected in random number generators
      random = min(am*iy,rnmx)
     
      return
      end



c-------------------------------------------------------------------------------------------------c
      function randgauss(idum,mean,sd)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This function returns a random number consistent with a gaussian distribution that has a user  c
c  input mean and standard deviation.                                                             c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

c input parameters
      integer idum
      double precision mean, sd

c parameters for this subroutine only
      integer iset
      double precision random, randgauss, x1, x2
      double precision gset, pi      
      save iset
      data iset/0/ 

      pi = 4.d0*datan(1.d0)
      x1 = random(idum)
      x2 = random(idum)

      if (iset .eq. 0) then
        randgauss = mean + sd*dsqrt(-2.d0*dlog(x1))*dcos(2.d0*pi*x2)
        iset = 1
      else
        randgauss = mean + sd*dsqrt(-2.d0*dlog(x1))*dsin(2.d0*pi*x2)
        iset = 0
      endif

      return
      end      