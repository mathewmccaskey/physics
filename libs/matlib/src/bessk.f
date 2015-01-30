c-------------------------------------------------------------------------------------------------c
      function bessk(n,x)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Returns the modified Bessel function K_{n}(x) for positive x and n >= 2.                       c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      
      integer n
      double precision bessk, x

c uses bessk0, bessk1
      integer j
      double precision bk, bkm, bkp, tox, bessk0, bessk1

      if (n.lt.2) pause 'bad argument n in bessk'
      tox = 2.d0/x
      bkm = bessk0(x)
      bk = bessk1(x)
      do j=1,n-1
        bkp = bkm + j*tox*bk
        bkm = bk
        bk = bkp
      enddo
      bessk = bk

      return
      end