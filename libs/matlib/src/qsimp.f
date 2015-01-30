c-------------------------------------------------------------------------------------------------c
      subroutine qsimp(func,a,b,s)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Returns as s the integral of the function func from a to b.  The parameters eps can be set to  c
c  the desired fractional accuracy and jmax so that 2 to the power jmax-1 is the maximjum         c
c  allowed number of steps.  Integration is performed by Simpson's rule.                          c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

      integer jmax
      double precision a,b,func,s,eps
      external func
      parameter (eps=1.0d-6, jmax=20)

c uses trapzd
      integer j
      double precision os,ost,st

      ost = -1.0d30
      os = -1.0d30
      do j=1,jmax
        call trapzd(func,a,b,st,j)
        s = (4.d0*st-ost)
        if (dabs(s-os).lt.eps*dabs(os)) return
        if (s.eq.0.d0.and.os.eq.0.d0.and.j.gt.6) return
        os = s
        ost = st
      enddo

      pause 'too many steps in qsimp'
      end
