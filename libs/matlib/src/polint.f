c-------------------------------------------------------------------------------------------------c
      subroutine polint(xa,ya,n,x,y,dy)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Given arrays xa and ya, each of length n, and given a value x, this routine returns a value    c
c  y, and an error estimate dy.  If P(x) is the polynomial of degree N=1 such that                c
c  P(xa_i) = ya_i, i=1...n then the returned value y = P(x).                                      c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      integer n,nmax
      double precision dy,x,y,xa(n),ya(n)
      parameter (nmax=10)
      integer i,m,ns
      double precision den,dif,dift,ho,hp,w,c(nmax),d(nmax)
      ns = 1
      dif = dabs(x-xa(1))

c Here we find the index ns of the closest table entry
      do i=1,n
        dift = dabs(x-xa(i))
        if (dift.lt.dif) then
          ns = i
          dif = dift
        endif

c and initialize the tableau of c's and d's
        c(i) = ya(i)
        d(i) = ya(i)
      enddo

c This is the initial approximation to y
      y = ya(ns)
      ns = ns-1

c For each column of the tableau,
c we loop over the current c's and d's and update them
      do m=1,n-1
        do i=1,n-m
          ho = xa(i)-x
          hp = xa(i+m)-x
          w = c(i+1)-d(i)
          den = ho-hp

c This error can occur only if two input xa's are (to within roundoff) identical.
          if (den .eq. 0.d0) pause 'failure in polint'
          den = w/den

c Here the c's and d's are updated
          d(i) = hp*den
          c(i) = ho*den
        enddo
        if(2*ns .lt. n-m) then
          dy = c(ns+1)
        else
          dy = d(ns)
          ns = ns-1
        endif
        y = y+dy
      enddo

      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine polint2(x1a,x2a,ya,m,n,x1,x2,y,dy)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Given arrays x1a(1:m) and x2a(1:n) of independant variables, and an m by m array of function   c
c  values ya(1:m c,1:n), tabulated at the grid points defined by x1a and x2a; and given values    c
c  x1 and x2 of the independent variablesl this routine returns an interpolated function value    c
c  y, and an accuract indication dy (based only on the interpolation in the x1 direction,         c
c  however).                                                                                      c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

      integer m,n,nmax,mmax
      double precision dy,x1,x2,y,x1a(m),x2a(n),ya(m,n)
      parameter(nmax=20,mmax=20)

c uses polint
      integer j,k
      double precision ymtmp(mmax), yntmp(nmax)

c loop over rows
      do j=1,m

c Copy row into tempoarary storage
        do k=1,n
          yntmp(k) = ya(j,k)
        enddo

c Interpolate answer into temporary storage
        call polint(x2a,yntmp,n,x2,ymtmp(j),dy)
      enddo

c Do the final interpolation
      call polint(x1a,ymtmp,m,x1,y,dy)

      return
      end