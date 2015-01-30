c-------------------------------------------------------------------------------------------------c
      subroutine simpson(func,a,b,s,n)
c-------------------------------------------------------------------------------------------------c
c						                                                                                      c
c  This is my own subroutine for the simpson's rule for numerical integration of a one            c
c  dimensional integral.  apparently the simpson's rule numerical recipe has difficulty           c
c  integrating the function I give for some reason so boo on that.				                        c
c 							                                                                                  c
c-------------------------------------------------------------------------------------------------c
      implicit none

c input parameteres
      double precision func,a,b,s
      integer n
      external func

c parameters used in this subroutine only   
      double precision dx, x
      integer i

c checking to see that n is even
      if (mod(n,2).eq.1) n = n+1
      dx = (b-a)/dble(n)
      s = (func(a)+func(b))
      do i=2,n-2,2
        x = a + dble(i)*dx
        s = s + 2.d0*func(x)
      enddo
      do i=1,n-1,2
        x = a + dble(i)*dx
        s = s + 4.d0*func(x)
      enddo
      s = (dx/3.d0)*s

      return
      end