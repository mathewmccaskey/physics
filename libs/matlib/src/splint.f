c-------------------------------------------------------------------------------------------------c
      subroutine splint(xa,ya,y2a,n,x,y)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the xa's in  c
c  order), and given the array y2a(1:n), which is the output from spline and given a value of x,  c
c  this routine returns a cubic spline interpolated value y.                                      c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

      integer n
      double precision x,y,xa(n),y2a(n),ya(n)
      integer k,khi,klo
      double precision a,b,h

c We find the right place in the table by means of bisection
      klo = 1
      khi = n
 1    if (khi-klo.gt.1) then
        k = (khi+klo)/2
        if (xa(k).gt.x) then
         khi = k
        else
         klo = k
        endif
      goto 1
      endif
      h = xa(khi)-xa(klo)
      if (h.eq.0.d0) pause 'bad xa input in splint'
      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      y = a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*
     .  (h**2)/6.d0

      return
      end