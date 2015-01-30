c-------------------------------------------------------------------------------------------------c
      function brent(ax,bx,cx,f,tol,xmin)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Given a function f and a triplet of numbers ax, bx, and cx, with ax < bx < cx and              c
c  f(bx) < f(ax), f(cx) this function calculates the minimum of f to a precision of tol.  brent   c
c  returns the value of the function while xmin returns the value of x where the minimum occurs.  c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

      integer itmax
      double precision brent, ax, bx, cx, tol, xmin, f, cgold, zeps
      external f
      parameter (itmax=100,cgold=0.3819660d0,zeps=1.d-10)
      integer iter
      double precision a, b, d, e, etemp, fu, fv, fw, fx, p, q, r, tol1,
     . tol2, u, v, w, x, xm

c a and b must be in ascending order though ax and cx need not be
      a = min(ax,cx)
      b = max(ax,cx)

c Initializations
      v = bx
      w = v
      x = v

c This will be the distance moved on the step before last
      e = 0.d0
      fx = f(x)
      fv = fx
      fw = fx

c Main program loop
      do iter=1,itmax
       xm = 0.5d0*(a+b)
       tol1 = tol*dabs(x)+zeps
       tol2 = 2.d0*tol1

c Test to see if function is done
       if (dabs(x-xm).le.(tol2-0.5d0*(b-a))) goto 30
       if (dabs(e).gt.tol1) then

c Construct a trial parabolic fit
        r = (x-w)*(fx-fv)
        q = (x-v)*(fx-fw)
        p = (x-v)*q-(x-w)*r
        q = 2.d0*(q-r)
        if (q.gt.0.d0) p = -p
        q = dabs(q)
        etemp = e
        e = d

c Determine the acceptability of the parabolic fit.
        if ((dabs(p).ge.dabs(0.5d0*q*etemp)).or.(p.le.q*(a-x)).or.
     .      (p.ge.q*(b-x))) goto 10

c Take a parabolic step
        d = p/q
        u = x+d
        if ((u-a.lt.tol2).or.(b-u.lt.tol2)) d = dsign(tol1,xm-x)

c Skip over golden section step
        goto 20
       endif

c We arrive here for a golden section step, which we take the into
c the larger of the two segments
 10    if (x.ge.xm) then
        e = a-x
       else
        e = b-x
       endif

c Take the golden section step
       d = cgold*e

c Arrive here with d computed eiter from parabolic fit, or else from 
c golden section
 20    if (dabs(d).ge.tol1) then
        u = x+d
       else
        u = x+dsign(tol1,d)
       endif

c This is the one function evaluation per iteration
       fu = f(u)

c Housekeeping
       if (fu.le.fx) then
        if (u.ge.x) then
         a = x
        else
         b = x
        endif
        v = w
        fv = fw
        w = x
        fw = fx
        x = u
        fx = fu
       else
        if (u.lt.x) then
         a = u
        else
         b = u
        endif
        if ((fu.le.fw).or.(w.eq.x)) then
         v = w
         fv = fw
         w = u
         fw = fu
        else if ((fu.le.fv).or.(v.eq.x).or.(v.eq.w)) then
         v = u
         fv = fu
        endif
       endif
c Done with housekeeping. Back for another iteration

      enddo
      pause 'brent exceed maximum iterations'

c Arrive here ready to exit with best values
 30   xmin = x
      brent = fx
      
      return
      end