c-------------------------------------------------------------------------------------------------c
      subroutine odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqs)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Runge-Kutta driver with adaptive stepsize control.  Integrate the starting values              c
c  ystart(1:nvar) from x1 to x2 with accuracy eps, storing intermediate results in the common     c
c  block /path/.  h1 should be set as a guessed first stepsize, hmin as the minimum allowed       c
c  stepsize (can be zero). On output nok and nbad are the number of good and bad (but retried     c
c  and fixed) steps taken, and ystart is replaced by values at the end of the integration         c
c  interval.  derivs is the user-supplied subroutine to be used. /path/ contains its own info     c
c  about how often an intermediate value is stored.                                               c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none 
      integer nbad,nok,nvar,kmaxx,maxstp,nmax
      double precision eps,h1,hmin,x1,x2,ystart(nvar),tiny
      external derivs,rkqs
      parameter (maxstp=100000,nmax=50,kmaxx=200,tiny=1.0d-30)
      integer i,kmax,kount,nstp
      double precision dxsav,h,hdid,hnext,x,xsav,dydx(nmax),xp(kmaxx),y(nmax),
     & yp(nmax,kmaxx),yscal(nmax)
      common /path/ kmax,kount,dxsav,xp,yp
      x = x1
      h = dsign(h1,x2-x1)
      nok = 0
      nbad = 0
      kount = 0
      do i=1,nvar
        y(i)=ystart(i)
      enddo
      if (kmax .gt. 0) xsav = x-2.d0*dxsav
      do nstp=1,maxstp
        call derivs(x,y,dydx)
        do i=1,nvar
          yscal(i)=dabs(y(i))+dabs(h*dydx(i))+TINY
        enddo
        if (kmax .gt. 0) then
          if (dabs(x-xsav) .gt. dabs(dxsav)) then
            if(kount .lt. kmax-1) then
              kount = kount+1
              xp(kount) = x
              do i=1,nvar
                yp(i,kount) = y(i)
              enddo
              xsav = x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1) .gt. 0.0d0) h = x2-x
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        if(hdid.eq.h) then
          nok = nok+1
        else
          nbad = nbad+1
        endif
        if((x-x2)*(x2-x1) .ge. 0.0d0) then
          do i=1,nvar
            ystart(i) = y(i)
          enddo
          if(kmax.ne.0) then
            kount = kount+1
            xp(kount) = x
            do i = 1,nvar
              yp(i,kount) = y(i)
            enddo
          endif
          return
        endif
        if(dabs(hnext) .lt. hmin) pause
     &     'stepsize smaller than minimum in odeint'
        h = hnext
        enddo
      pause 'too many steps in odeint'

      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure accuracy and  c
c  adjust stepsize.  Input are the dependent variable y vector(1:n) and its derivative dydx(1:n)  c
c  at the starting value of the independent variable x.  Also input are the stepsize to be        c
c  attempted htry, the required accuracy eps, and the vector yscal(1:n) against which the error   c
c  is scaled.  On output, y and x are replaced by their new values, hdid is the stepsize that was c
c  actually accomplished, and hnext is the estimated next stepsize.  derivs is the user supplied  c
c  subroutine that computes the right-hand side derivatives.                                      c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      integer n,nmax
      double precision eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      external derivs
      parameter (nmax=50)

c     Maximum number of equations
      integer i
      double precision errmax,h,htemp,xnew,yerr(nmax),ytemp(nmax),safety,pgrow,
     &  pshrink,errcon
      parameter (safety=0.9d0,pgrow=-.2d0,pshrink=-.25d0,errcon=1.89d-4)
      h = htry
    1 call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.0d0
      do i=1,n
        errmax = max(errmax,dabs(yerr(i)/yscal(i)))
      enddo
      errmax = errmax/eps
      if(errmax .gt. 1.0d0) then
        htemp = safety*h*(errmax**pshrink)
        h = dsign(max(dabs(htemp), 0.1d0*dabs(h)), h)
        xnew = x+h
        if (xnew.eq.x) pause 'stepsize underflow in rkqs'
        goto 1
      else
        if(errmax.gt.ERRCON) then
          hnext = safety*h*(errmax**pgrow)
        else
          hnext = 5.0d0*h
        endif
        hdid = h
        x = x+h
        do i=1,n
          y(i) = ytemp(i)
        enddo
        return
      endif
      end



c-------------------------------------------------------------------------------------------------c
      subroutine rkck(y,k1,n,x,h,yout,yerr,derivs)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This subroutine calculates a step in solving an ODE using the Runga-Kutta-Cash-Karp method.    c                 
c                                                                                                 c
c  derivs = external subroutine for calculating derivatives                                       c
c  n = number of ODEs being solved                                                                c
c  y = inital point y in the step                                                                 c
c  x = initial point x in the step                                                                c
c  h = stepsize                                                                                   c
c  yout = output of the new y-value                                                               c
c  yerr = estimated error of the step                                                             c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

c input parameters
      integer n
      double precision x,y(n),h,yout(n),yerr(n)
      external derivs

c parameters for this subroutine only
      integer i
      double precision k1(n),k2(n),k3(n),k4(n),k5(n),k6(n),ytemp(n)
      double precision c2,c3,c4,c5,c6
      double precision a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65
      double precision b1,b2,b3,b4,b5,b6
      double precision b1s,b2s,b3s,b4s,b5s,b6s
      parameter (c2=1.d0/5.d0,c3=3.d0/10.d0,c4=3.d0/5.d0,c5=1.d0,c6=7.d0/8.d0,a21=1.d0/5.d0,
     . a31=3.d0/40.d0,a32=9.d0/40.d0,a41=3.d0/10.d0,a42=-9.d0/10.d0,a43=6.d0/5.d0,
     . a51=-11.d0/54.d0,a52=5.d0/2.d0,a53=-70.d0/27.d0,a54=-35.d0/27.d0,
     . a61=1631.d0/55296.d0,a62=175.d0/512.d0,a63=575.d0/13824.d0,a64=44275.d0/110592.d0,
     . a65=253.d0/4096.d0,b1=37.d0/378.d0,b2=0.d0,b3=250.d0/621.d0,b4=125.d0/594.d0,
     . b5=0.d0,b6=512.d0/1771.d0,b1s=2825.d0/27648.d0,b2s=0.d0,b3s=18575.d0/48384.d0,
     . b4s=13525.d0/55296.d0,b5s=277.d0/14336.d0,b6s=1.d0/4.d0)

c many setup calculations are required
c      call derivs(x,y,k1)
      do i=1,n
        ytemp(i) = y(i) + h*a21*k1(i)
      enddo
      call derivs(x+h*c2,ytemp,k2)
      do i=1,n
        ytemp(i) = y(i) + h*(a31*k1(i) + a32*k2(i))
      enddo
      call derivs(x+h*c3,ytemp,k3)
      do i=1,n
        ytemp(i) = y(i) + h*(a41*k1(i) + a42*k2(i) + a43*k3(i))
      enddo
      call derivs(x+h*c4,ytemp,k4)
      do i=1,n
        ytemp(i) = y(i) + h*(a51*k1(i) + a52*k2(i) + a53*k3(i) + a54*k4(i))
      enddo
      call derivs(x+h*c5,ytemp,k5)
      do i=1,n
        ytemp(i) = y(i) + h*(a61*k1(i) + a62*k2(i) + a63*k3(i) + a64*k4(i) + a65*k5(i))
      enddo
      call derivs(x+h*c6,ytemp,k6)

c time to calculate the result and error
      do i=1,n
        yout(i) = y(i) + h*(b1s*k1(i) + b2s*k2(i) + b3s*k3(i) + b4s*k4(i) + b5s*k5(i) + b6s*k6(i))
        yerr(i) = h*((b1-b1s)*k1(i) + (b2-b2s)*k2(i) + (b3-b3s)*k3(i) + (b4-b4s)*k4(i) + (b5-b5s)*k5(i) + (b6-b6s)*k6(i))
      enddo

      return
      end