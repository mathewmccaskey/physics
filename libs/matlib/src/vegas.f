c-------------------------------------------------------------------------------------------------c
      subroutine vegas(region,ndim,fxn,init,ncall,itmx,nprn,tgral,sd,chi2a)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Monte Carlo integration of external function fxn                                               c
c  region - 2*ndim vector of the boundaries of integration                                        c 
c  ndim - number of dimensions of the integration                                                 c
c  fxn - external function to integrate over                                                      c
c  init - new grid (0), prev grid w/o (1) or w/ (2) results                                       c
c  ncall - number of times to call the function                                                   c
c  itmx - number of iterations (5)                                                                c
c  nprn - flag for controlling output: nprn < 0 for nothing                                       c
c  tgral - best estimate of the integral                                                          c
c  sd - standard deviation                                                                        c
c  chi2a - chi squared per degree of freedom                                                      c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

      integer init, itmx, ncall, ndim, nprn, ndmx, mxdim
      double precision tgral, chi2a, sd, region(2*ndim), fxn, alph, tiny
      parameter (alph = 1.5d0, ndmx = 50, mxdim = 10, tiny = 1.d-30)
      external fxn

c uses fxn, random, rebin
      integer i,idum,it,j,k,mds,nd,ndo,ng,npg,ia(mxdim),kg(mxdim)
      double precision calls,dv2g,dxg,f,f2,f2b,fb,rc,ti,tsi,
     . wgt,xjac,xn,xnd,xo,d(ndmx,mxdim),di(ndmx,mxdim),dt(mxdim),
     . dx(mxdim),r(ndmx),x(mxdim),xi(ndmx,mxdim),xin(ndmx),random
      double precision schi,si,swgt

c Means for random number initialization
      common/randseed/ idum
      save

c Normal entry. Enter here on cold start
      if (init.le.0) then

c Change mds = 0 to disable stratified sampling (i.e. use importance 
c sampling only)
        mds = 1
        ndo = 1
        do j=1,ndim
          xi(1,j) = 1.d0
        enddo
      endif

c Enter here to inherit the grid from a previous call, but not its answers
      if (init.le.1) then
        si = 0.d0
        swgt = 0.d0
        schi = 0.d0
      endif

c Enter here to inherit the previous grid and its answers
      if (init.le.2) then
        nd = ndmx
        ng = 1

c Set up for stratification
        if (mds.ne.0) then
          ng = (dble(ncall)/2.d0+0.25d0)**(1.d0/dble(ndim))
          mds = 1
          if ((2*ng-ndmx).ge.0) then
            mds = -1
            npg = ng/ndmx + 1
            nd = ng/npg
            ng = npg*nd
          endif
        endif
        k = ng**ndim
        npg = max(ncall/k,2)
        calls = dble(npg)*dble(k)
        dxg = 1.d0/dble(ng)
        dv2g = (calls*dxg**ndim)**2/dble(npg)/dble(npg)/(dble(npg)-1.d0)
        xnd = dble(nd)
        dxg = dxg*xnd
        xjac = 1.d0/calls
        do j=1,ndim
          dx(j) = region(j+ndim)-region(j)
          xjac = xjac*dx(j)
        enddo

c Do binning if necessary
        if (nd.ne.ndo) then
          do i=1,max(nd,ndo)
            r(i) = 1.d0
          enddo
          do j=1,ndim
            call rebin(ndo/xnd,nd,r,xin,xi(1,j))
          enddo
          ndo = nd
        endif
        if (nprn.ge.0) write(*,200) ndim,calls,it,itmx,nprn,ALPH,mds,nd,
     *(j,region(j),j,region(j+ndim),j=1,ndim)
      endif

c Main iteration loop
      do it=1,itmx
        ti = 0.d0
        tsi = 0.d0
        do j=1,ndim
          kg(j) = 1
          do i=1,nd
            d(i,j) = 0.d0
            di(i,j) = 0.d0
          enddo
        enddo
10      continue
        fb = 0.d0
        f2b = 0.d0
        do k=1,npg
          wgt = xjac
          do j=1,ndim
            xn = (kg(j)-random(idum))*dxg+1.d0
            ia(j) = max(min(int(xn),NDMX),1)
            if (ia(j).gt.1) then
              xo = xi(ia(j),j)-xi(ia(j)-1,j)
              rc = xi(ia(j)-1,j)+(xn-ia(j))*xo
            else
              xo = xi(ia(j),j)
              rc = (xn-ia(j))*xo
            endif
            x(j) = region(j)+rc*dx(j)
            wgt = wgt*xo*xnd
          enddo
          f = wgt*fxn(x,wgt)
          f2 = f*f
          fb = fb+f
          f2b = f2b+f2
          do j=1,ndim
            di(ia(j),j) = di(ia(j),j)+f
            if (mds.ge.0) d(ia(j),j)=d(ia(j),j)+f2
          enddo
        enddo
        f2b = dsqrt(f2b*npg)
        f2b = (f2b-fb)*(f2b+fb)
        if (f2b.le.0.d0) f2b=tiny
        ti = ti+fb
        tsi = tsi+f2b
c Use stratified sampling
        if (mds.lt.0) then
          do j=1,ndim
            d(ia(j),j) = d(ia(j),j)+f2b
          enddo
        endif
        do k=ndim,1,-1
          kg(k) = mod(kg(k),ng)+1
          if (kg(k).ne.1) goto 10
        enddo

c Compute final results with this iteration
        tsi = tsi*dv2g
        wgt = 1.d0/tsi
        si = si+dble(wgt)*dble(ti)
        schi = schi+dble(wgt)*dble(ti)**2
        swgt = swgt+dble(wgt)
        tgral = si/swgt
        chi2a = max((schi-si*tgral)/(it-.99d0),0.d0)
        sd=dsqrt(1.d0/swgt)
        tsi=dsqrt(tsi)
        if (nprn.ge.0) then
          write(*,201) it,ti,tsi,tgral,sd,chi2a
          if (nprn.ne.0) then
            do j=1,ndim
              write(*,202) j,(xi(i,j),di(i,j),i=1+nprn/2,nd,nprn)
            enddo
          endif
        endif

c Refine the grid
        do j=1,ndim
          xo = d(1,j)
          xn = d(2,j)
          d(1,j) = (xo+xn)/2.d0
          dt(j) = d(1,j)
          do i=2,nd-1
            rc = xo+xn
            xo = xn
            xn = d(i+1,j)
            d(i,j) = (rc+xn)/3.d0
            dt(j) = dt(j)+d(i,j)
          enddo
          d(nd,j) = (xo+xn)/2.d0
          dt(j) = dt(j)+d(nd,j)
        enddo
        do j=1,ndim
          rc = 0.d0
          do i=1,nd
            if (d(i,j).lt.tiny) d(i,j)=tiny
            r(i) = ((1.d0-d(i,j)/dt(j))/(dlog(dt(j))-dlog(d(i,j))))**ALPH
            rc = rc+r(i)
          enddo
          call rebin(rc/xnd,nd,r,xin,xi(1,j))
        enddo
      enddo

200   format(/' input parameters for vegas:  ndim=',i3,'  ncall=',
     *f8.0/28x,'  it=',i5,'  itmx=',i5/28x,'  nprn=',i3,'  alph=',
     *f5.2/28x,'  mds=',i3,'   nd=',i4/(30x,'xl(',i2,')= ',g11.4,' xu(',
     *i2,')= ',g11.4))
201   format(/' iteration no.',I3,': ','integral =',g14.7,'+/- ',g9.2/
     *' all iterations:   integral =',g14.7,'+/- ',g9.2,
     *' chi**2/it''n =',g9.2)
202   format(/' data for axis ',I2/'    x       delta i       ',
     *'   x       delta i       ','    x       delta i       ',/(1x,
     *f7.5,1x,g11.4,5x,f7.5,1x,g11.4,5x,f7.5,1x,g11.4))

      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine rebin(rc,nd,r,xin,xi)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Utility routine used by vegas to rebin a vector of densities xi into new bins defined by a     c
c  vector r                                                                                       c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

      integer nd
      double precision rc, r(*), xi(*), xin(*)
      integer i, k
      double precision dr, xn, xo

      k = 0
      xo = 0.d0
      dr = 0.d0
      do i=1,nd-1
1      if (rc.gt.dr) then
        k = k+1
        dr = dr+r(k)
        goto 1
       endif
       if (k.gt.1) xo = xi(k-1)
       xn = xi(k)
       dr = dr-rc
       xin(i) = xn-(xn-xo)*dr/r(k)
      enddo
      do i=1,nd-1
       xi(i) = xin(i)
      enddo
      xi(nd) = 1.d0

      return
      end