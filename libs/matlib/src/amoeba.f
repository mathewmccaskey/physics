c-------------------------------------------------------------------------------------------------c
      subroutine amoeba(p,y,mp,np,ndim,ftol,func,iter)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Multidimensional minimization of the function func(x) where x(1:ndim) is a vector in ndim      c
c  dimensions, by the downhill simplex method of Nelder and Mead.  The matrix p(1:ndim+1,1:ndim)  c
c  is input.  Its ndim+1 rows are ndim-dimensional vectors which are the vertices of the          c
c  starting simplex.  Also input is the vector y(1:ndim+1) whose components must be               c
c  pre-initialized to the values of func evaluated at the ndim+1 vertices (rows) of p: and ftol   c
c  the fractional convergence tolerance to be achieved in the function value (n.b.!).  On         c
c  output, p and y will have been reset to ndom+1 new points all within ftol of a minimum         c
c  function value, and iter gives the number of function evaluations taken.                       c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
       
c input parameters
      integer iter,mp,ndim,np,nmax,itmax
      double precision ftol,p(mp,np),y(mp),func,tiny
      parameter (nmax=20,itmax=10000,tiny=1.d-10)
      external func
c uses amotry,func

c parameters for this subroutine only
      integer i,ihi,ilo,inhi,j,m,n
      double precision rtol,sum,swap,ysave,ytry,psum(nmax),amotry
      iter = 0
       
c Enter here when starting or have just overall contracted
 1    do n=1,ndim

c Recompute psum
        sum = 0.d0
        do m=1,ndim+1
          sum = sum + p(m,n)
        enddo
        psum(n) = sum
      enddo
      
c Enter here when have just changed a single point.
 2    ilo = 1
 
c Determine which point is the highest (worst). next-highest,
c and lowest (best)
      if (y(1).gt.y(2)) then
        ihi = 1
        inhi = 2
      else
        ihi = 2
        inhi = 1
      endif
      do i=1,ndim+1
        if (y(i).le.y(ilo)) ilo = i
        if (y(i).gt.y(ihi)) then
          inhi = ihi
          ihi = i
        else if(y(i).gt.y(inhi)) then
          if(i.ne.ihi) inhi = i
        endif
      enddo

c Compute the fractional range from highest to lowest and return if satisfactory
      rtol = 2.d0*dabs(y(ihi)-y(ilo))/(dabs(y(ihi))+dabs(y(ilo))+tiny)
      if (rtol.lt.ftol) then
        
c If returning, put the best point and value in slot 1.
        swap = y(1)
        y(1) = y(ilo)
        y(ilo) = swap
        do n=1,ndim
          swap = p(1,n)
          p(1,n) = p(ilo,n)
          p(ilo,n) = swap
        enddo
        return
      endif
      if (iter.ge.itmax) then
        pause 'itmax exceeded in amoeba'
        return
      endif

c Begin a new iteration. First extrapolate by a factor -1 through the face of
c the simplex across from the high point, i.e. reflect the simplex from the 
c high point
      iter = iter + 2
      ytry = amotry(p,y,psum,mp,np,ndim,func,ihi,-1.d0)

c Gives a result better than the best point, so try an additional extrapolation
c by a factor of 2.
      if (ytry.le.y(ilo)) then
        ytry = amotry(p,y,psum,mp,np,ndim,func,ihi,2.d0)

c The reflected point is worse than the second-highest, so look for an intermediate
c lower point, i.e., do a one-dimensional contraction
      else if (ytry.ge.y(inhi)) then
        ysave = y(ihi)
        ytry = amotry(p,y,psum,mp,np,ndim,func,ihi,0.5d0)
        
c Can't seem to get rid of that high point. Better contract
c around the lowest (best) point
        if (ytry.ge.ysave) then
          do i=1,ndim+1
            if(i.ne.ilo) then
              do j=1,ndim
                psum(j) = 0.5d0*(p(i,j)+p(ilo,j))
                p(i,j) = psum(j)
              enddo
              y(i) = func(psum)
            endif
          enddo
          iter = iter + ndim
          goto 1
        endif
      else
        iter = iter - 1
      endif
      goto 2

      END



c-------------------------------------------------------------------------------------------------c
      function amotry(p,y,psum,mp,np,ndim,func,ihi,fac)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Extrapolates by a factor fac through the face of the simplex across from the high point,       c
c  tried it, and replaces the high point if the new point is better.                              c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

c input parameters
      integer ihi,mp,ndim,np,nmax
      double precision amotry,fac,p(mp,np),psum(np),y(mp),func
      parameter (nmax=20)
      external func
      
c parameters used in this subroutine only
      integer i
      double precision fac1,fac2,ytry,ptry(nmax)
      fac1 = (1.d0-fac)/ndim
      fac2 = fac1 - fac
      do i=1,ndim
        ptry(i) = psum(i)*fac1-p(ihi,i)*fac2
      enddo

c Evaluate the function at the trial point.
c If it's better than the highest, then replace the highest.
      ytry = func(ptry)
      if (ytry.lt.y(ihi)) then
        y(ihi) = ytry
        do i=1,ndim
          psum(i) = psum(i) - p(ihi,i) + ptry(i)
          p(ihi,i) = ptry(i)
        enddo
      endif
      amotry = ytry
      return

      end