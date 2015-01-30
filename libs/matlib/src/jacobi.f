c-------------------------------------------------------------------------------------------------c
      subroutine jacobi(a,n,np,d,v,nrot)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Computes all eigenvalues and eigenvectors of a real symmetric matrix a, which is of size       c
c  n by n, stored in a physical np by np array.  On ouput, elements of a above the diagonal are   c
c  destroyed. d returns the eigenvalues of a in its first n elements.  v is a matrix with the     c
c  same logical and physical dimensions as a, whose columns contain, on output, the normalized    c
c  eigenvectors of a.  nrot reutrns the number of Jacobi rotations that were required.            c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

c Input parameters
      integer n,np,nrot,nmax
      double precision a(np,np),d(np),v(np,np)

c Parameters appearing in this subroutine only
      parameter (nmax=500)
      integer i,ip,iq,j
      double precision c,g,h,s,sm,t,tau,theta,tresh
      double precision b(nmax),z(nmax)

c Initialize to the identity matrix
      do ip=1,n
       do iq=1,n
        v(ip,iq) = 0.d0
       enddo
       v(ip,ip) = 1.d0
      enddo

c Initialize b and d to the diagonal of a
      do ip=1,n
       b(ip) = a(ip,ip)
       d(ip) = b(ip)
       z(ip) = 0.d0
      enddo
      nrot = 0
      do i=1,50

c Sum of off-diagonal elements
       sm = 0.d0
       do ip=1,n-1
        do iq=ip+1,n
         sm = sm+dabs(a(ip,iq))
        enddo
       enddo

c The normal return, which relies on quadratic convergence to
c machine underflow
       if(sm.eq.0.d0) return

c .. on the first three sweeps
       if(i.lt.4) then
        tresh = 0.2d0*sm/n**2

c .. thereafter
       else
        tresh = 0.d0
       endif
       do ip=1,n-1
        do iq=ip+1,n
         g = 100.d0*dabs(a(ip,iq))

c After four sweeps, skip the rotation if the off-diagonal
c element is small.
         if((i.gt.4).and.(dabs(d(ip))+g.eq.dabs(d(ip)))
     .     .and.(dabs(d(iq))+g.eq.dabs(d(iq)))) then
          a(ip,iq) = 0.d0
         else if(dabs(a(ip,iq)).gt.tresh) then
          h = d(iq)-d(ip)
          if (dabs(h)+g.eq.dabs(h)) then
           t = a(ip,iq)/h
          else
           theta = 0.5d0*h/a(ip,iq)
           t = 1.d0/(dabs(theta)+dsqrt(1.d0+theta**2))
           if(theta.lt.0.d0) t=-t
          endif
          c = 1.d0/dsqrt(1.d0+t**2)
          s = t*c
          tau = s/(1.d0+c)
          h = t*a(ip,iq)
          z(ip) = z(ip)-h
          z(iq) = z(iq)+h
          d(ip) = d(ip)-h
          d(iq) = d(iq)+h
          a(ip,iq) = 0.d0

c Case of rotations 1 <= j < p
          do j=1,ip-1
           g = a(j,ip)
           h = a(j,iq)
           a(j,ip) = g-s*(h+g*tau)
           a(j,iq) = h+s*(g-h*tau)
          enddo

c Case of rotations p < j < q
          do j=ip+1,iq-1
           g = a(ip,j)
           h = a(j,iq)
           a(ip,j) = g-s*(h+g*tau)
           a(j,iq) = h+s*(g-h*tau)
          enddo

c Case of rotations q < j <= n
          do j=iq+1,n
           g = a(ip,j)
           h = a(iq,j) 
           a(ip,j) = g-s*(h+g*tau)
           a(iq,j) = h+s*(g-h*tau)
          enddo
          do j=1,n
           g = v(j,ip)
           h = v(j,iq)
           v(j,ip) = g-s*(h+g*tau)
           v(j,iq) = h+s*(g-h*tau)
          enddo
          nrot = nrot+1
         endif
        enddo
       enddo
       do ip=1,n
        b(ip) = b(ip)+z(ip)

c Update d with the sum of ta_pq and reinitialize z
        d(ip) = b(ip)
        z(ip) = 0.d0
       enddo
      enddo
      pause 'too many iterations in jacobi'
      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine jacobicmplx(a,n,np,d,v,nrot)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Computes all eigenvalues and eigenvectors of a hermitian matrix of size n x n stored in a      c
c  physical size np x np.  The complex matrix is broken into real and imaginary parts and the     c
c  2n x 2n problem is fed into jacobi.                                                            c
c-------------------------------------------------------------------------------------------------c
      implicit none

c Input parameters
      integer n,np,nrot
      double complex a(np,np),d(np),v(np,np)

c Parameters appearing in this subroutine only
      double precision matreal(n,n), matimag(n,n)
      double precision aug(2*n,2*n),daug(2*n), vaug(2*n,2*n)
      double complex minusi
      integer i,j,k

c Complex number used for extracting the imaginary part of a c-number
      minusi = dcmplx(0.d0,-1.d0)

c Setting up the real and imaginary martrices matreal and matimag
      do i=1,n
       do j=1,n
        matreal(i,j) = dble(a(i,j)+conjg(a(i,j)))/2.d0
        matimag(i,j) = dble(minusi*(a(i,j)-conjg(a(i,j))))/2.d0
       enddo
      enddo

c Setting up the augmented matrix aug
      do i=1,n
       do j=1,n
        aug(i,j) = matreal(i,j)
        aug(i,j+2) = -1.d0*matimag(i,j)
        aug(i+2,j) = matimag(i,j)
        aug(i+2,j+2) = matreal(i,j)
       enddo
      enddo

c Calling jacobi with the augmented matrix
      call jacobi(aug,2*n,2*n,daug,vaug,nrot)

c Just in case jacobi orders the eigenvalues strangely
      call eigsrt(daug,vaug,2*n,2*n)

      do i=1,n
       d(i) = daug(2*i-1)
       do j=1,n
        v(j,i) = dcmplx(vaug(j,2*i-1),vaug(j+n,2*i-1))
       enddo
      enddo

      return
      end
