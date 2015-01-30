c-------------------------------------------------------------------------------------------------c
      function mt(n,p_vis,m_invisible,pT_missing)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This function calculates mt2 for n number of visible particles (with corresponding 4-momenta   c
c  given in p_vis) and the invisible particles mass, m_inv, and missing pT.                       c
c                                                                                                 c
c  Generally for mT the idea is to sum up all the visible particles before projecting them into   c
c  the transverse plane.  Hence ideally we should only have n = 1.                                c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'
      
c input parameters
      integer n
      double precision p_vis(0:3,nmax), m_invisible, pt_missing(0:3)

c parameters used in this subroutine only
      double precision E_total, p_dummy(0:3), pt_total(0:3)
      
c initializing variables
      E_total = 0.d0
      do x1=0,3
        pt_total(x1) = 0.d0
      enddo

c gets all the E_t and p_t information from the visible particles      
      do x1=1,n
        do x2=0,3
          p_dummy(x2) = p_vis(x2,x1)
        enddo
        E_total = E_total + Et(p_dummy)
        pt_total(1) = pt_total(1) + p_vis(1,x1)
        pt_total(2) = pt_total(2) + p_vis(2,x1)
      enddo
      
c same thing for the invisible particles
      E_total = E_total + dsqrt(m_invisible**2 + p3dot(pt_missing,pt_missing))
      pt_total(1) = pt_total(1) + pt_missing(1)
      pt_total(2) = pt_total(2) + pt_missing(2)
      
      mt = dsqrt(E_total**2 - p3dot(pt_total,pt_total))
      
      return
      end
      
      
      
c-------------------------------------------------------------------------------------------------c
      function mt2()
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This function calculates the value of mt2 given two sets of visible particle momenta           c
c  (p_vis1, p_vis2), the mass of the invisible particles (mX), and the total missing pt.          c
c                                                                                                 c
c  Common block for parameters needed for mt2 in matlib.inc                                       c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'
      
c parameters used in this subroutine only
      double precision points(3,2), x_val(2), y_val(3), ftol
      parameter (ftol=1d-4)
      integer iter
      double precision mt2_func
      external mt2_func

c setting up 3 initial points with huge pt's for mt2 calculation
      points(1,1) = 1.d4
      points(1,2) = 1.d4
      points(2,1) = -1.d4
      points(2,2) = 1.d4
      points(3,1) = 0.d0
      points(3,2) = -1.d4

c calculating the initial values of y used for amoeba      
      do x1=1,3
        do x2=1,2
          x_val(x2) = points(x1,x2)
        enddo
        y_val(x1) = mt2_func(x_val)
      enddo
      
c call to the minimization subroutine amoeba
      call amoeba(points,y_val,3,2,2,ftol,mt2_func,iter)
      
      mt2 = y_val(1)
      
      return
      end
      
      
      
c-------------------------------------------------------------------------------------------------c
      function mt2_func(x)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This is the function that is to be minimized in order to calculate mt2.                        c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'
      
c input parameters
      double precision x(2), mt2_func
      
c parameters used in this subroutine only
      double precision p_inv1(0:3), p_inv2(0:3), mt_1, mt_2

c initializing missing momentum arrays      
      do x1=0,3
        p_inv1(x1) = 0.d0
        p_inv2(x1) = 0.d0
      enddo

c setting up the missing pt arrays to be consistent with the total missing pt (pt_miss)      
      do x1=1,2
        p_inv1(x1) = x(x1)
        p_inv2(x1) = pt_miss(x1) - x(x1)
      enddo
      
      mt_1 = mt(n_vis(1),p_vis1,m_inv(1),p_inv1)
      mt_2 = mt(n_vis(2),p_vis2,m_inv(2),p_inv2)
      
      mt2_func = max(mt_1,mt_2)
      
      return
      end
