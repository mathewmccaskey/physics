c-------------------------------------------------------------------------------------------------c
      program entropy
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This program calculates the total entropy density given the mass spectrum of particles.  We    c
c  can use this to compare the g_star quantity that is just written all willy nilly in textbooks  c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'
        
c parameters for this program only
      integer nsmparticles
      double precision s_total, s_integral, g_star, entropy_func
      double precision xinput(600), boson(600), fermion(600)
      double precision masses(100), forb(100), dof(100)
      double precision temp, x_particle, s_particle
      
c common blocks
      double precision x_pass, isfermion
      common/passtemp/ x_pass, isfermion

c external functions
      external entropy_func

c counters
      integer i, j

c setting the parameters
      nsmparticles = 17
        
c SM masses: G, gamma, Z, W, t, b, c, s, d, u, e, mu, tau, nu_e, nu_mu, nu_tau, h
      masses(1) = 0.d0
      masses(2) = 0.d0
      masses(3) = 91.188d0
      masses(4) = 80.4d0
      masses(5) = 125.d0
      masses(6) = 174.3d0
      masses(7) = 4.2d0
      masses(8) = 1.27d0
      masses(9) = 0.104d0
      masses(10) = 4.8d-3
      masses(11) = 2.4d-3
      masses(12) = 0.511d-3
      masses(10) = 0.d0
      masses(11) = 0.d0
      masses(12) = 0.d0
      masses(13) = 0.1057d0
      masses(14) = 1.777d0
      masses(15) = 0.d0
      masses(16) = 0.d0
      masses(17) = 0.d0

c SM degrees of freedom (color charge * spin/polarizations * particle/anti-particle)
      dof(1) = 16.d0
      dof(2) = 2.d0
      dof(3) = 3.d0
      dof(4) = 6.d0
      dof(5) = 1.d0
      dof(6) = 12.d0
      dof(7) = 12.d0
      dof(8) = 12.d0
      dof(9) = 12.d0
      dof(10) = 12.d0
      dof(11) = 12.d0
      dof(12) = 4.d0
      dof(13) = 4.d0
      dof(14) = 4.d0
      dof(15) = 2.d0
      dof(16) = 2.d0
      dof(17) = 2.d0

c setting the isfermion parameter
      do i=1,5
        forb(i) = 0.d0
      enddo
      do i=6,nsmparticles
        forb(i) = 1.d0
      enddo

      open(unit=42, file='boson.txt', status='unknown')
      open(unit=43, file='fermion.txt', status='unknown')

c        read(42,*) xinput(i), boson(i)
c        read(43,*) xinput(i), fermion(i)

      do i=1,401
        xinput(i) = 10**(dble(i-1)/100-2.d0)
        
        x_pass = xinput(i)

        isfermion = 1.d0
        call romberg(entropy_func, 0.d0, 1.d0, s_integral, 1.0d-12, 0, 23)
        if (s_integral.gt.1.d-100) then
          boson(i) = s_integral
        else
          boson(i) = 0.d0
        endif
        write(42,fmt='(2(E16.8))') xinput(i), boson(i)
        
        isfermion = 0.d0
        call romberg(entropy_func, 0.d0, 1.d0, s_integral, 1.0d-12, 0, 23)
        if (s_integral.gt.1.d-100) then
          fermion(i) = s_integral
        else
          fermion(i) = 0.d0
        endif
        write(43,fmt='(2(E16.8))') xinput(i), fermion(i)
        
        write(*,fmt='(3(E16.8))') xinput(i), boson(i), fermion(i)
      enddo      

      close(42)
      close(43)
      
c      do i=1,500
c        temp = 10**(4.d0-dble(i)/100.d0)
c        s_total = 0.d0
c        
c        do j=1,nsmparticles
c          x_particle = masses(j)/temp
c          
c          if (x_particle < 0.01d0) then
c            s_particle = 2.d0*pi**2*dof(j)*temp**3/45.d0 * (7.d0/8.d0)**forb(j)
c
c          else if (x_particle < 1000.d0) then
c            if (forb(j).eq.0.d0) then
c              s_particle = cubeint(x_particle,xinput,boson,500)
c            else if (forb(j).eq.1.d0) then
c              s_particle = cubeint(x_particle,xinput,fermion,500)
c            endif
c
c            s_particle = s_particle*dof(j)*temp**3/(6.d0*pi**2)
c            
c          else
c            s_particle = 0.d0
c          endif
c          
c          s_total = s_total + s_particle
c        enddo
c        
c        g_star = s_total/(2.d0*pi**2/45.d0*temp**3)
c        write(*,fmt='(2(E16.8,1X))') temp, g_star
c
c      enddo

      end



c-------------------------------------------------------------------------------------------------c
      function entropy_func(beta)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This subroutine is the integrand of the entropy density contribution of any particle in the    c
c  model.  The temperature, mass, number of internal degrees of freedom, and whether or not the   c
c  particle is a fermion or boson are passed from the main program so that the function can be    c
c  integrated properly.                                                                           c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

c input parameters
      double precision beta

c function definition
      double precision entropy_func

c common blocks
      double precision x_pass, isfermion
      common/passtemp/ x_pass, isfermion

c just to make sure that we don't hit any of those beta = 1 singularities
      if (1.d0-beta**2 .eq. 0.d0) then
        entropy_func = 0.d0
        return
      endif
      
      entropy_func = x_pass**(4.d0)*beta**2*(3.d0+beta**2)/(1.d0-beta**2)**3
     .      /(dexp(x_pass/dsqrt(1.d0-beta**2)) + (-1.d0)**(isfermion))
     .      *15.d0/(4.d0*pi**4)
     
      return
      end