c-------------------------------------------------------------------------------------------------c
      program adiabatic
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This program finds the adiabatic relationship between temperature and volume.                  c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
        
c parameters for this program only
      double precision entropy_func, xinput, T(601), entropy(601)
      double precision T_log(601), entropy_log(601)
      double precision T_slope(600), entropy_slope(600)
      integer iter
      
c common blocks
      double precision m_pass, T_pass
      common/passvars/ m_pass, T_pass

c external functions
      external entropy_func

c counters
      integer i, j

c Get the mass from the user
      write(*,*) 'Enter particle mass'
      read(*,*) m_pass

      open(unit=42,file='entropy_v_temp.txt')
      do i=1,601

        xinput = 10**(dble(601-i)/100.d0-3.d0)
        T(i) = m_pass/xinput
        T_pass = T(i)
        call romberg(entropy_func, 0.d0, 1.d0, entropy(i), 1.0d-12, 20)
        write(42,fmt='(2(SE16.8))') T(i), entropy(i)
        write(*,*) T(i), entropy(i)
        
      enddo      
      close(42)

c Now take the log of both T and the entropy
      open(unit=42,file='log_entropy_v_log_temp.txt')
      do i=1,601
        T_log(i) = dlog(T(i))
        entropy_log(i) = dlog(entropy(i))
        write(42,fmt='(2(SE16.8))') T_log(i), entropy_log(i)
      enddo
      close(42)
      
c Now take the slope of log(entropy) v. log(T)
      open(unit=42,file='slope.txt')
      do i=1,600
        T_slope(i) = (T(i+1)+T(i))/2.d0
        entropy_slope = (entropy_log(i+1)-entropy_log(i))/(T_log(i+1)-T_log(i))
        write(42,fmt='(2(SE16.8))') T_slope(i), entropy_slope(i)
      enddo
      close(42)

      end



c-------------------------------------------------------------------------------------------------c
      function entropy_func(beta)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This function gives the integrand of the entropy density for any given mass and termperature.  c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

c input parameters
      double precision beta

c function definition
      double precision entropy_func

c common blocks
      double precision m_pass, T_pass
      common/passvars/ m_pass, T_pass

c just to make sure that we don't hit any of those beta = 1 singularities
      if (1.d0-beta**2 .eq. 0.d0) then
        entropy_func = 0.d0
        return
      endif

C c Entropy function that included the rest mass energy in the energy density
C       entropy_func = m_pass**4/(6.d0*3.14159265359d0**2*T_pass)
C      .      *beta**2*(3.d0+beta**2)/(1.d0-beta**2)**3
C      .      /(1.d0+dexp((m_pass/T_pass)*(1.d0/dsqrt(1.d0-beta**2)-1.d0)))

c Entropy function that does not include the rest mass energy in the energy density
      entropy_func = m_pass**4/(6.d0*3.14159265359d0**2*T_pass)
     .      *beta**2*(3.d0*(1.d0-dsqrt(1.d0-beta**2))+beta**2)/(1.d0-beta**2)**3
     .      /(1.d0+dexp((m_pass/T_pass)*(1.d0/dsqrt(1.d0-beta**2)-1.d0)))
  
      return
      end