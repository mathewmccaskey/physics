c-------------------------------------------------------------------------------------------------c
      function besselk(alpha,x,eps)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This subroutine calculates the value of the modified bessel function of the second kind        c
c  K_{\alpha}(x) given the user supplied values of alpha, x, and the required precision.          c
c                                                                                                 c
c  The formula comes from the infinite series expansion taken from Wikipedia.                     c
c                                                                                                 c
c  K_{\alpha}(x) = \sqrt{\frac{\pi}{2x}}e^{-x}(1+\frac{4\alpha^{2}-1}{8x}                         c
c    + \frac{(4\alpha^{2}-1)(4\alpha^{2}-9)}{2!(8x)^{2}}                                          c
c    + \frac{(4\alpha^{2}-1)(4\alpha^{2}-9)(4\alpha^{2}-25)}{3!(8x)^{3}} + ...)                   c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

c function definitions and input parameters
      double precision besselk, besselk_factor
      double precision alpha, x, eps
      
c parameters for this function only
      double precision pre_factor, sum_factor, pi
      parameter (pi=4.d0*datan(1.d0))
      
      pre_factor = dsqrt(pi/(2.d0*x))*dexp(-x)
      sum_factor = besselk_factor(alpha,x,eps)
            
      besselk = pre_factor*sum_factor
      
      return
      end



c-------------------------------------------------------------------------------------------------c
      function besselk_factor(alpha,x,eps)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This subroutine calculates the infinite sum associated with the modified bessel function of    c
c  the second kind given the user supplies values of alpha, x, and the required precision.        c
c                                                                                                 c
c  K_{\alpha}(x) sum = (1+\frac{4\alpha^{2}-1}{8x}                                                c
c    + \frac{(4\alpha^{2}-1)(4\alpha^{2}-9)}{2!(8x)^{2}}                                          c
c    + \frac{(4\alpha^{2}-1)(4\alpha^{2}-9)(4\alpha^{2}-25)}{3!(8x)^{3}} + ...)                   c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

c function definition and input parameters
      double precision besselk_factor
      double precision alpha, x, eps

c parameters for this function only
      integer i, n, factorial
      double precision result, result_old, factor

      result = 1.d0
      result_old = 0.d0
      n = 0
      
      do while(dabs((result-result_old)/result).ge.eps)
        result_old = result
        n = n + 1
        
        factor = 1.d0
        do i=1,n
          factor = factor*(4.d0*alpha**2-i**2)/(dble(factorial(i))*(8.d0*x)**i)
        enddo
        
        result = result + factor
      enddo
      
      besselk_factor = result
      
      return
      end