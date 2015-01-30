c-------------------------------------------------------------------------------------------------c
      subroutine setgamma(gamma0,gamma1,gamma2,gamma3)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This subroutine sets the 4 complex gamma matrices. This subourtine is intended to be called    c
c  at the beginning of any program so that the gamma matrices can be used anywhere else in the    c
c  program.                                                                                       c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

c gamma^0
      gamma0(0,0) = dcmplx(0.d0,0.d0)
      gamma0(0,1) = dcmplx(0.d0,0.d0)
      gamma0(0,2) = dcmplx(1.d0,0.d0)
      gamma0(0,3) = dcmplx(0.d0,0.d0)
      gamma0(1,0) = dcmplx(0.d0,0.d0)
      gamma0(1,1) = dcmplx(0.d0,0.d0)
      gamma0(1,2) = dcmplx(0.d0,0.d0)
      gamma0(1,3) = dcmplx(1.d0,0.d0)
      gamma0(2,0) = dcmplx(1.d0,0.d0)
      gamma0(2,1) = dcmplx(0.d0,0.d0)
      gamma0(2,2) = dcmplx(0.d0,0.d0)
      gamma0(2,3) = dcmplx(0.d0,0.d0)
      gamma0(3,0) = dcmplx(0.d0,0.d0)
      gamma0(3,1) = dcmplx(1.d0,0.d0)
      gamma0(3,2) = dcmplx(0.d0,0.d0)
      gamma0(3,3) = dcmplx(0.d0,0.d0)

c gamma^1
      gamma1(0,0) = dcmplx(0.d0,0.d0)
      gamma1(0,1) = dcmplx(0.d0,0.d0)
      gamma1(0,2) = dcmplx(0.d0,0.d0)
      gamma1(0,3) = dcmplx(1.d0,0.d0)
      gamma1(1,0) = dcmplx(0.d0,0.d0)
      gamma1(1,1) = dcmplx(0.d0,0.d0)
      gamma1(1,2) = dcmplx(1.d0,0.d0)
      gamma1(1,3) = dcmplx(0.d0,0.d0)
      gamma1(2,0) = dcmplx(0.d0,0.d0)
      gamma1(2,1) = dcmplx(-1.d0,0.d0)
      gamma1(2,2) = dcmplx(0.d0,0.d0)
      gamma1(2,3) = dcmplx(0.d0,0.d0)
      gamma1(3,0) = dcmplx(-1.d0,0.d0)
      gamma1(3,1) = dcmplx(0.d0,0.d0)
      gamma1(3,2) = dcmplx(0.d0,0.d0)
      gamma1(3,3) = dcmplx(0.d0,0.d0)

c gamma^2
      gamma2(0,0) = dcmplx(0.d0,0.d0)
      gamma2(0,1) = dcmplx(0.d0,0.d0)
      gamma2(0,2) = dcmplx(0.d0,0.d0)
      gamma2(0,3) = dcmplx(0.d0,-1.d0)
      gamma2(1,0) = dcmplx(0.d0,0.d0)
      gamma2(1,1) = dcmplx(0.d0,0.d0)
      gamma2(1,2) = dcmplx(0.d0,1.d0)
      gamma2(1,3) = dcmplx(0.d0,0.d0)
      gamma2(2,0) = dcmplx(0.d0,0.d0)
      gamma2(2,1) = dcmplx(0.d0,1.d0)
      gamma2(2,2) = dcmplx(0.d0,0.d0)
      gamma2(2,3) = dcmplx(0.d0,0.d0)
      gamma2(3,0) = dcmplx(0.d0,-1.d0)
      gamma2(3,1) = dcmplx(0.d0,0.d0)
      gamma2(3,2) = dcmplx(0.d0,0.d0)
      gamma2(3,3) = dcmplx(0.d0,0.d0)

c gamma^3
      gamma3(0,0) = dcmplx(0.d0,0.d0)
      gamma3(0,1) = dcmplx(0.d0,0.d0)
      gamma3(0,2) = dcmplx(1.d0,0.d0)
      gamma3(0,3) = dcmplx(0.d0,0.d0)
      gamma3(1,0) = dcmplx(0.d0,0.d0)
      gamma3(1,1) = dcmplx(0.d0,0.d0)
      gamma3(1,2) = dcmplx(0.d0,0.d0)
      gamma3(1,3) = dcmplx(-1.d0,0.d0)
      gamma3(2,0) = dcmplx(-1.d0,0.d0)
      gamma3(2,1) = dcmplx(0.d0,0.d0)
      gamma3(2,2) = dcmplx(0.d0,0.d0)
      gamma3(2,3) = dcmplx(0.d0,0.d0)
      gamma3(3,0) = dcmplx(0.d0,0.d0)
      gamma3(3,1) = dcmplx(1.d0,0.d0)
      gamma3(3,2) = dcmplx(0.d0,0.d0)
      gamma3(3,3) = dcmplx(0.d0,0.d0)

      return
      end