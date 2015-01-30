c-------------------------------------------------------------------------------------------------c
      subroutine cmplxconjugate(n,m,mat,output)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Gives the complex conjugate of a complex n x m matrix                                          c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

c input parameters
      integer n,m
      double complex mat(n,m),output(n,m)

c parameters used in this subroutine only
      integer i,j

      do i=1,n
       do j=1,m
        output(i,j) = conjg(mat(i,j))
       enddo
      enddo

      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine dagger(n,m,mat,output)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Gives the Hermitian conjugate of a complex n x m matrix                                        c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

c input parameters
      integer n,m
      double complex mat(n,m),output(m,n)

c parameters used in this subroutine only
      double complex dummy(n,m)

c first get the complex conjugate of the matrix
      call cmplxconjugate(n,m,mat,dummy)

c then transpose the dummy matrix
      call transposecmplx(n,m,dummy,output)

      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine matrixmultiply(n1,m1,n2,m2,mat1,mat2,output)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Calculates the matrix multiplication of matrices m1 and m2 (of dimension r1xc1 and r2xc2       c
c  respectively).  The result is given in matrix output (dim r1xc2).                              c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

c input parameters
      integer n1,m1,n2,m2
      double precision mat1(n1,m1), mat2(n2,m2), output(n1,m2)

c parameters used in this subroutine only
      integer i,j,k

      if(m1.ne.n2) then
       write(*,*) 'Cannot multiply these matrices'
       stop
      endif

c initialize the output matrix
      do i=1,n1
       do j=1,m2
        output(n1,m2) = 0.d0
       enddo
      enddo

      do i=1,n1
       do j=1,m2
        do k=1,m1
         output(i,j) = output(i,j) + mat1(i,k)*mat2(k,j)
        enddo
       enddo
      enddo

      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine matrixmultiplycmplx(n1,m1,n2,m2,mat1,mat2,output)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Calculates the complex matrix multiplication of matrices m1 and m2 (of dimension r1xc1 and     c
c  r2xc2 respectively).  The result is given in matrix output (dim r1xc2).                        c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

c input parameters
      integer n1,m1,n2,m2
      double complex mat1(n1,m1), mat2(n2,m2), output(n1,m2)

c parameters used in this subroutine only
      integer i,j,k

      if(m1.ne.n2) then
       write(*,*) 'Cannot multiply these matrices'
       stop
      endif

c initialize the output matrix
      do i=1,n1
       do j=1,m2
        output(n1,m2) = dcmplx(0.d0,0.d0)
       enddo
      enddo

      do i=1,n1
       do j=1,m2
        do k=1,m1
         output(i,j) = output(i,j) + mat1(i,k)*mat2(k,j)
        enddo
       enddo
      enddo

      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine transposedp(n,m,mat,output)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Transposes a double precision n by m matrix mat into an m by n matrix output                   c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

c input parameters
      integer n,m
      double precision mat(n,m),output(m,n)

c parameters used in this subroutine only
      integer i,j

      do i=1,n
       do j=1,m
        output(j,i) = mat(i,j)
       enddo
      enddo

      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine transposecmplx(n,m,mat,output)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Transposes a double complex n by m matrix mat into an m by n matrix output                     c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

c input parameters
      integer n,m
      double complex mat(n,m),output(m,n)

c parameters used in this subroutine only
      integer i,j

      do i=1,n
       do j=1,m
        output(j,i) = mat(i,j)
       enddo
      enddo

      return
      end
