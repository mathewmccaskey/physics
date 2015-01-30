c-------------------------------------------------------------------------------------------------c
      subroutine eigsrt(d,v,n,np,option)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Given the eigenvalues d and eigenvectors v as output from jacobi this routine sorts the        c
c  eigenvalues into descending order, and rearranges the columns of v correspondingly.  The       c
c  method is straight insertion                                                                   c
c                                                                                                 c
c  option: 1-eigenvalues in descending order                                                      c
c          2-eigenvalues in ascending order                                                       c
c-------------------------------------------------------------------------------------------------c
      implicit none

      integer n,np,option
      double precision d(np),v(np,np)
      integer i,j,k
      double precision p

      do i=1,n-1
       k = i
       p = d(i)
       do j=i+1,n
        if(option.eq.1) then
         if(d(j).ge.p) then
          k = j
          p = d(j)
         endif
        else if(option.eq.2) then
         if(d(j).le.p) then
          k = j
          p = d(j)
         endif
        endif
       enddo
       if(k.ne.i) then
        d(k) = d(i)
        d(i) = p
        do j=1,n
         p = v(j,i)
         v(j,i) = v(j,k)
         v(j,k) = p
        enddo
       endif
      enddo
      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine eigsrt2(n,eigen,sort,type)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Sorts out the eigenvalues in a different way than the numerical recipe eigensort.  This        c
c  allows a clearer way to see what is the composition of each eigenvalue in terms of the         c
c  original mixing matrix.                                                                        c
c                                                                                                 c
c  n - the number of eigenvalues                                                                  c
c  eigen - n-dimensional array of eigenvalues in double precision                                 c
c  sort - n-dimensional array of integers that lists the index in order of increasing eigenvalue  c
c  type - 1 for ascending order 2 for descending order.                                           c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

c input parameters
      integer n,sort(n),type
      double precision eigen(n)

c parameters used in this subroutine only
      integer i,j,k,tmp
      double precision p

      if ((type.gt.2).or.(type.lt.1)) then
        write(*,*) 'sorting type not chosen'
        return
      endif

      do i=1,n
       sort(i) = i
      enddo

      do i=1,n-1
       k = i
       p = eigen(sort(k))
       do j=i+1,n
        if (((eigen(sort(j)).le.p).and.(type.eq.1)).or.((eigen(sort(j)).ge.p).and.(type.eq.2))) then
         k = j
         p = eigen(sort(k))
        endif
       enddo
       if (k.ne.i) then
        tmp = sort(i)
        sort(i) = sort(k)
        sort(k) = tmp
       endif
      enddo

      return
      end