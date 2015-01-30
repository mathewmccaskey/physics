c-------------------------------------------------------------------------------------------------c
      function cubeint(xvalue,xarray,yarray,n)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This function finds the cubic interpolation given an input value xvalue and arrays xarray and  c
c  yarray (both arrays are length n).                                                             c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'
      
c input parameters
      integer n
      double precision xvalue, xarray(n), yarray(n)
      
c parameters used in this subroutine only
      integer index, imin, i, j
      double precision coeff(4)
      
c find the index of the array we need
      index = getindex(xvalue, xarray, n)
      
c set the minimum of the four indices needed for the cubic interpolation
      if (index .eq. 1) then
        imin = 1
      else if (index .eq. n-1) then
        imin = n-3
      else
        imin = index-1
      endif
      
c set up the coefficients for the cubic interpolation
      do i=1,4
        coeff(i) = yarray(imin-1+i)
        do j=1,4
          if (i.ne.j) then
            coeff(i) = coeff(i)/(xarray(imin-1+i)-xarray(imin-1+j))
          endif
        enddo
      enddo

c and now to calculate the cubic interpolation
      cubeint = 0.d0
      do i=1,4 
        do j=1,4
          if (i.ne.j) then
            coeff(i) = coeff(i)*(xvalue-xarray(imin-1+j))
          endif
        enddo
        cubeint = cubeint + coeff(i)
      enddo
      
      return
      end
      
      
      
c-------------------------------------------------------------------------------------------------c
      function getindex(xvalue,xarray,n)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This function returns the index of the array, i, such that xarray(i) < xvalue < xarray(i+1).   c
c  If xvalue is outside of the array then getindex returns the boundary.                          c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

c input parameters
      integer n
      double precision xvalue, xarray(n)

c parameters used in this subroutine only      
      integer lowindex, highindex, mid
      double precision low, high
      
c If the xvalue is outside the given array then set the index to the boundary
      if (xvalue .lt. xarray(1)) then
        getindex = 1
        return
      else if (xvalue .gt. xarray(n)) then
        getindex = n - 1
        return
      endif
      
c Binary search
      lowindex = 1
      highindex = n

      do 100 while ((highindex-lowindex).gt.1)
        mid = (highindex+lowindex)/2
        if (xarray(mid).le.xvalue) lowindex = mid
        if (xarray(mid).gt.xvalue) highindex = mid
 100  continue

c Now we have our index
      getindex = lowindex

      return
      end



c-------------------------------------------------------------------------------------------------c
      function lineint(x,xarray,farray,nx)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This function finds the linear interpolation yvalue given an input value x and arrays xarray   c
c  and yarray (both arrays are length nx).                                                        c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'
      
c input parameters
      integer nx
      double precision x, xarray(nx), farray(nx)

c paramters used in this subroutine only
      integer xindex
      double precision x_1, x_2
      
c find the index of the array we need    
      xindex = getindex(x, xarray, nx)

      x_1 = xarray(xindex)
      x_2 = xarray(xindex+1)

c interpolate the array with a line      
      lineint = farray(xindex  )*((x_2-x)/(x_2-x_1))
     .        + farray(xindex+1)*((x-x_1)/(x_2-x_1))
      
      return
      end
      
      
      
c-------------------------------------------------------------------------------------------------c
      function lineint2D(x,y,xarray,yarray,farray,nx,ny)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This function finds the 2D linear interpolation lineint2D given the xy coordinates             c
c  (xvalue,yvalue), arrays of the x and y grid (1,...,nx) and (1,...,ny) respectively, and the    c
c  values of the function in farray (1,...,nx)(1,...,ny).                                         c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'
      
c input parameters
      integer nx, ny
      double precision x, y, xarray(nx), yarray(ny), farray(nx,ny)

c parameters used in this subroutine only
      integer xindex, yindex
      double precision x_1, x_2, y_1, y_2
      
      xindex = getindex(x,xarray,nx)
      yindex = getindex(y,yarray,ny)
      
      x_1 = xarray(xindex)
      x_2 = xarray(xindex+1)
      y_1 = yarray(yindex)
      y_2 = yarray(yindex+1)
      
      lineint2D = farray(xindex  ,yindex  )*((x_2-x)/(x_2-x_1))*((y_2-y)/(y_2-y_1))
     .          + farray(xindex  ,yindex+1)*((x_2-x)/(x_2-x_1))*((y-y_1)/(y_2-y_1))
     .          + farray(xindex+1,yindex  )*((x-x_1)/(x_2-x_1))*((y_2-y)/(y_2-y_1))
     .          + farray(xindex+1,yindex+1)*((x-x_1)/(x_2-x_1))*((y-y_1)/(y_2-y_1))

      return
      end
    
    
    
c-------------------------------------------------------------------------------------------------c
      function lineint3D(x,y,z,xarray,yarray,zarray,farray,nx,ny,nz)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This function finds the 3D linear interpolation lineint3D given the xyz coordinates (x,y,z),   c
c  arrays of the x, y, and z grid (1,..,nx), (1,..,ny), and (1,..,nz) respectively, and the       c
c  values of the function in farray (1,..,nx)(1,..,ny)(1,..,nz).                                  c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'
      
c input parameters
      integer nx, ny, nz
      double precision x, y, z, xarray(nx), yarray(ny), zarray(nz), farray(nx,ny,nz)
      
c parameters used in this subroutine only
      integer xindex, yindex, zindex
      double precision x_1, x_2, y_1, y_2, z_1, z_2
      
      xindex = getindex(x,xarray,nx)
      yindex = getindex(y,yarray,ny)
      zindex = getindex(z,zarray,nz)
      
      x_1 = xarray(xindex)
      x_2 = xarray(xindex+1)
      y_1 = yarray(yindex)
      y_2 = yarray(yindex+1)
      z_1 = zarray(zindex)
      z_2 = zarray(zindex+1)
      
      lineint3D = farray(xindex  ,yindex  ,zindex  )*((x_2-x)/(x_2-x_1))*((y_2-y)/(y_2-y_1))*((z_2-z)/(z_2-z_1))
     .          + farray(xindex  ,yindex  ,zindex+1)*((x_2-x)/(x_2-x_1))*((y_2-y)/(y_2-y_1))*((z-z_1)/(z_2-z_1))
     .          + farray(xindex  ,yindex+1,zindex  )*((x_2-x)/(x_2-x_1))*((y-y_1)/(y_2-y_1))*((z_2-z)/(z_2-z_1))
     .          + farray(xindex  ,yindex+1,zindex+1)*((x_2-x)/(x_2-x_1))*((y-y_1)/(y_2-y_1))*((z-z_1)/(z_2-z_1))
     .          + farray(xindex+1,yindex  ,zindex  )*((x-x_1)/(x_2-x_1))*((y_2-y)/(y_2-y_1))*((z_2-z)/(z_2-z_1))
     .          + farray(xindex+1,yindex  ,zindex+1)*((x-x_1)/(x_2-x_1))*((y_2-y)/(y_2-y_1))*((z-z_1)/(z_2-z_1))
     .          + farray(xindex+1,yindex+1,zindex  )*((x-x_1)/(x_2-x_1))*((y-y_1)/(y_2-y_1))*((z_2-z)/(z_2-z_1))
     .          + farray(xindex+1,yindex+1,zindex+1)*((x-x_1)/(x_2-x_1))*((y-y_1)/(y_2-y_1))*((z-z_1)/(z_2-z_1))
            
      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine readinfile(filename,n,xarray,yarray)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This subroutine is designed to make it very easy to read in a file consisting of two columns   c
c  of double precision data.  In conjunction with lineint or cubeint we can interpolate these     c
c  curves at any point.                                                                           c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none

      integer n, i
      double precision xarray(n), yarray(n) 
      character*100 filename

      open(unit=10,file=filename,status='old')
      
      do i=1,n
       read(10,*) xarray(i),yarray(i)
      enddo

      close(10)

      return
      end