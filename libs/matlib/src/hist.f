c-------------------------------------------------------------------------------------------------c
      subroutine histinit()
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  Clears out all the histogram information and sets data to some default information that will   c
c  make it easy to see if something has gone horribly wrong.                                      c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

      nhist = 0
      do x1=1,nhistmax
        nbins(x1) = 0
        label(x1) = 'Error'
        xmins(x1) = -1.d0
        xmaxs(x1) = -1.d0
        do x2=1,nbinmax
          hist(x1,x2) = -1.d0
        enddo
      enddo

      return
      end
              


c-------------------------------------------------------------------------------------------------c
      subroutine hbook(id,inlabel,nx,xmin,xmax,zinitial)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This subroutine initializes a one independent variable histogram shamelessly based on similar  c
c  subroutines from hanlib.                                                                       c
c                                                                                                 c
c  id = output integer used to identify histogram                                                 c
c  inlabel = label to be written on the output by the plotting program (character of len <= 30)   c
c  nx = number of x bins                                                                          c
c  xmin = min x value                                                                             c
c  xmax = max x value                                                                             c
c  zinitial = initial value for each bin                                                          c
c                                                                                                 c
c  Internal common blocks for the hbook routines:                                                 c
c                                                                                                 c
c  labels(i) = label for histogram i                                                              c
c  nhist = number of histograms (starts at 0, max is 200)                                         c
c  idnumber(i) = code number to identify histogram in hfill                                       c
c  hist(i,j) = data for histogram i, bin j                                                        c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

c input parameters
      integer id, nx
      character* (*) inlabel
      double precision xmin, xmax, zinitial      

      if (nhist .eq. nhistmax) then
        write(*,*) 'Maximum number of histograms exceeded'
      else
        nhist = nhist+1
        id = nhist
        label(id) = inlabel

        if (nx .gt. nbinmax) then
          write(*,*) 'Number of bins too large: setting to nbinmax'
          nbins(id) = nbinmax
        else
          nbins(id) = nx
        endif

        xmins(id) = xmin
        xmaxs(id) = xmax
        do x1=1,nx
          hist(id,x1) = zinitial
        enddo
      endif

      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine hfill(id,x,zincrement)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This subroutine adds zincrement to the bin that contains the value x in histogram id.          c
c                                                                                                 c
c  id - identification number of the hitogram                                                     c
c  x - value to find out which bit to fill                                                        c
c  zincrement - value to fill in the bin                                                          c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

c input parameters
      integer id
      double precision x, zincrement

c parameters for this subroutine only
      integer lowindex, highindex, mid
      double precision xvalue

c checks to see if the histogram has been initialized
      if (id.gt.nhist) then
        write(*,*) 'histogram not initialized'
        return
      endif

c if the value is outside of the range of the histogram it just exits.
      if ((x.gt.xmaxs(id)).or.(x.lt.xmins(id))) then
        return
      endif
  
c Binary search bitches!
      lowindex = 1
      highindex = nbins(id)+1

      do 100 while ((highindex-lowindex).gt.1)
        mid = (highindex+lowindex)/2
        xvalue = (xmaxs(id)-xmins(id))*dble(mid-1)/dble(nbins(id))+xmins(id)

        if (xvalue.le.x) lowindex = mid
        if (xvalue.gt.x) highindex = mid

 100  continue

c fill in the histogram
      hist(id,lowindex) = hist(id,lowindex) + zincrement

      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine hist2file(id,outscale)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This subroutine outputs the histogram data using the id to a file named 'label.dat'.           c
c  Outscale is an input scale that will scale all the bins by that value (i.e. scaling by         c
c  1/number of events will give you a differential cross section.                                 c
c                                                                                                 c
c  id = number to identify histogram                                                              c
c  outscale = scaling factor for the output                                                       c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

c input parameters
      integer id
      double precision outscale

c parameters for this subroutine only
      character*50 filename
      integer indexfile
      double precision xvalue

c checks to see if the histogram has been initialized
      if (id.gt.nhist) then
        write(*,*) 'histogram not initialized'
        return
      endif      

      indexfile = index(label(id),' ')-1
      filename = 'output/'//label(id)(1:indexfile)//'.txt'

c File for opening
      open(314,file=filename)
      write(314,*) '# Histogram label : ',label(id)(1:indexfile)
      do x1=1,nbins(id)
        xvalue = (xmaxs(id)-xmins(id))*dble(x1-1)/dble(nbins(id))+xmins(id)  
        write(314,fmt='(2E16.8)') xvalue, hist(id,x1)*outscale
      enddo
      write(314,fmt='(2E16.8)') xmaxs(id), 0.d0

      close(314)
      return
      end
