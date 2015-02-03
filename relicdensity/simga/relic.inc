c relativistic degrees of freedom
      double precision tempgstar(100), gstar(100), tempgstarS(100), gstarS(100)
      common/gstarstuff/ tempgstar, gstar, tempgstarS, gstarS

c taacs
      double precision xtaacs(100),taacs(100)
      common/taacsstuff/ xtaacs, taacs
      integer npoints
      common/moretaacsstuff/ npoints

c more variables for the relativistic degrees of freedom
      double precision g_star, g_star_S, temp

c passing x variable
      double precision x_pass
      common/passing/ x_pass

c model parameters
      double precision g, m
      common/modelparms/ g, m

c planck mass
      double precision mpl, tau
      parameter (mpl=1.2209d19, tau=8.d0*datan(1.d0))      
