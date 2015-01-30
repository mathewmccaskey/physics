c function variable names
      double precision besselk, bessi0, bessi1, bessk0, bessk1, bessk, brent, cubeint,
     . dabscmplx, deltaphi, deltaR, Et, factorial, gmunu, invmass, lambda,
     . lineint, lineint2D, lineint3D, mt, mt2, p3dot, p4dot, pt, random,
     . randgauss, rapidity, theta3
      integer getindex

c generic input momenta for the functions
      double precision p(0:3), q(0:3)

c generic output momenta for decay subroutines
      double precision pf1(0:3), pf2(0:3), pf3(0:3), pf4(0:3)

c angles 
      double precision theta, phi

c wave functions
      double complex uin(4), uout(4), vin(4), vout(4), ubarout(1,4)

c counters
      integer x1, x2, x3, x4

c constants
      double precision pi, a_lep, a_jet, b_lep, b_jet
      parameter (pi=4.d0*datan(1.d0), a_lep=0.1d0, a_jet=0.5d0, b_lep=0.007d0, b_jet=0.03d0)

c gamma matrices
      double complex gamma0(0:3,0:3), gamma1(0:3,0:3), gamma2(0:3,0:3), gamma3(0:3,0:3)
c      common/gammas/ gamma0, gamma1, gamma2, gamma3

c parameters used for reading/manipulating MadGraph events
      integer nmax
      parameter (nmax=50)

      integer nevents, nparticles, indexevent, pid(nmax), piif(nmax), mother1(nmax), mother2(nmax), 
     .  color1(nmax), color2(nmax)
      common/readMGeventsint/ nevents, nparticles, indexevent, pid, piif, mother1, mother2, color1, color2 

      double precision cross, truncweight, unitweight, weight, energyscale, alphaqed, alphaqcd,
     .  pall(0:3,nmax), pmass(nmax), helicity(nmax)
      common/readMGeventsdble/ cross, truncweight, unitweight, weight, energyscale, alphaqed, alphaqcd, 
     .  pall, pmass, helicity

c parameters used for analyzing MadGraph events
      double precision ajet, bjet, alep, blep, pt_jet_cut, eta_jet_cut, pt_lep_cut, eta_lep_cut, dRjj, dRjl, dRll
      common/analyzecuts/ ajet, bjet, alep, blep, pt_jet_cut, eta_jet_cut, pt_lep_cut, eta_lep_cut, dRjj, dRjl, dRll

      double precision missing_pt(0:3), etas(nmax), pts(nmax), thetas(nmax), phis(nmax), pt_jet(nmax), pt_lep(nmax)
      common/analyzedble/ missing_pt, etas, pts, thetas, phis, pt_jet, pt_lep
      
      integer nlep, njet, nlep_pass, njet_pass, nfsp, ptype(nmax), fsp(nmax), fsp_jet(nmax), fsp_lep(nmax), 
     . pt_jet_order(nmax), pt_lep_order(nmax)
      common/analyzeint/ nlep, njet, nlep_pass, njet_pass, nfsp, ptype, fsp, fsp_jet, fsp_lep, 
     . pt_jet_order, pt_lep_order

      logical accepted(nmax)
      common/analyzebool/ accepted
      
c histogram parameters
      integer nhistmax, nbinmax
      parameter (nhistmax=100, nbinmax=1000)

      integer nhist, nbins(nhistmax)
      common/histint/ nhist, nbins

      character*30 label(nhistmax)
      common/histchar/ label

      double precision xmins(nhistmax), xmaxs(nhistmax), hist(nhistmax,nbinmax)
      common/histdble/ xmins, xmaxs, hist

c random seed
      integer idum, idum2
      common/randseed/ idum, idum2

c mt2 variables
      integer n_vis(2)
      common/mt2int/ n_vis
      
      double precision p_vis1(0:3), p_vis2(0:3), m_inv(2), pt_miss(0:3)
      common/mt2dble/ p_vis1, p_vis2, m_inv, pt_miss
