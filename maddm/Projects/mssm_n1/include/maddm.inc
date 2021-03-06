c function names
c init.f
      double precision get_Y_eq, get_gstar, get_Wij_ann, Wij_ann_integrand
c interpolate.f
      integer getindex
      double precision get_Wij_value, get_taacs_value, lineint, bessk, bessk1
c phasespace.f
      double precision gmunu, invmass, lambda, p3dot, p4dot, random
c smatrix.f
      double precision smatrix_ann, smatrix_dm2dm, smatrix_scattering
c relicdensity_canon.f
      double precision relic_canon, taacs_canon, taacs_integrand_canon
c relic_coupled.f
      double precision relic_coupled, Wij_integrand_other, get_taacs_dm2dm, taacs_dm2dm_integrand
      double precision taacs_coupled, taacs_integrand
c tests.f
      double precision crossijk, cross_check_process


c random seed
      integer idum
      common/randint/ idum

c counting parameters
      integer x1, x2, x3, x4

c generic momentum 4-vectors
      double precision p(0:3), q(0:3)

c model parameters
      integer nparticles, isfermion(100)
      common/modelints/ nparticles, isfermion
      double precision masses(100), dofs(100)
      common/modeldbls/ masses, dofs

c relativistic degrees of freedom
      double precision gstar_x(401), gstar_boson(401), gstar_fermion(401)
      common/gstarstuff/ gstar_x, gstar_boson, gstar_fermion

c more variables for the relativistic degrees of freedom
      double precision g_star, g_star_S, temperature

c old gstar and gstarS stuff
      double precision gstar_temperature(100), gstar(100)
      double precision gstarS_temperature(100), gstarS(100)
      common/oldgstarstuff/ gstar_temperature, gstar, gstarS_temperature, gstarS

c variables needed for the phase space routines
      double precision pf1(0:3), pf2(0:3), pf3(0:3), pf4(0:3), pmass(4)

c dm paraeters
      integer maxdms, max_dm2dm
      parameter (maxdms=10, max_dm2dm=100)
      integer ndmparticles, nvar, dm_index(maxdms), dm_antiindex(maxdms)
      common/dmparms_int/ ndmparticles, nvar, dm_index, dm_antiindex
      integer num_xbar, xbar_index(maxdms), num_not_xbar, not_xbar_index(maxdms), rev_index(maxdms)
      common/xbar_int/ num_xbar, xbar_index, num_not_xbar, not_xbar_index, rev_index
      double precision mdm(maxdms), dof_dm(maxdms), dof_total(maxdms), dof_SM(maxdms,maxdms,max_dm2dm)
      double precision dof_SM_total(maxdms,maxdms,max_dm2dm), xf(maxdms), min_xbar_mass
      common/dmparms/ mdm, dof_dm, dof_total, dof_SM, dof_SM_total, xf, min_xbar_mass
      character(10) dm_names(maxdms), dm_antinames(maxdms)
      common/dmnames/ dm_names, dm_antinames
      logical solved_ODE(maxdms), found_xf(maxdms), dm_sm_scattering(maxdms), print_out, dmsm_warning
      common/dmparms_logical/ solved_ODE, found_xf, dm_sm_scattering, print_out, dmsm_warning

c relic density integrals
      integer nWij
      common/num_Wijs/ nWij
      double precision betas(1000), beta_step_min, beta_step_max
      common/beta_steps/ betas, beta_step_min,beta_step_max
      double precision Wij_ann(maxdms,maxdms,1000), Wij_dm2dm(maxdms,maxdms,max_dm2dm,1000)
      double precision Wij_scattering(maxdms,maxdms,max_dm2dm,1000), x_taacs(100), taacs_ann(maxdms,maxdms,100)
      double precision taacs_dm2dm(maxdms,maxdms,max_dm2dm,100), taacs_scattering(maxdms,maxdms,max_dm2dm,100)
      common/dmintegrals/ Wij_ann, Wij_dm2dm, Wij_scattering, x_taacs, taacs_ann, taacs_dm2dm, taacs_scattering

c number of annihilation particles and associated diagrams
      integer num_processes, ann_num_processes, dm2dm_num_processes, scattering_num_processes
      integer ann_nprocesses(maxdms,maxdms), dm2dm_nprocesses(maxdms,maxdms)
      integer scattering_nprocesses(maxdms,maxdms), dm2dm_fs(maxdms,maxdms,max_dm2dm,2)
      common/process_int/ num_processes, ann_num_processes, dm2dm_num_processes,
     . scattering_num_processes, ann_nprocesses, dm2dm_nprocesses, scattering_nprocesses, dm2dm_fs
      character(30) process_names(1000)
      common/process_names/ process_names
      logical ann_process_iden_init(maxdms,maxdms,100), dm2dm_process_iden_init(maxdms,maxdms,100)
      logical relic_canonical, calc_taacs_ann_array, calc_taacs_dm2dm_array, calc_taacs_scattering_array
      common/process_logical/ ann_process_iden_init, dm2dm_process_iden_init, relic_canonical,
     . calc_taacs_ann_array, calc_taacs_dm2dm_array, calc_taacs_scattering_array

c flags for the test subroutines to see if certain quantities have already been calculated
      logical Wij_ann_calc, Wij_dm2dm_calc, Wij_scattering_calc
      logical taacs_ann_calc, taacs_dm2dm_calc, taacs_scattering_calc, dm_sm_scattering_warn
      common/calc_checks/ Wij_ann_calc, Wij_dm2dm_calc, Wij_scattering_calc, 
     . taacs_ann_calc, taacs_dm2dm_calc, taacs_scattering_calc, dm_sm_scattering_warn

c passing variables
      double precision x_pass, beta_pass
      common/passing_dble/ x_pass, beta_pass
      integer dmi, dmj, p_k
      common/passing_int/ dmi, dmj, p_k

c precisions and number of iterations for intergration routines
      double precision eps_ode, eps_taacs, eps_wij
      common/epses/ eps_ode, eps_taacs, eps_wij
      integer iter_taacs, iter_wij
      common/iterations/ iter_taacs, iter_wij

c path common block for the ODE integrator
      integer kmax, kount
      double precision dxsav, xp(200), yp(50,200)
      common /path/ kmax, kount, dxsav, xp, yp

c constants
      double precision pi, mpl, ZERO
      parameter (pi=3.14159265359d0, mpl=1.2209d19, ZERO=0.d0)