c-------------------------------------------------------------------------c
      program maddm
c-------------------------------------------------------------------------c
c                                                                         c
c  This is the main driver for maddm.  It can be used to call all the     c
c  desired subroutines and functions.                                     c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'
      include 'input.inc'
      include 'coupl.inc'

c parameters used in this routine only
      double precision Oh2, Oh2_old

c include files generated my the python side that contains necessary information
      include 'maddm_card.inc'
      include 'diagrams.inc'

c Sets all the parameters used in MG/ME to the values in the param card
      call setpara('Cards/param_card.dat')

c If the user wants to change any of the parameters from within the fortran code then
c just change the variables given in 'input.inc' and simply run the subroutine coup()
c
c NOTE: If there are any other model parameters that are dependent on the parameters
c being changed then make sure those parameters are also changed as well
c 
c Example:
c      del1 = 0.11d0
c      Mx1 = 500.d0
c      call coup()


c subroutine to initialize the necessary parameters for the calculation
      call init_relic()

c calculate the Wijs for all the annihilation channels
      call calculate_Wij_ann()

c calculate the relic abundance of the DM candidate using either the canonical or coupled method
      if (relic_canonical) then
        nvar = 1
        Oh2 = relic_canon()
	write(*,*) Oh2
C c In order to call odeint_check make sure dxsav, kmax, and kount are set in relicdensity.f
C         call odeint_check(0)

      else
        nvar = ndmparticles
        Oh2 =  relic_coupled()
	write(*,*) Oh2

C c In order to call odeint_check make sure dxsav, kmax, and kout are set in relic_old.f
C         call odeint_check(1)

      endif
      
C       write(*,*) 'The relic density calculated in maddm is ',Oh2
      
c-------------------------------------------------------------------------c
c  Other test functions
c-------------------------------------------------------------------------c

C c bessel function check (min_x, max_x, nsteps, logscale)
C       call bessel_check(0.1d0,100.d0,201,1)
C 

C c relativisitic degrees of freedom checks (min_temp, max_temp, nsteps, logscale)
C       call dof_check(0.001d0,100.d0,201,1)
C 

C c Wij check (prints out all array of Wij's calculated in relicdensity())
C       call Wij_check()
C 

C c taacs test (min_x, max_x, nsteps, logscale)
C       call taacs_check(1.00d0,100d0,101,1)
C 

C c check all cross sections for a range of inital state momenta (xinit_min, x_init_max, nsteps, p_or_E, logscale)
C       call cross_check_scan_all(1.d0,1.d3,201,1,1)
C 

C c cross section for an individual process over range of inital state momenta
C c (dm_i, dm_j, process_k, xinit_min, xinit_max, nsteps, p_or_E, logscale)
C       call cross_check_scan(1,1,1,1.d0,1.0d3,201,1,1)
C 

C c cross-section for all processes at a single initial state momentum or energy (x1init, x2init, p_or_E)
C       call cross_check_all_process(1.d3,1.d3,1)
C 

C c cross-section check for a single process at a single initial state momentum or energy
C c (dm_i, dm_j, process_k, x1init, x2init, p_or_E)
C       write(*,*) cross_check_process(1,1,1,1.d3,1.d3,1)
C

C c matrix element check for all processes at a single initial state momentum or energy (x1init, x2init, p_or_E)
C       call matrix_element_check_all(1.0d3,1.0d3,0)
C

C c matrix element check for a single process at a single initial state momentum or energy
C c (dm_i, dm_j, process_k, x1init, x2init, p_or_E)
C       call matrix_element_check_process(1,1,1,1.0d3,1.0d3,0)

      end
