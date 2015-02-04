c print out some status updates of the various calculations
      print_out = .false.

c flag to use the canonical calculation or the coupled relic density calculation
      relic_canonical = .true.

c flags for calculating an array of taacs values or calculate them on the spot
      calc_taacs_ann_array = .true.
      calc_taacs_dm2dm_array = .true.
      calc_taacs_scattering_array = .true.

c precision of the ODE integrator
      eps_ode = 1.0d-3

c precision of the romberg integration for Wij and taacs
      eps_wij = 1.0d-4
      eps_taacs = 1.0d-4

c minimum number of iterations in the romberg integration algorithm for both Wij and taacs
      iter_wij = 10
      iter_taacs = 10

c In the cacluation of the wij's these set the min and max of the adaptive step size
      beta_step_min = 0.0001d0
      beta_step_max = 0.01d0
