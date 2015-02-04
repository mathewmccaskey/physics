c List of the process names in order of dmi, dmj, nprocesses(dmi, dmj)
c Annihilation process names
      process_names(1) = 'n1n1_uux'
      process_names(2) = 'n1n1_ddx'
      process_names(3) = 'n1n1_vevex'
      process_names(4) = 'n1n1_vtvtx'
      process_names(5) = 'n1n1_emep'
      process_names(6) = 'n1n1_ttx'
      process_names(7) = 'n1n1_bbx'
      process_names(8) = 'n1n1_zz'
      process_names(9) = 'n1n1_zh1'
      process_names(10) = 'n1n1_wpwm'
      process_names(11) = 'n1n1_h1h1'
      process_names(12) = 'n1n1_tamtap'

c DM -> DM process names

c DM/SM scattering process names

c Total number of processes for each category
      num_processes = 12
      ann_num_processes = 12
      dm2dm_num_processes = 0
      scattering_num_processes = 0

c Boolean operators for identical particles in the initial state
c Annihialtion diagrams
      ann_process_iden_init(1,1,1) = .true.
      ann_process_iden_init(1,1,2) = .true.
      ann_process_iden_init(1,1,3) = .true.
      ann_process_iden_init(1,1,4) = .true.
      ann_process_iden_init(1,1,5) = .true.
      ann_process_iden_init(1,1,6) = .true.
      ann_process_iden_init(1,1,7) = .true.
      ann_process_iden_init(1,1,8) = .true.
      ann_process_iden_init(1,1,9) = .true.
      ann_process_iden_init(1,1,10) = .true.
      ann_process_iden_init(1,1,11) = .true.
      ann_process_iden_init(1,1,12) = .true.

c DM -> DM diagrams

c Final state information for all the DM -> DM processes

c Initial state degrees of freedom for all the DM/SM scattering processes

c Total initial state degrees of freedom for all the DM/SM scattering processes
