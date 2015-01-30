c-------------------------------------------------------------------------------------------------c
      subroutine readeventMG(fileunit,pass)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This subroutine reads in an event that is output from MadGraph.  The event folder should       c
c  already be open before calling this subroutine.                                                c
c                                                                                                 c
c  fileunit - fileunit of the MG .lhe file                                                        c
c  nparticles - number of particles                                                               c
c  indexevent - index of event (not very important)                                               c
c  weight - weight of the evernt                                                                  c
c  energyscale - energy scale of the event                                                        c
c  alphaqed - alpha QED at the energy scale escale                                                c
c  alphaqcd - alpha QCD at the energy scale escale                                                c
c  pid - particle ID (PDG code)                                                                   c
c  piif- initial/intermediate/final (-1,2,1) state particle                                       c
c  mother1/2 - mother particle(s) listed by particle number                                       c
c  color1/2 - color/anti-color of particle                                                        c
c  pall - 4-momenta of the particles in the event                                                 c
c  mass - it's the mass.... duh                                                                   c
c  helicity - it's the helicity of the particle                                                   c
c                                                                                                 c
c  IMPORTANT NOTE: Make sure to initialize the pass variable to true before sending it to         c
c  readevenMG.                                                                                    c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'
      
c parameters used in this subroutine only
      character*40 flag
      double precision junk

c input parameters
      integer fileunit
      logical pass

      flag = 'initialize'

      do while ((pass).and.(flag.ne.'<event>'))
        read(unit=fileunit,fmt='(A20)') flag
        if (flag.eq.'</LesHouchesEvents>') then
          pass = .false.
          return
        endif
      enddo

      read(fileunit,*) nparticles, indexevent, weight, energyscale, alphaqed, alphaqcd

      do x1=1, nparticles
        read(fileunit,*) pid(x1), piif(x1), mother1(x1), mother2(x1), color1(x1), 
     .    color2(x1), pall(1,x1), pall(2,x1), pall(3,x1), pall(0,x1), pmass(x1), junk, helicity(x1)
      enddo

      flag = 'reset'

      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine readgeninfoMG(fileunit)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This subroutine reads the generation information that MadGraph uses to generate events.        c
c                                                                                                 c
c  fileunit = user supplied file unit of MG .lhe file                                             c
c  numevents = number of generated events                                                         c
c  intweight = integrated weight                                                                  c
c  truncweight = truncated weight                                                                 c
c  unitweight = intweight/numevents                                                               c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

c input parameters
      integer fileunit

c parameters for this subroutine only
      character*20 flag, junk

      flag = 'initialize'

      do while (flag.ne.'<MGGenerationInfo>')
c        write(*,*) flag
        read(fileunit,fmt='(A20)') flag
      enddo

      read(fileunit,*) junk,junk,junk,junk,junk,nevents
      read(fileunit,*) junk,junk,junk,junk,junk,cross
      read(fileunit,*) junk,junk,junk,junk,junk,truncweight
      read(fileunit,*) junk,junk,junk,junk,unitweight

      return
      end



c-------------------------------------------------------------------------------------------------c
      subroutine writeeventMG()
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This subroutine simply outputs a MG event that is read using the subroutine readMGevent.       c
c  This will be helpful in debugging purposes.                                                    c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

      write(*,*) '---------------------------------------------------------'
      write(*,fmt='(2(A12),4(A16))') '# particles','index','weight','E scale','alpha qed','alpha qcd'
      
c writes out the MG event info      
      write(*,fmt='(2I12,4F16.8)') nparticles, indexevent, weight, energyscale, alphaqed, alphaqcd

      write(*,fmt='(6(A12),6(A16))') 'PDG code','i/i/f','mother1','mother2','color1','color2','px','py',
     .  'pz','E','mass','helicity'

c writes out the particle information      
      do x1=1,nparticles
        write(*,fmt='(6I12,6F16.8)') pid(x1), piif(x1), mother1(x1), mother2(x1),
     .    color1(x1), color2(x1), pall(1,x1), pall(2,x1), pall(3,x1), pall(0,x1),
     .    pmass(x1), helicity(x1)
      enddo
      
      write(*,*) '---------------------------------------------------------'
      
      return
      end
      
      
      
c-------------------------------------------------------------------------------------------------c
      subroutine writeeventMGfull()
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This subroutine not only writes out the full event record via the subroutine writeeventMG it   c
c  will output all the information calculated in analyzeevent.                                    c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

      call writeeventMG()

      write(*,fmt='(A20,I3,A30,I3)') 'Number of jets: ',njet,'Number of isolated jets :',njet_pass
      write(*,*) 'pT ordered jets: '
      write(*,fmt='(2A7,4A16,A12)') 'index','pid','pT','eta','phi','theta','accepted?'
      do x1=1,njet
        write(*,fmt='(2I7,4F16.10,L12)') fsp(fsp_jet(pt_jet_order(x1))), pid(fsp(fsp_jet(pt_jet_order(x1)))), 
     .      pt_jet(pt_jet_order(x1)), etas(fsp_jet(pt_jet_order(x1))), phis(fsp_jet(pt_jet_order(x1))), 
     .      thetas(fsp_jet(pt_jet_order(x1))), accepted(fsp_jet(pt_jet_order(x1)))
      enddo

      write(*,*) '--------------------------------------------------------------------------'
      write(*,fmt='(A20,I3,A30,I3)') 'Number of leptons: ',nlep,'Number of isolated leptons :',nlep_pass          
      write(*,*) 'pT ordered leptons: '
      write(*,fmt='(2A7,4A16,A12)') 'index','pid','pT','eta','phi','theta','accepted?'
      do x1=1,nlep
        write(*,fmt='(2I7,4F16.10,L12)') fsp(fsp_lep(pt_lep_order(x1))), pid(fsp(fsp_lep(pt_lep_order(x1)))),
     .      pt_lep(pt_lep_order(x1)), etas(fsp_lep(pt_lep_order(x1))), phis(fsp_lep(pt_lep_order(x1))),
     .      thetas(fsp_lep(pt_lep_order(x1))), accepted(fsp_lep(pt_lep_order(x1)))
      enddo      

      return
      end