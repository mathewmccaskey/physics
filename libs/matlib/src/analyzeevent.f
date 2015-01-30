c-------------------------------------------------------------------------------------------------c
      subroutine analyzeeventMG()
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This suroutine will analyze the event that is currently in memory.  It will calculate the pt,  c
c  eta, phi or all the final state particles.  The jets and leptons will be smeared, separated    c
c  and ordered with respect to pt.  pt, eta, and isolation cuts are all applied.  The final       c
c  number of accepted jets and leptons are tabulated and the total missing pt is calculated.      c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

c variables used in this subroutine only
      double precision pdummy1(0:3), pdummy2(0:3)
      
c reset the various particle counters
      nlep = 0
      njet = 0
      nlep_pass = 0
      njet_pass = 0
      nfsp = 0
      do x1=0,3
        missing_pt(x1) = 0.d0
      enddo
      
c filters out the initial and intermediate particles leaving only the final state        
      do x1=1,nparticles
        if (piif(x1).eq.1) then
          nfsp = nfsp + 1
          fsp(nfsp) = x1

c smears the final state particles according to alep, blep, ajet, and bjet            
          do x2=0,3
            pdummy1(x2) = pall(x2,x1)
          enddo

          if (abs(pid(x1)).le.5) then
            call smearp(pdummy1,ajet,bjet)
          else if ((abs(pid(x1)).eq.11).or.(abs(pid(x1)).eq.13).or.(abs(pid(x1)).eq.15)) then
            call smearp(pdummy1,alep,blep)
          endif

          do x2=0,3
            pall(x2,x1) = pdummy1(x2)
          enddo
          
c calculates the pt and rapidity of all the final state particles                    
          etas(nfsp) = rapidity(pdummy1)
          pts(nfsp) = pt(pdummy1)
          call getangles(pdummy1,thetas(nfsp),phis(nfsp))
          accepted(nfsp) = .true.

c separating the jets, leptons, and neutrinos            
          if (abs(pid(x1)).le.5) then
            ptype(nfsp) = 1
            njet = njet + 1
            fsp_jet(njet) = nfsp
            pt_jet(njet) = pts(nfsp)

c adding to the missing pt vector            
            do x2=1,2
              missing_pt(x2) = missing_pt(x2) - pdummy1(x2)
            enddo
            
          else if ((abs(pid(x1)).eq.11).or.(abs(pid(x1)).eq.13).or.(abs(pid(x1)).eq.15)) then
            ptype(nfsp) = 2
            nlep = nlep + 1
            fsp_lep(nlep) = nfsp
            pt_lep(nlep) = pts(nfsp)

c adding to the missing pt vector
            do x2=1,2
              missing_pt(x2) = missing_pt(x2) - pdummy1(x2)
            enddo
            
          else if ((abs(pid(x1)).eq.12).or.(abs(pid(x1)).eq.14).or.(abs(pid(x1)).eq.16)) then
            ptype(nfsp) = 3
          endif
        endif
      enddo

c pt sorting for the jets and leptons
      call eigsrt2(njet,pt_jet,pt_jet_order,2)
      call eigsrt2(nlep,pt_lep,pt_lep_order,2)

c acceptance cuts for the final state particles        
      do x1=1,nfsp

c acceptance cuts for jets
        if (ptype(x1).eq.1) then
          if ((pts(x1).le.pt_jet_cut).or.(dabs(etas(x1)).ge.eta_jet_cut)) then
            accepted(x1) = .false.
          endif
            
          do x2=0,3
            pdummy1(x2) = pall(x2,fsp(x1))
          enddo

c isolation cuts for jets            
          do x2=1,nfsp
            if (x1.ne.x2) then
              do x3=0,3
                pdummy2(x3) = pall(x3,fsp(x2))
              enddo  

c deltaR must be greater than the cut value (cut can be different between jet-jet and jet-lepton)
              if ((deltaR(pdummy1,pdummy2).le.dRjj).and.(ptype(x2).eq.1)) then
                accepted(x1) = .false.
              else if ((deltaR(pdummy1,pdummy2).le.dRjl).and.(ptype(x2).eq.2)) then
                accepted(x1) = .false.
              endif
            endif
          enddo

c acceptance cuts for leptons            
        else if (ptype(x1).eq.2) then
          if ((pts(x1).le.pt_lep_cut).or.(dabs(etas(1)).ge.eta_lep_cut)) then
            accepted(x1) = .false.
          endif
            
          do x2=0,3
            pdummy1(x2) = pall(x2,fsp(x1))
          enddo

c isolation cuts for leptons            
          do x2=1,nfsp
            if (x1.ne.x2) then
              do x3=0,3
                pdummy2(x3) = pall(x3,fsp(x2))
              enddo  
                
c deltaR must be greater than the cut value (cut can be different between lepton-lepton and jet-lepton)
              if ((ptype(x2).eq.1).and.(deltaR(pdummy1,pdummy2).le.dRjl)) then
                accepted(x1) = .false.
              else if ((ptype(x2).eq.2).and.(deltaR(pdummy1,pdummy2).le.dRll)) then
                accepted(x1) = .false.
              endif    
            endif
          enddo
        endif
      enddo

      do x1=1,nfsp
        if (accepted(x1)) then
          if (ptype(x1).eq.1) then
            njet_pass = njet_pass + 1
          else if (ptype(x1).eq.2) then
            nlep_pass = nlep_pass + 1
          endif
        endif
      enddo
      
      return
      end
      
      
      
c-------------------------------------------------------------------------------------------------c
      subroutine reseteventMG()
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This subroutine will reset the entire event record that is currently stored in memory.  This   c
c  will be a good way to make sure that nothing from an old event can mess up any new event       c
c  read.                                                                                          c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'
      
      nevents = 0
      nparticles = 0
      indexevent = 0
      cross = 0.d0
      truncweight = 0.d0
      unitweight = 0.d0
      weight = 0.d0
      energyscale = 0.d0
      alphaqed = 0.d0
      alphaqcd = 0.d0
            
      nlep = 0
      njet = 0
      nlep_pass = 0
      njet_pass = 0
      nfsp = 0
      do x1=0,3
        missing_pt(x1) = 0.d0
      enddo
      
      do x1=1,50
        do x2=0,3
          pall(x2,x1) = 0.d0
        enddo
        pmass(x1) = 0.d0
        helicity(x1) = 0.d0

        pid(x1) = 0
        piif(x1) = 0
        mother1(x1) = 0
        mother2(x1) = 0 
        color1(x1) = 0
        color2(x1) = 0
                
        etas(x1) = 0.d0
        pts(x1) = 0.d0
        thetas(x1) = 0.d0
        phis(x1) = 0.d0
        pt_jet(x1) = 0.d0
        pt_lep(x1) = 0.d0
        ptype(x1) = 0
        fsp(x1) = 0
        fsp_jet(x1) = 0
        fsp_lep(x1) = 0 
        pt_jet_order(x1) = 0
        pt_lep_order(x1) = 0
        accepted(x1) = .true.
      enddo
      
      return
      end