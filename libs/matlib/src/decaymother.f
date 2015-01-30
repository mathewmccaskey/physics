c-------------------------------------------------------------------------------------------------c
      subroutine decay_mother(decayunit,pdg_mother)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This subroutine will decay any mother particle in a MG/ME event record.  Assuming that         c
c  readeventMG has already been called the user supplies the fileunit for the decay event record  c
c  and the pdg code for the particle to be decayed.  The event record held in memory is updated   c
c  to include the new set of final state particles.                                               c
c                                                                                                 c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'
      
c input parameters
      integer decayunit, pdg_mother
      
c parameters for this subroutine only
      integer nparticles_old, pid_old(nmax), piif_old(nmax), mother1_old(nmax), mother2_old(nmax), 
     .  color1_old(nmax), color2_old(nmax)
      double precision pall_old(0:3,nmax), pmass_old(nmax), helicity_old(nmax)
      
      integer nparticles_new, pid_new(nmax), piif_new(nmax), mother1_new(nmax), mother2_new(nmax),
     .  color1_new(nmax), color2_new(nmax)
      double precision pall_new(0:3,nmax), pmass_new(nmax), helicity_new(nmax)
      
      double precision pframe(0:3), pboost(0:3)
      integer newparticles
      logical passed

c copy the old event record into corresponding old variables      
      nparticles_old = nparticles
        
      do x1=1,nparticles_old
        pid_old(x1) = pid(x1)
        piif_old(x1) = piif(x1)
        mother1_old(x1) = mother1(x1)
        mother2_old(x1) = mother2(x1)
        color1_old(x1) = color1(x1)
        color2_old(x1) = color2(x1)
        do x2=0,3
          pall_old(x2,x1) = pall(x2,x1)
        enddo
        pmass_old(x1) = pmass(x1)
        helicity_old(x1) = helicity(x1)
      enddo
        
      newparticles = 0

      do x1=1,nparticles_old
c looks for the desired mother particle in the event record
        if ((pid_old(x1).eq.pdg_mother).and.(piif_old(x1).eq.1)) then

c sets the particle to intermediate state
          piif_old(x1) = 2

c setting up the rest frame of the mother particle
          do x2=0,3
            pframe(x2) = pall_old(x2,x1)*gmunu(x2,x2)
          enddo

c reads in a decay event of the mother particle
          call readeventMG(decayunit,passed)

          if (passed) then
c excludes the mother particle in the decay event record
            do x2=2,nparticles
              newparticles = newparticles + 1
                
              pid_new(newparticles) = pid(x2)
              piif_new(newparticles) = piif(x2)
              if(mother1(x2).eq.1) then
                mother1_new(newparticles) = x1 + mother1(x2) - 1
              else
                mother1_new(newparticles) = nparticles_old + mother1(x2) - 1
              endif
              color1_new(newparticles) = color1(x2)
              color2_new(newparticles) = color2(x2)
c boosts all the final state particles in the rest frome of the mother
              do x3=0,3
                pboost(x3) = pall(x3,x2)
              enddo
                
              call boost(pframe,pboost)
                
              do x3=0,3
                pall_new(x3,newparticles) = pboost(x3)
              enddo
              pmass_new(newparticles) = pmass(x2)
              helicity_new(newparticles) = helicity(x2)
            enddo
          endif
        endif
      enddo

c total number of particles now in the event record
      nparticles = nparticles_old + nparticles_new
      
c resetting the old
      do x1=1,nparticles_old
        pid(x1) = pid_old(x1)
        piif(x1) = piif_old(x1)
        mother1(x1) = mother1_old(x1)
        mother2(x1) = mother2_old(x1)
        color1(x1) = color1_old(x1)
        color2(x1) = color2_old(x1)
        do x2=0,3
          pall(x2,x1) = pall_old(x2,x1)
        enddo
        pmass(x1) = pmass_old(x1)
        helicity(x1) = helicity_old(x1)
      enddo
      
      do x1=nparticles_old,nparticles
        pid(x1) = pid_new(x1)
        piif(x1) = piif_new(x1)
        mother1(x1) = mother1_new(x1)
        mother2(x1) = mother2_new(x1)
        color1(x1) = color1_new(x1)
        color2(x1) = color2_new(x1)
        do x2=0,3
          pall(x2,x1) = pall_new(x2,x1)
        enddo
        pmass(x1) = pmass_new(x1)
        helicity(x1) = helicity_new(x1)
      enddo         
      
      return
      end


  
c-------------------------------------------------------------------------------------------------c
      subroutine decay_mother_file(inputunit,decayunit,outputunit,pdg_mother,init)
c-------------------------------------------------------------------------------------------------c
c                                                                                                 c
c  This subroutine will decay any mother particle in a MG/ME event record.  User supplies         c
c  fileunits for the original event record, the decay event record, and the final event record.   c
c  The pdg code for the desired decay particle is provided.  Every instance of that mother        c
c  particle is decayed with the final state particles boosted in the rest frame of the mother     c
c  particle.  A new event record is created with the new final state particles.                   c
c                                                                                                 c
c init variable:                                                                                  c
c  -1 : Complete a full file in one run                                                           c
c   0 : Start a new output file in multiple runs                                                  c
c   1 : Continue writing to the output file                                                       c
c   2 : Complete output file                                                                      c
c-------------------------------------------------------------------------------------------------c
      implicit none
      include 'matlib.inc'

c input parameters
      integer inputunit, decayunit, outputunit, pdg_mother, init

c parameters for this subroutine only
      integer nparticles_old, indexevent_old
      double precision weight_old, energyscale_old, alphaqed_old, alphaqcd_old
      integer pid_old(nmax), piif_old(nmax), mother1_old(nmax), mother2_old(nmax), color1_old(nmax), color2_old(nmax)
      double precision pall_old(0:3,nmax), pmass_old(nmax), helicity_old(nmax)
      
      integer nparticles_new, indexevent_new
      double precision weight_new, energyscale_new, alphaqed_new, alphaqcd_new
      integer pid_new(nmax), piif_new(nmax), mother1_new(nmax), mother2_new(nmax), color1_new(nmax), color2_new(nmax)
      double precision pall_new(0:3,nmax), pmass_new(nmax), helicity_new(nmax)
      
      double precision pframe(0:3), pboost(0:3)
      integer i, j, k, newparticles, count
      logical passed, passed2


c Begin the output file with flag (for new file or full file only)
      if (init.le.0) then

c here is where one would write the information in the original lhe file to the new lhe file
        write(outputunit,fmt='(A18)') '<LesHouchesEvents>'
      endif
      
      passed = .true.
      passed2 = .true.
      count = 0

c here we would write the information in the decay file to the new lhe file      
      call readeventMG(inputunit,passed)
      
      do while(passed)
        count = count + 1
cc every 1000 events this will print out number of events processed
c        if (mod(count,1000).eq.0) write(*,*) 'Event #: ',count

c storing all the input particle information        
        nparticles_old = nparticles
        indexevent_old = indexevent
        weight_old = weight
        energyscale_old = energyscale
        alphaqed_old = alphaqed
        alphaqcd_old = alphaqcd
        
        do i=1,nparticles_old
          pid_old(i) = pid(i)
          piif_old(i) = piif(i)
          mother1_old(i) = mother1(i)
          mother2_old(i) = mother2(i)
          color1_old(i) = color1(i)
          color2_old(i) = color2(i)
          do j=0,3
            pall_old(j,i) = pall(j,i)
          enddo
          pmass_old(i) = pmass(i)
          helicity_old(i) = helicity(i)
        enddo
        
        newparticles = 0
        
        do i=1,nparticles_old
c looks for the desired mother particle in the event record
          if ((pid_old(i).eq.pdg_mother).and.(piif_old(i).eq.1)) then

c sets the particle to intermediate state
            piif_old(i) = 2

c setting up the rest frame of the mother particle
            do j=0,3
              pframe(j) = pall_old(j,i)*gmunu(j,j)
            enddo

c reads in a decay event of the mother particle
            call readeventMG(decayunit,passed2)

            if (passed2) then
c excludes the mother particle in the decay event record
              do j=2,nparticles
                newparticles = newparticles + 1
                
                pid_new(newparticles) = pid(j)
                piif_new(newparticles) = piif(j)
                if(mother1(j).eq.1) then
                  mother1_new(newparticles) = i + mother1(j) - 1
                else
                  mother1_new(newparticles) = nparticles_old + mother1(j) - 1
                endif
                color1_new(newparticles) = color1(j)
                color2_new(newparticles) = color2(j)

c boosts all the final state particles in the rest frome of the mother
                do k=0,3
                  pboost(k) = pall(k,j)
                enddo
                
                call boost(pframe,pboost)
                
                do k=0,3
                  pall_new(k,newparticles) = pboost(k)
                enddo
                pmass_new(newparticles) = pmass(j)
                helicity_new(newparticles) = helicity(j)
              enddo
            endif
          endif
        enddo

c writing out the full event record        
        write(outputunit,fmt='(A7)') '<event>'
        write(outputunit,fmt='(I4,I4,4E15.7)') nparticles_old+newparticles,indexevent_old,weight_old,
     .    energyscale_old,alphaqed_old,alphaqcd_old

c old particles first
        do i=1,nparticles_old
          write(outputunit,fmt='(I11,5I6,5E15.7,2F4.0)') pid_old(i),piif_old(i),mother1_old(i),mother2_old(i),color1_old(i),
     .      color2_old(i),pall_old(1,i),pall_old(2,i),pall_old(3,i),pall_old(0,i),pmass_old(i),0.0,helicity_old(i)
        enddo

c new final state particles
        do i=1,newparticles
          write(outputunit,fmt='(I11,5I6,5E15.7,2F4.0)') pid_new(i),piif_new(i),mother1_new(i),mother1_new(i),color1_new(i),
     .      color2_new(i),pall_new(1,i),pall_new(2,i),pall_new(3,i),pall_new(0,i),pmass_new(i),0.0,helicity_new(i)
        enddo
        write(outputunit,fmt='(A8)') '</event>'

c onto the next event      
        call readeventMG(inputunit,passed)  
      enddo
      
c finishes up the event record
      if ((init.eq.-1).or.(init.eq.2)) then
        write(outputunit,fmt='(A19)') '</LesHouchesEvents>'
        close(outputunit)
      endif
      
c close the input and decay files
      close(inputunit)
      close(decayunit)
      
      return
      end