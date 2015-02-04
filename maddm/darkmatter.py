import os
import shutil
import sys
import subprocess
import fileinput
import glob
import pickle

# python files from the madgraph source directory
import madgraph.core.base_objects as base_objects
import madgraph.interface.madgraph_interface as madgraph_interface
import models.model_reader as model_reader
import models.check_param_card as check_param_card
import madgraph.iolibs.export_v4 as export_v4
import madgraph.core.helas_objects as helas_objects
import madgraph.iolibs.file_writers as writers
import madgraph.various.misc as misc

# Output routines for MadDM
import MGoutput

# Root path
rp = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]
sys.path.append(rp)

# Current directory path
pp = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[1]

sys.path.append(rp+'/'+pp)





#-------------------------------------------------------------------------#
class darkmatter(base_objects.Particle):
#-------------------------------------------------------------------------#
#                                                                         #
#  This class conains all the routines and functions that are needed to   #
#  do the following actions.                                              #
#                                                                         #
#  1. Determine the Dm candidate and coannihilation particles of a        #
#     given model.                                                        #
#  2. Generate all the relevant annihilation diagrams that are needed to  #
#     calculate the overall relic abundance of DM.                        #
#  3. Create a fortran project where all the numerical calculations will  #
#     take place.                                                         #
#                                                                         #
#-------------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def __init__(self):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  This routine initializes some of the lists that are used to find     #
    #  the dark matter candidate and coannihilation partners needed for     #
    #  the relic abundance calculation.                                     #
    #                                                                       #
    #-----------------------------------------------------------------------#

      # Initializes useful variables
      self._bsm_particles = []
      self._bsm_masses = []
      self._bsm_final_states = []
      self._dm_particles = []
      self._coann_particles = []
      self._dm_names_list = []
      self._dm_antinames_list = []
      self._dm_thermal_scattering = []
      self._wanted_lorentz = []
      self._wanted_couplings = []
      self._new_proj = True
      self._modelname = ''
      self._projectname = ''
      self._paramcard = ''

            
    #-----------------------------------------------------------------------#
    def init_from_model(self,modelname, projectname, new_proj = True):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  Given a user input model this routine initializes the project        #
    #  folder and initializes several class variables used in other         #
    #  routines.                                                            #  
    #                                                                       #
    #-----------------------------------------------------------------------#

      # Set up the model and project name
      self._modelname = modelname
      self._projectname = projectname
      self._new_proj = new_proj

      #print "------ Importing Model ------"

      # MadGraph 5 command interface 
      self._MG5Cmd = madgraph_interface.MadGraphCmd()
      self._mgme_dir = self._MG5Cmd._mgme_dir
      self._MG5Cmd.do_import('model %s' % self._modelname)
      self._Model = self._MG5Cmd._curr_model

      # Gets a list of all the particles in the user-supplied model
      self._particles = self._Model.get('particles')

      #print "----- Creating Default Param Card ------"

      # Copy over the param_card.dat from v4 model files. If the project already exists
      # do nothing but set the member variables to point to the project files.
      if self._new_proj:
      	if self._MG5Cmd._model_v4_path:
        	shutil.copy2(self._mgme_dir+'/Models/'+modelname+'/param_card.dat','Projects/param_card.dat')
     	 # If the model is v5, generate the default param card from the python files
      	else:
        	write_dir = 'Projects/'
        	model_builder = export_v4.UFO_model_to_mg4(self._Model, write_dir)
        	model_builder.create_param_card()

      	# Initialize the paramcard variable
      	self._paramcard = 'Projects/param_card.dat'
      	self._project_exists = False

      	if ((self._modelname == 'mssm') or (self._modelname[0:5] == 'mssm-')):
        	print "Warning: You are using a restricted model. MadDM will automatically change the parameter card format."
        	check_param_card.convert_to_mg5card(self._paramcard, self._paramcard)

      else:
	  	project_list = glob.glob('Projects/'+self._projectname+'_*')
	  	self._projectpath = project_list[0]
	  	self._paramcard = self._projectpath+'/Cards/param_card.dat'

      #print "------ Initialization Complete! ------\n"

    #-----------------------------------------------------------------------#
    def FindDMCandidate(self, prompts = True):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  This routine finds the dark matter candidate and assigns it to the   #
    #  self._dm_particles list.  There are two ways in which this is done.  #
    #  Either the user can input the DM candidate or let the following      #
    #  algorithm find the dm candidates.                                    #
    #                                                                       #
    #  1. The DM particle must be a BSM particle (pdg_code > 25)            #
    #  2. The particle should have no charge (electric or color)            #
    #  3. The particle's width should be 0 or 'ZERO'                        #
    #  4. The assigned DM candidate is the lightest of all the particles    #
    #     that meet the above criteria.                                     #
    #                                                                       #
    #-----------------------------------------------------------------------#

      # initialize boolean flags for the multiple ways of inputting the DM particles
      found_DM_particle = False
      self._ask_param_card = False

      # All the particles with pdg_codes > 25 are added to the array of BSM particles
      for i in range(len(self._particles)):
        if (self._particles[i]['pdg_code'] > 25):
          self._bsm_particles.append(i)

      # When manually entering in the DM and coannihilation particles, both particle and anti-particle
      # are not both required to be listed.  If one is not included then it is added later.
      while (not found_DM_particle):
        
        if (prompts == True):
      
          print "Enter DM candidate (press Enter to automatically find the DM candidate): "
          print " (Enter 'particles' to see a list of particles)"
          dm_answer = raw_input()
        else:
          dm_answer = ''

        # Here is the algorithm for automatically finding the DM particle and coannihilation candidates
        if (dm_answer == ''):

          # Ask if the param card needs to be changed
          if (not self._ask_param_card):
            self.ChangeParamCard(prompts)
            self._ask_param_card = True

          # Looping over all the BSM particles we check the criteria described above
          for i in range(len(self._bsm_particles)):
            if (self._particles[self._bsm_particles[i]]['charge'] == 0.0):

              # We separate the particles that have 'ZERO' for the width since we don't need to numerically evaluate it
              if (self._particles[self._bsm_particles[i]]['width'] == 'ZERO'):

                # If nothing has been found so far then set the first DM candidate
                if (len(self._dm_particles) == 0):
                  self._dm_particles.append(self._particles[self._bsm_particles[i]])
                  self._dm_mass = self._bsm_masses[i]
                # If we already found a candidate, comare the masses and keep the one that's lighter
                elif (self._bsm_masses[i] < self._dm_mass):
                  self._dm_particles[0] = self._particles[self._bsm_particles[i]]
                  self._dm_mass = self._bsm_masses[i]

              # If the width is not set to 'ZERO' then we have to get the width from the param card
              elif (self.GetWidth(self._particles[self._bsm_particles[i]]['pdg_code']) == 0.0):
                if (len(self._dm_particles) == 0):
                  self._dm_particles.append(self._particles[self._bsm_particles[i]])
                  self._dm_mass = self._bsm_masses[i]
                elif (self._bsm_masses[i] < dm_mass):
                  self._dm_particles[0] = self._particles[self._bsm_particles[i]]
                  self._dm_mass = self._bsm_masses[i]

          # Check to see if we actually found a DM candidate
          if (self._dm_particles == []):
            print "No dark matter candidates in the model!"
            sys.exit(0) 
          else:
            found_DM_particle = True

        elif (dm_answer != 'particles'):
          # We loop over all the particles and check to see if the desired DM candidate is indeed in the model
          for particle in self._particles:
            if ((particle['name'] == dm_answer) or (particle['antiname'] == dm_answer)):
              self._dm_particles.append(particle)
              found_DM_particle = True
              self._dm_mass = -1.0

          # Check to see if we found the desired DM candidate in the model.
          if (self._dm_particles == []):
            print "Dark Matter candidate not present in the model! Try again."

        else:
          print ''
          self._MG5Cmd.do_display('particles')
          print ''

      # Print out the DM candidate
      if prompts == True:
     	print "-------------------------------"
      	print "DARK MATTER CANDIDATE:"
      	print "-------------------------------"
      	print self._dm_particles[0]



    #-----------------------------------------------------------------------#
    def FindCoannParticles(self, prompts = True, coann_eps = 0.1):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  This routine finds the coannihilation particles for the relic        #
    #  density calculation.  Either the user can manually input the desired #
    #  particles, or the code can search for all the BSM particles that are #
    #  within an input mass difference with the DM candidate.  All          #
    #  coannihilation particles are then added to the self._dm_particles    #
    #  list.                                                                #
    #                                                                       #
    #-----------------------------------------------------------------------#

      # initialize boolean flags for the multiple ways of inputting the Coannihilation particles
      found_DM_particle = False
      legit_coann_answer = False

      # Time to find the coannihilation particles
      while (not legit_coann_answer):

        if prompts == True:
            print "Enter the coannihilation particles:"
            print "(press Enter to automatically find the coannihilation particles)"
            print "(Enter 'particles' to see a list of particles)"
            coann_answer = raw_input()

        else:
          coann_answer =''

        # Automatically find the coannihilation candidates
        if (coann_answer == ''):
          legit_coann_answer = True

          ratio_answer = ''
          if prompts == True:
             print "Enter the mass difference ratio desired for coannihilating particles [0.1]:"
             print "(Enter 0 for no coannihilating particles)"
             ratio_answer = raw_input()
        
          if ratio_answer == '':
            self._coann_eps = coann_eps
          else:
            self._coann_eps = float(ratio_answer)

          # If the param card hasn't been changed then as if the param card needs to be changed.
          if (not self._ask_param_card):
            self.ChangeParamCard(prompts)
            self._ask_param_card = True
            if (self._dm_mass < 0.0):
              self._dm_mass = abs(self.GetMass(self._dm_particles[0]['pdg_code']))

          # If the user wishes to find coannihilation candidates we simply loop over the rest of the BSM particles
          # and see which particles have a mass within the input fractional mass difference.
          if (self._coann_eps > 0.0):

            # Loop over BSM particles
            for i in range(len(self._bsm_particles)):
              if (self._particles[self._bsm_particles[i]] != self._dm_particles[0]):
                if (abs(self._dm_mass-self._bsm_masses[i])/self._dm_mass <= self._coann_eps):
                  self._coann_particles.append(self._particles[self._bsm_particles[i]])

                # If there are BSM particles that are too small to be included in the coannihilation they
                # are still tabulated to include in the final state particles with the SM particles.
                elif (self._bsm_masses[i] < (1.0-self._coann_eps)*self._dm_mass):
                  self._bsm_final_states.append(self._particles[self._bsm_particles[i]])


        # This is the case where the user inputs their own set of coannihilation particles
        elif (coann_answer != 'particles'):
          legit_coann_answer = True

          # break up the string into a list of particles
          input_coann_names_list = coann_answer.split(" ")

          # Loop over the model particles so we can add them to the self._dm_particles list
          for i in range(len(input_coann_names_list)):
            for particle in self._particles:
              # Checks to see if either the particle or anti-particle is in the model as well as if the
              # particle isn't already included in the list of DM particles (to avoid double counting)
              if ((particle['name'] == input_coann_names_list[i]) and (not (particle in self._dm_particles)) \
                  and (not (particle in self._coann_particles))):
                self._coann_particles.append(particle)
              elif ((particle['antiname'] == input_coann_names_list[i]) and (not (particle in self._dm_particles)) \
                  and (not (particle in self._coann_particles))):
                self._coann_particles.append(particle)

          # If the param card hasn't been changed then as if the param card needs to be changed.
          if (not self._ask_param_card):
            self.ChangeParamCard(prompts)
            self._ask_param_card = True
            if (self._dm_mass < 0.0):
              self._dm_mass = abs(self.GetMass(self._dm_particles[0]['pdg_code']))

          # Tabulates all the BSM particles that are lighter than the DM candidate so they can be included in the
          # final state particles along with the SM particles.
          for i in range(len(self._bsm_particles)):
            if ((not (self._particles[self._bsm_particles[i]] in self._dm_particles)) and \
                (not (self._particles[self._bsm_particles[i]] in self._coann_particles))):
              if (self._bsm_masses[i] < self._dm_mass):
                self._bsm_final_states.append(self._particles[self._bsm_particles[i]])

        else:
          print ''
          self._MG5Cmd.do_display('particles')
          print ''

      # For organizational purposes we put the coannihilation particles by alphabetical order by name
      coann_ordered = sorted(self._coann_particles, key=lambda k: k['name']) 
      self._dm_particles += coann_ordered

      # If we found any coannihilation particles then print out the candidates
      if (len(self._coann_particles) > 0 and prompts == True):
        print "-------------------------------"
        print "COANNIHILATION PARTNERS:"
        print "-------------------------------"
        for i in range(len(self._coann_particles)):
          print self._coann_particles[i]



    #-----------------------------------------------------------------------#
    def GetProjectName(self):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  Once we find the DM and Coannihilation particles we can append the   #
    #  the project name and see if the project has already been generated   #
    #                                                                       #
    #-----------------------------------------------------------------------#

      project_suffix = ''

      # For all the DM particles we create lists that contain all the names and append to the project suffix
      for dm_particle in self._dm_particles:
        self._dm_names_list.append(dm_particle['name'])
        self._dm_antinames_list.append(dm_particle['antiname'])
        project_suffix += dm_particle['name']

      # Convert the '+', '-', and '~' to 'p', 'm', and 'x' respectively
      project_suffix = self.Convertname(project_suffix)

      # Now we can set the project name and project path
      self._projectname += '_'+project_suffix
      self._projectpath = os.path.join('Projects', self._projectname)

      if os.path.isdir(self._projectpath):
        self._project_exists = True



    #-----------------------------------------------------------------------#
    def ChangeParamCard(self, prompts = True):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  This routine allows the user to either edit the default param card   #
    #  or enter in the location of the param card that they wish to use.    #
    #  if the param card entered doesn't exist we'll just repeat until a    #
    #  valid paramcard or answer is given.                                  #
    #                                                                       #
    #-----------------------------------------------------------------------#

      # Loop until we have a vaid result to change (or not change) the param card
      change_param_card = False
      while (not change_param_card):

        # Ask if the param_card needs to be changed
        answer = ''
        if prompts == True:
            print "Would you like to edit the default param_card?[n] (y/n):"
            print "(or enter the location of param_card to be used)"
            answer = raw_input()

        # Edit the default param card
        if ((answer == 'y') or (answer == 'Y')):
          subprocess.call('vi %s' % self._paramcard, shell= True)
          change_param_card = True
        # Do nothing
        elif ((answer == 'n') or (answer == 'N') or (answer == '')):
          change_param_card = True
        # Check to see if the entered file exists.  If so, copy it to the param card
        else:
          if (os.path.isfile(answer)):
            shutil.copy2(answer,self._paramcard)
            # If the model is based off of the mssm, then we may need to convert the param card
            # from the slha1 format to the slha2 format
            if ((self._modelname == 'mssm') or (self._modelname[0:5] == 'mssm-')):
              check_param_card.convert_to_mg5card(self._paramcard, self._paramcard)
            change_param_card = True
          else:
            print "\nThe file entered does not exist. Try Again.\n"

      # Now that we have the param card set, we can set up the model reader
      # This is necessary to get the numerical value for the particle masses
      self._fullmodel = model_reader.ModelReader(self._Model)
      self._fullmodel.set_parameters_and_couplings(param_card=self._paramcard)

      # For all the bsm particles we create the array of each particle's mass.
      for i in self._bsm_particles:
        self._bsm_masses.append(abs(self.GetMass(self._particles[i]['pdg_code'])))



    #-----------------------------------------------------------------------#
    def GetMass(self, pdg_id):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  Finds the mass of a particle from a Model object given a PDG code    #
    #  and returns it. CAUTION: Masses in MadGraph are stored as complex    #
    #  numbers.
    #                                                                       #
    #-----------------------------------------------------------------------#
    
      mass = self._fullmodel.get('parameter_dict')[self._Model.get_particle(pdg_id).get('mass')].real
                    
      return mass



    #-----------------------------------------------------------------------#
    def GetWidth(self, pdg_id):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  Same as the previous method but for the width.                       #
    #                                                                       #
    #-----------------------------------------------------------------------#

      width = self._fullmodel.get('parameter_dict')[self._Model.get_particle(pdg_id).get('width')].real
                    
      return width



    #-----------------------------------------------------------------------#
    def Convertname(self, name):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  This routine takes a character sting, and for every '+', '-', and    #
    #  '~' in the name it is replaced with a 'p', 'm', and 'x'              #
    #  respectively.  This is used for both project and process names.      #
    #                                                                       #
    #-----------------------------------------------------------------------#

      # Gets a list of the individual characters
      name_list = list(name)
      # Loops over the characters and makes appropriate changes
      for i in range(len(name_list)):
        if (name_list[i] == '+'):
          name_list[i] = 'p'
        elif (name_list[i] == '-'):
          name_list[i] = 'm'
        elif (name_list[i] == '~'):
          name_list[i] = 'x'
      # Combines everything to the new name
      name = ''.join(name_list)
      return name



    #-----------------------------------------------------------------------#
    def GenerateDiagrams(self):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  This routine generates all the 2 -> 2 annihilation matrix elements.  #
    #  The SM particles, BSM particles in self._bsm_final_states, and the   #
    #  DM particles are used as final states.  The initial states are       #
    #  looped over all the different combinations of DM particles.          #
    #                                                                       #
    #-----------------------------------------------------------------------#

      # If the project already exists we don't need to generate diagrams
      if (self._project_exists):
        return

      # Else generate the matrix elements
      #print "----- Generating Matrix Elements -----\n",

      # Initialize the arrays that will store all the matrix element information
      self._annihilation_me = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
      self._dm2dm_me = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
      self._scattering_me = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]

      # Set up the initial state multiparticles that contain the particle and antiparticle
      self._DM_all_names = list(set(self._dm_names_list + self._dm_antinames_list))
      self._DM_init_state = []
      for i in range(len(self._dm_particles)):
        if (self._dm_names_list[i] != self._dm_antinames_list[i]):
          self._DM_init_state.append(self._dm_names_list[i]+' '+self._dm_antinames_list[i])
        else:
          self._DM_init_state.append(self._dm_names_list[i])
        self._MG5Cmd.do_define('DM_particle'+str(i+1)+' = '+self._DM_init_state[i])

      # Set up the SM pdg codes
      leptons = range(11, 17)
      quarks = range(1, 7)
      bosons = range(21, 26)
      sm_pdgs = leptons + quarks + bosons

      # Get the names of the SM particles from the pdg codes
      sm_names_list = []
      for particle in self._particles:
        if particle['pdg_code'] in sm_pdgs:
          sm_names_list.append(particle['name'])
          if (not particle['self_antipart']):
            sm_names_list.append(particle['antiname'])

      # With the list of BSM particles tabulated in FindDMCandidate we can get the list of bsm names
      bsm_names_list = []
      for bsm_particle in self._bsm_final_states:        
        bsm_names_list.append(bsm_particle['name'])
        if (not bsm_particle['self_antipart']):
          bsm_names_list.append(bsm_particle['antiname'])

      # Create a MadGraph multiparticle that conatins all three groups of particles
      sm_names = ' '.join(sm_names_list)
      bsm_names = ' '.join(bsm_names_list)
      dm_names = ' '.join(self._DM_all_names)
      self._MG5Cmd.do_define('dm_particles = '+dm_names)
      self._MG5Cmd.do_define('fs_particles = '+sm_names+' '+bsm_names)


      # Generate the annihilation diagrams by going through all the combinations of
      # initial state particles (including coannihilations)
      for i in range(len(self._dm_particles)):
        for j in range(i,len(self._dm_particles)):

          # Create the appropriate string of initial to final state particles
          try:
            proc = 'DM_particle'+str(i+1)+' DM_particle'+str(j+1)+' > fs_particles fs_particles'
            self._MG5Cmd.do_generate(proc)

            # Once we generate the diagrams, we then immediately get the matrix elements for this process
            # so we can keep track of the number of processes as well as generate the next set
            curr_matrix_elements = helas_objects.HelasMultiProcess(self._MG5Cmd._curr_amps)
            self._annihilation_me[i][j] = curr_matrix_elements.get_matrix_elements()
            self._wanted_lorentz += curr_matrix_elements.get_used_lorentz()
            self._wanted_couplings += curr_matrix_elements.get_used_couplings()
          except:
            continue

      # Generate all the DM -> DM processes
      for i in range(len(self._dm_particles)):
        for j in range(i,len(self._dm_particles)):

          # Create the appropriate string of initial to final state particles
          try:
            proc = 'DM_particle'+str(i+1)+' DM_particle'+str(j+1)+' > dm_particles dm_particles'
            self._MG5Cmd.do_generate(proc)

            # Get the matrix elements and make sure that we don't have any pure scattering processes
            curr_matrix_elements = helas_objects.HelasMultiProcess(self._MG5Cmd._curr_amps)
            self._dm2dm_me[i][j] = curr_matrix_elements.get_matrix_elements()
            self._wanted_lorentz += curr_matrix_elements.get_used_lorentz()
            self._wanted_couplings += curr_matrix_elements.get_used_couplings()
            for me in self._dm2dm_me[i][j]:
              if (set(me.get('processes')[0].get_initial_ids()) == (set(me.get('processes')[0].get_final_ids()))):
                self._dm2dm_me[i][j].remove(me)
          except:
            continue

      # Generate all the DM particles scattering off of the thermal background and 
      # change to a different DM species (again, no pure scatterings)
      for i in range(len(self._dm_particles)):
        for j in range(len(self._dm_particles)):

          # We do not want to genereate the purely scattering processes
          if (i != j):

            # Create the appropriate strong of initial to final state particles
            try:
              proc = 'DM_particle'+str(i+1)+' fs_particles > DM_particle'+str(j+1)+' fs_particles'
              self._MG5Cmd.do_generate(proc)

              # Get the matrix elements
              curr_matrix_elements = helas_objects.HelasMultiProcess(self._MG5Cmd._curr_amps)
              self._scattering_me[i][j] = curr_matrix_elements.get_matrix_elements()
              self._wanted_lorentz += curr_matrix_elements.get_used_lorentz()
              self._wanted_couplings += curr_matrix_elements.get_used_couplings()
            except:
              continue

      # Look at which particles have the thermal scattering diagrams and those that do
      # get flagged so that the fortran side will know how to properly handle the numerical code
      for i in range(len(self._dm_particles)):
        sum_thermal_scattering = 0
        for j in range(len(self._dm_particles)):
          sum_thermal_scattering += len(self._scattering_me[i][j])
        if (sum_thermal_scattering == 0):
          self._dm_thermal_scattering.append(False)
        else:
          self._dm_thermal_scattering.append(True)

      # Since we are done generating the matrix elements we convert the names so the 
      # particle names can be used in the fortran code. 
      for i in range(len(self._dm_particles)):
        self._dm_names_list[i] = self.Convertname(self._dm_names_list[i])
        self._dm_antinames_list[i] = self.Convertname(self._dm_antinames_list[i])
      for i in range(len(self._DM_all_names)):
        self._DM_all_names[i] = self.Convertname(self._DM_all_names[i])

      #print "Finished!\n"



    #-----------------------------------------------------------------------#
    def CreateNumericalSession(self, prompts = True):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  This routine creates several final that are needed to run the        #
    #  FORTRAN side of the code.  The files that are created include:       #
    #                                                                       #
    #  - all the individual matrix element fortran files.                   #
    #  - all the individual pmass include files associated with each        #
    #    matrix element                                                     #
    #  - 
    #  code.                                                                #
    #                                                                       #
    #-----------------------------------------------------------------------#

      # If the project already exists just move over the param card
      if (self._project_exists):
        # Move over the param_card to the project directory
        shutil.move('Projects/param_card.dat',self._projectpath+'/Cards/param_card.dat')
        self._paramcard = os.path.join(self._projectpath,'Cards/param_card.dat') #MIHAILO ADDED THIS
        
    	#Dump the darkmatter object so it can be read back in by the scrip
        pickle.dump(self, open( self._projectpath+'/dm_object.pik', 'wb' ), -1)     
        
        return

      # Otherwise create a new numerical session    

      # Set up the export class
      self._exporter = MGoutput.ProcessExporterFortranMadDM(self._mgme_dir, self._projectpath)
      self._exporter.opt['model'] = self._modelname

      # Copy over the template directory 
      self._exporter.copy_template(self._projectpath)
      
      # Move over the param_card to the project directory
      shutil.move('Projects/param_card.dat',self._projectpath+'/Cards/param_card.dat')
      self._paramcard = os.path.join(self._projectpath,'Cards/param_card.dat') #MIHAILO ADDED THIS
      
      #Copy the parameter scan default script 
      shutil.copy('param_scan_default.py',self._projectpath+'/param_scan_default.py') #MIHAILO ADDED THIS
      #Dump the darkmatter object so it can be read back in by other Python scripts. 
      print "PP: "+self._projectpath
      pickle.dump(self, open( self._projectpath+'/dm_object.pik', 'wb' ), -1)  #MIHAILO ADDED THIS

      # Arrays that are needed to store the information for the annihilation diagrams
      self._ann_nprocesses = [[0 for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
      self._ann_process_names = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
      self._ann_process_iden_init = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]

      # Arrays that are needed to store the information for the DM -> DM diagrams
      self._dm2dm_nprocesses = [[0 for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
      self._dm2dm_process_names = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
      self._dm2dm_process_iden_init = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
      self._dm2dm_final_states = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]

      # Arrays that are needed to store the information for the DM SM -> DM SM processes
      self._scattering_nprocesses = [[0 for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
      self._scattering_process_names = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
      self._scattering_initial_dofs = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
      self._scattering_initial_dofs_total = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]

      # writing the matrix elements
      path_matrix = os.path.join(self._projectpath, 'matrix_elements')

      # Annihilation matrix elements
      for i in range(len(self._dm_particles)):
        for j in range(i,len(self._dm_particles)):

          # Get the total number of annihilation processes
          self._ann_nprocesses[i][j] = len(self._annihilation_me[i][j])

          for me in self._annihilation_me[i][j]:

            # Get the name of the process (we disregard the first two characters in the name '0_')
            process_name = me.get('processes')[0].shell_string()            
            process_name = process_name[2:len(process_name)]
            self._ann_process_names[i][j].append(process_name)

            # Check to see if the initial state particles are identical
            initial_state = me.get('processes')[0].get_initial_ids()
            if (initial_state[0] == initial_state[1]):
              self._ann_process_iden_init[i][j].append(True)
            else:
              self._ann_process_iden_init[i][j].append(False)

            # Using the proess name we create the filename for each process and export the fortran file
            filename_matrix = os.path.join(path_matrix, 'matrix_' + process_name + '.f')
            self._exporter.write_matrix_element(\
                    writers.FortranWriter(filename_matrix),\
                    me, self._MG5Cmd._curr_fortran_model)

            # We also need the external mass include file for each process.
            filename_pmass = os.path.join(self._projectpath, 'include', 'pmass_' + process_name + '.inc')
            self._exporter.write_pmass_file(\
                    writers.FortranWriter(filename_pmass), me)


      # DM -> DM matrix elements
      for i in range(len(self._dm_particles)):
        for j in range(i,len(self._dm_particles)):

          # Get the total number of dm2dm processes
          self._dm2dm_nprocesses[i][j] = len(self._dm2dm_me[i][j])

          for me in self._dm2dm_me[i][j]:

            # Get the name of the process (we disregard the first two characters in the name '0_')
            process_name = me.get('processes')[0].shell_string()            
            process_name = process_name[2:len(process_name)]
            self._dm2dm_process_names[i][j].append(process_name)

            # Check to see if the initial state particles are identical
            initial_state = me.get('processes')[0].get_initial_ids()
            if (initial_state[0] == initial_state[1]):
              self._dm2dm_process_iden_init[i][j].append(True)
            else:
              self._dm2dm_process_iden_init[i][j].append(False)

            # Get the final states of the DM -> DM processes
            final_state = me.get('processes')[0].get_final_ids()
            # Change the PDG codes to the index in the DM particle list
            for k in range(len(self._dm_particles)):
              if abs(final_state[0]) == abs(self._dm_particles[k]['pdg_code']):
                final_state[0] = k+1
              if abs(final_state[1]) == abs(self._dm_particles[k]['pdg_code']):
                final_state[1] = k+1
            self._dm2dm_final_states[i][j].append(final_state)            

            # Using the proess name we create the filename for each process and export the fortran file
            filename_matrix = os.path.join(path_matrix, 'matrix_' + process_name + '.f')
            self._exporter.write_matrix_element(\
                    writers.FortranWriter(filename_matrix),\
                    me, self._MG5Cmd._curr_fortran_model)

            # We also need the external mass include file for each process.
            filename_pmass = os.path.join(self._projectpath, 'include', 'pmass_' + process_name + '.inc')
            self._exporter.write_pmass_file(\
                    writers.FortranWriter(filename_pmass), me)


      # SM DM -> SM DM matrix elements
      for i in range(len(self._dm_particles)):
        for j in range(len(self._dm_particles)):

          # Add to the total number of dm2dm processes
          self._scattering_nprocesses[i][j] = len(self._scattering_me[i][j])

          for me in self._scattering_me[i][j]:

            # Get the name of the process (we disregard the first two characters in the name '0_')
            process_name = me.get('processes')[0].shell_string()            
            process_name = process_name[2:len(process_name)]
            self._scattering_process_names[i][j].append(process_name)

            # Get the number of degrees of freedom for the SM particle in the initial state
            initial_state = me.get('processes')[0].get_initial_ids()
            for k in range(len(self._particles)):
              if (abs(initial_state[1]) == self._particles[k]['pdg_code']):
                initial_SM_dof = self._particles[k]['spin']*self._particles[k]['color']
                self._scattering_initial_dofs[i][j].append(initial_SM_dof)
                if (self._particles[k]['self_antipart']):
                  self._scattering_initial_dofs_total[i][j].append(initial_SM_dof)
                else:
                  self._scattering_initial_dofs_total[i][j].append(2*initial_SM_dof)

            # Using the proess name we create the filename for each process and export the fortran file
            filename_matrix = os.path.join(path_matrix, 'matrix_' + process_name + '.f')
            self._exporter.write_matrix_element(\
                    writers.FortranWriter(filename_matrix),\
                    me, self._MG5Cmd._curr_fortran_model)

            # We also need the external mass include file for each process.
            filename_pmass = os.path.join(self._projectpath, 'include', 'pmass_' + process_name + '.inc')
            self._exporter.write_pmass_file(\
                    writers.FortranWriter(filename_pmass), me)


      # Create the dm_info.inc file that contains all the DM particles as well as spin and mass information
      self.WriteDMInfo()

      # Create the model_info.txt file that is used to calculate the number of relativistic degrees of freedom.
      self.WriteModelInfo()

      # Create the diagrams.inc file that contains the number of processes for each pair of initial state particles
      self.WriteDiagramInfo()

      # Create the process_names.inc file that contains the names of all the individual processes
      self.WriteProcessNames()
 
      # Create the smatrix.f file
      self.Write_smatrix()

      # Create the makefile for compiling all the matrix elements
      self.Write_makefile()

      # This creates the fortran files for the model
      # v4 model
      if self._MG5Cmd._model_v4_path:
        print 'Copy %s model files to directory %s' % (os.path.basename(self._model_v4_path), self._projectpath)
        self._exporter.export_model_files(self._model_v4_path)
        self._exporter.export_helas(pjoin(self._mgme_dir,'HELAS'))
      # v5 model
      else:
        print 'Export UFO model to MG4 format'
        self._exporter.convert_model_to_mg4(self._Model, self._wanted_lorentz, self._wanted_couplings)

      # Copy over the coupl.inc and input.inc files to the project's include directory.
      shutil.copy2(self._projectpath+'/Source/MODEL/coupl.inc', self._projectpath + '/include/coupl.inc')
      shutil.copy2(self._projectpath+'/Source/MODEL/input.inc', self._projectpath + '/include/input.inc')

	  #Edit the maddm_card.inc file
      if prompts == True:
         print "Would you like to edit maddm_card.inc?[n] (y/n):"
         answer = raw_input()
         if answer == 'y' or answer == 'Y':
            subprocess.call('vi '+self._projectpath+'/include/maddm_card.inc', shell=True)
      
      # Now that everything is created we can compile the whole code            
      #print "------Compiling the numerical code------"
      curr_dir = os.getcwd() # get current directory
      makedir = os.path.join(curr_dir, self._projectpath)        
      subprocess.call(['make'],cwd = makedir)

      #print "Finished!\n"


    #-----------------------------------------------------------------------#
    def WriteDMInfo(self):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  This function writes the information about the inital state particle #
    #  dof to a file for use by the FORTRAN part of the code.               #
    #                                                                       #
    #-----------------------------------------------------------------------#

      filename_dof = os.path.join('Projects',self._projectname,'include','dm_info.inc')
      file_writer = open(filename_dof, 'w')

      # Write out some header comments for the file
      file_writer.write("c------------------------------------------------------------------------------c\n")
      file_writer.write("c This file contains the needed information about all the DM particles.\n")

      # The counter is used to write the correct array index needed for the fortran code
      counter = 1
      for dm_particle in self._dm_particles:

        dof = float(dm_particle['spin'])*float(dm_particle['color'])

        # Write out the mass parameter, degrees of freedom, particle name, and the index for the end of the name
        file_writer.write("c------------------------------------------------------------------------------c\n")
        file_writer.write("      mdm(" + str(counter) + ") = abs(" + dm_particle['mass'] + ")\n")
        file_writer.write("      dof_dm(" + str(counter) +") = " + str(dof) + "\n")
        if (self._dm_thermal_scattering[counter-1]):
          file_writer.write("      dm_sm_scattering(" + str(counter) + ") = .true.\n")
        else:
          file_writer.write("      dm_sm_scattering(" + str(counter) + ") = .false.\n")
        file_writer.write("      dm_names(" + str(counter) + ") = \'" + self._dm_names_list[counter-1] + "\'\n")
        file_writer.write("      dm_index(" + str(counter) + ") = " + str(len(dm_particle['name'])) + "\n")
        file_writer.write("      dm_antinames(" + str(counter) + ") = \'" + self._dm_antinames_list[counter-1] + "\'\n")
        file_writer.write("      dm_antiindex(" + str(counter) + ") = " + str(len(dm_particle['antiname'])) + "\n") 

        counter = counter + 1

      file_writer.write("c------------------------------------------------------------------------------c\n")
      file_writer.close()



    #-----------------------------------------------------------------------#
    def WriteModelInfo(self):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  This routine writes the mass, degrees of freedom and boson/fermion   #
    #  information for all the relevant particles in the model so that the  #
    #  number of relativistic degrees of freedom can be calculated.         #
    #                                                                       #
    #  The SM particles are in the template by default so this would just   #
    #  add the relevant BSM particles.                                      #
    #                                                                       #
    #-----------------------------------------------------------------------#

      # Flags that are needed to write the model info file
      __flag1 = '#NUM_PARTICLES'
      __flag2 = '#BSM_INFO'

      # Open up the template and output file
      model_info_file = open('Projects/'+self._projectname+'/include/model_info.txt','w')
      model_info_template = open('Projects/'+self._projectname+'/include/model_info_template.txt','r')
      model_info_lines = model_info_template.readlines()

      # Read all the lines from the template file and insert appropriate parts
      # which are flagged with __flag1 and __flag2 
      for line in model_info_lines:

        # write out the number of DM paticles
        if __flag1 in line:
          model_info_file.write(str(17 + len(self._dm_particles) + len(self._bsm_final_states))+'\n')       

        # write out the mass, degrees of freedom and fermion/boson info for the bsm particles
        elif __flag2 in line:
          bsm_particles = self._dm_particles + self._bsm_final_states
          for bsm_particle in bsm_particles:

            # Get all the necessary information
            mass = self.GetMass(bsm_particle['pdg_code'])
            dof = float(bsm_particle['spin'])*float(bsm_particle['color'])

            # If the particle has an anti particle double the number of degrees of freedom
            if (not bsm_particle['self_antipart']):
              dof *= 2.0

            # The spin information is written in 2s+1 format (i.e. 1, 3, etc = boosn; 2, 4, etc = fermion)
            if (int(bsm_particle['spin']) == (int(bsm_particle['spin']/2)*2)):
              boson_or_fermion = 1
            else:
              boson_or_fermion = 0

            # write the line to the file
            new_line = str(mass)+'   '+str(dof)+'   '+str(boson_or_fermion)+'\n'
            model_info_file.write(new_line)

        # If there is no flag then it's a SM particle which we just write to the output file
        else:
          model_info_file.write(line)

      model_info_file.close()



    #-----------------------------------------------------------------------#
    def WriteDiagramInfo(self):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  This routine creates the diagrams.inc file which contains the number #
    #  of annihilation diagrams for each pair of DM particles.              #
    #                                                                       #
    #-----------------------------------------------------------------------#

      # open up the file and first print out the number of dm particles
      diagramsfile = open('Projects/'+self._projectname+'/include/diagrams.inc','w')
      stringtowrite = 'c Total number of IS particles participating in the coannihilations \n'+\
                          '      ndmparticles = '+str(len(self._dm_particles))+'\n\n' +\
                          'c Processes by class'
      diagramsfile.write(stringtowrite)

      # Write out how many annihlation processes for each combinations of initial states
      diagramsfile.write('\nc Number of annihilation processes for each DM pair\n')
      for i in range(len(self._dm_particles)):
        for j in range(i,len(self._dm_particles)):
          diagramsfile.write('      ann_nprocesses('+str(i+1)+','+str(j+1)+') = '+str(self._ann_nprocesses[i][j])+'\n')

      # Write out how many DM -> DM processes for each combinations of initial states
      diagramsfile.write('\nc Number of DM -> DM processes for each DM pair\n')
      for i in range(len(self._dm_particles)):
        for j in range(i,len(self._dm_particles)):
          diagramsfile.write('      dm2dm_nprocesses('+str(i+1)+','+str(j+1)+') = '+str(self._dm2dm_nprocesses[i][j])+'\n')

      # Write out how many DM/SM scattering processes for each combinations of initial states
      diagramsfile.write('\nc Number of DM/SM scattering processes for each DM pair\n')
      for i in range(len(self._dm_particles)):
        for j in range(len(self._dm_particles)):
          stringtowrite = '      scattering_nprocesses('+str(i+1) +','+str(j+1) +') = '+\
              str(self._scattering_nprocesses[i][j])+'\n'
          diagramsfile.write(stringtowrite)

      diagramsfile.close()



    #-----------------------------------------------------------------------#
    def WriteProcessNames(self):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  This routine creates the process_names.inc file which contains the   #
    #  names of all the individual processes.  These names are primairly    #
    #  used in the test subroutines on the fortran side.                    #
    #                                                                       #
    #-----------------------------------------------------------------------#

      # This creates the file process_names.inc which will have a list of all the individual process names.
      process_names_file = open('Projects/'+self._projectname+'/include/process_names.inc','w')

      # Write out the list of process names
      process_names_file.write('c List of the process names in order of dmi, dmj, nprocesses(dmi, dmj)\n')
      process_counter = 0

      # Annihilation process names
      process_names_file.write('c Annihilation process names\n')
      self._total_annihilation_nprocesses = 0
      for i in range(len(self._dm_particles)):
        for j in range(i,len(self._dm_particles)):
          self._total_annihilation_nprocesses += self._ann_nprocesses[i][j]
          for k in range(self._ann_nprocesses[i][j]):

            # Increment the process counters and write the process name
            process_counter += 1
            process_names_file.write('      process_names('+str(process_counter)+') = \''+self._ann_process_names[i][j][k]+'\'\n')

      # DM -> DM process names
      process_names_file.write('\nc DM -> DM process names\n')
      self._total_dm2dmscattering_nprocesses = 0
      for i in range(len(self._dm_particles)):
        for j in range(i,len(self._dm_particles)):
          self._total_dm2dmscattering_nprocesses += self._dm2dm_nprocesses[i][j]
          for k in range(self._dm2dm_nprocesses[i][j]):

            # Increment the process counter and write the process name
            process_counter += 1
            process_names_file.write('      process_names('+str(process_counter)+') = \''+self._dm2dm_process_names[i][j][k]+'\'\n')

      # DM/SM scattering process names
      process_names_file.write('\nc DM/SM scattering process names\n')
      self._total_scattering_nprocesses = 0
      for i in range(len(self._dm_particles)):
        for j in range(len(self._dm_particles)):
          self._total_scattering_nprocesses += self._scattering_nprocesses[i][j]
          for k in range(self._scattering_nprocesses[i][j]):

            # Increment the process counter and write the process name
            process_counter += 1
            process_names_file.write('      process_names('+str(process_counter)+') = \''+\
                  self._scattering_process_names[i][j][k]+'\'\n')


      # Write out the total number of process names and number of each group of processes
      process_names_file.write('\nc Total number of processes for each category\n')
      process_names_file.write('      num_processes = ' + str(process_counter) + '\n')
      process_names_file.write('      ann_num_processes = ' + str(self._total_annihilation_nprocesses) + '\n')
      process_names_file.write('      dm2dm_num_processes = ' + str(self._total_dm2dmscattering_nprocesses) + '\n')
      process_names_file.write('      scattering_num_processes = ' + str(self._total_scattering_nprocesses) + '\n')


      # Write out the boolean flags for the processes with identical initial state particles
      process_names_file.write('\nc Boolean operators for identical particles in the initial state\n')
      process_names_file.write('c Annihialtion diagrams\n')
      for i in range(len(self._dm_particles)):
        for j in range(i,len(self._dm_particles)):

          # Boolean operators for the annihilation processes
          for k in range(self._ann_nprocesses[i][j]):

            if (self._ann_process_iden_init[i][j][k]):
              process_names_file.write('      ann_process_iden_init('+str(i+1)+','+str(j+1)+','+str(k+1)+') = .true.\n')
            else:
              process_names_file.write('      ann_process_iden_init('+str(i+1)+','+str(j+1)+','+str(k+1)+') = .false.\n')

      process_names_file.write('\nc DM -> DM diagrams\n')
      for i in range(len(self._dm_particles)):
        for j in range(i, len(self._dm_particles)):

          # Boolean operators for the dm2dm processes
          for k in range(self._dm2dm_nprocesses[i][j]):

            if (self._dm2dm_process_iden_init[i][j][k]):
              process_names_file.write('      dm2dm_process_iden_init('+str(i+1)+','+str(j+1)+','+str(k+1)+') = .true.\n')
            else:
              process_names_file.write('      dm2dm_process_iden_init('+str(i+1)+','+str(j+1)+','+str(k+1)+') = .false.\n')


      # Write out final state information for all the DM -> DM processes
      process_names_file.write('\nc Final state information for all the DM -> DM processes\n')
      for i in range(len(self._dm_particles)):
        for j in range(i,len(self._dm_particles)):
          for k in range(self._dm2dm_nprocesses[i][j]):

            process_names_file.write('      dm2dm_fs('+str(i+1)+','+str(j+1)+','+str(k+1)+',1) = '+\
                  str(self._dm2dm_final_states[i][j][k][0])+'\n')
            process_names_file.write('      dm2dm_fs('+str(i+1)+','+str(j+1)+','+str(k+1)+',2) = '+\
                  str(self._dm2dm_final_states[i][j][k][1])+'\n')

      # Write out the inital state degrees of freedom for all the DM/SM scattering processes
      process_names_file.write('\nc Initial state degrees of freedom for all the DM/SM scattering processes\n')
      for i in range(len(self._dm_particles)):
        for j in range(len(self._dm_particles)):
          for k in range(self._scattering_nprocesses[i][j]):

            process_names_file.write('      dof_SM('+str(i+1)+','+str(j+1)+','+str(k+1)+') = '+\
                  str(self._scattering_initial_dofs[i][j][k])+'\n')

      # Write out the total inital state degrees of freedom for all the DM/SM scattering processes
      process_names_file.write('\nc Total initial state degrees of freedom for all the DM/SM scattering processes\n')
      for i in range(len(self._dm_particles)):
        for j in range(len(self._dm_particles)):
          for k in range(self._scattering_nprocesses[i][j]):

            process_names_file.write('      dof_SM_total('+str(i+1)+','+str(j+1)+','+str(k+1)+') = '+\
                  str(self._scattering_initial_dofs_total[i][j][k])+'\n')


    #-----------------------------------------------------------------------#
    def Write_smatrix(self):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  This routine creates the smatrix.f file based on the smatrix         #
    #  template.  This is used to call the appropriate matrix element as    #
    #  well as the appropriate pmass include file for each process.         #
    #                                                                       #
    #-----------------------------------------------------------------------#

      # Flags that are needed to create the matrix element makefile and smatrix.f files        
      __flag1 = '#PMASS_MADDM1'
      __flag2 = '#PMASS_MADDM2'
      __flag3 = '#PMASS_MADDM3'
      __flag4 = '#SMATRIX_MADDM1'
      __flag5 = '#SMATRIX_MADDM2'
      __flag6 = '#SMATRIX_MADDM3'

      # Edit smatrix.f file to incorporate all different subprocess
      smatrix_file = open('Projects/'+self._projectname+'/matrix_elements/smatrix.f','w')
      smatrix_template = open('Projects/'+self._projectname+'/matrix_elements/smatrix_template.f','r')
      smatrix_lines = smatrix_template.readlines()
        
      # Read all the lines from the smatrix_template.f and insert appropriate parts
      # which are flagged with __flag1 and __flag2 in the smatrix_template.f file.
      # Write the output into smatrix.f
      for line in smatrix_lines:
       
        # Include all the combinations of external pmass.inc files for the annihlation processes
        if __flag1 in line:
          firstentry = True
          for i in range(len(self._dm_particles)):
            for j in range(i,len(self._dm_particles)):
              for k in range(self._ann_nprocesses[i][j]):
                if (firstentry):
                  firstentry = False
                  new_line = '      if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                      '(k.eq.'+str(k+1)+')) then\n        include \'pmass_'+self._ann_process_names[i][j][k]+'.inc\' \n'
                  smatrix_file.write(new_line)
                else:
                  new_line = '      else if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                      '(k.eq.'+str(k+1)+')) then\n        include \'pmass_'+self._ann_process_names[i][j][k]+'.inc\' \n'
                  smatrix_file.write(new_line)
          if (self._total_annihilation_nprocesses > 0):
            smatrix_file.write('      endif')

        # Include all the combinations of external pmass.inc files for the DM -> DM processes
        elif __flag2 in line:
          firstentry = True
          for i in range(len(self._dm_particles)):
            for j in range(i,len(self._dm_particles)):
              for k in range(self._dm2dm_nprocesses[i][j]):
                if (firstentry):
                  firstentry = False
                  new_line = '      if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                      '(k.eq.'+str(k+1)+')) then\n        include \'pmass_'+self._dm2dm_process_names[i][j][k]+'.inc\' \n'
                  smatrix_file.write(new_line)
                else:
                  new_line = '      else if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                      '(k.eq.'+str(k+1)+')) then\n        include \'pmass_'+self._dm2dm_process_names[i][j][k]+'.inc\' \n'
                  smatrix_file.write(new_line)
          if (self._total_dm2dmscattering_nprocesses > 0):
            smatrix_file.write('      endif')

        # Include all the combinations of external pmass.inc files for the DM/SM scattering processes
        elif __flag3 in line:
          firstentry = True
          for i in range(len(self._dm_particles)):
            for j in range(len(self._dm_particles)):
              for k in range(self._scattering_nprocesses[i][j]):
                if (firstentry):
                  firstentry = False
                  new_line = '      if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                      '(k.eq.'+str(k+1)+')) then\n        include \'pmass_'+self._scattering_process_names[i][j][k]+'.inc\' \n'
                  smatrix_file.write(new_line)
                else:
                  new_line = '      else if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                      '(k.eq.'+str(k+1)+')) then\n        include \'pmass_'+self._scattering_process_names[i][j][k]+'.inc\' \n'
                  smatrix_file.write(new_line)
          if (self._total_scattering_nprocesses > 0):
            smatrix_file.write('      endif')

        # Creates all the appropriate calls to the individual annihilaton smatrix.f's
        elif __flag4 in line:
          firstentry = True
          for i in range(len(self._dm_particles)):
            for j in range(i,len(self._dm_particles)):
              for k in range(self._ann_nprocesses[i][j]):
                if (firstentry):
                  firstentry = False
                  new_line = '      if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                      '(k.eq.'+str(k+1)+')) then\n        call smatrix_'+\
                       self._ann_process_names[i][j][k]+'(p_ext,smatrix_ann)\n'
                  smatrix_file.write(new_line)
                else:
                  new_line = '      else if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                      '(k.eq.'+str(k+1)+')) then\n        call smatrix_'+\
                       self._ann_process_names[i][j][k]+'(p_ext,smatrix_ann)\n'
                  smatrix_file.write(new_line)
          if (self._total_annihilation_nprocesses > 0):
            smatrix_file.write('      endif')

        # Creates all the appropriate calls to the individual annihilaton smatrix.f's
        elif __flag5 in line:
          firstentry = True
          for i in range(len(self._dm_particles)):
            for j in range(i,len(self._dm_particles)):
              for k in range(self._dm2dm_nprocesses[i][j]):
                if (firstentry):
                  firstentry = False
                  new_line = '      if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                      '(k.eq.'+str(k+1)+')) then\n        call smatrix_'+\
                      self._dm2dm_process_names[i][j][k]+'(p_ext,smatrix_dm2dm)\n'
                  smatrix_file.write(new_line)
                else:
                  new_line = '      else if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                      '(k.eq.'+str(k+1)+')) then\n        call smatrix_'+\
                      self._dm2dm_process_names[i][j][k]+'(p_ext,smatrix_dm2dm)\n'
                  smatrix_file.write(new_line)
          if (self._total_dm2dmscattering_nprocesses > 0):
            smatrix_file.write('      endif')

        # Creates all the appropriate calls to the individual annihilaton smatrix.f's
        elif __flag6 in line:
          firstentry = True
          for i in range(len(self._dm_particles)):
            for j in range(len(self._dm_particles)):
              for k in range(self._scattering_nprocesses[i][j]):
                if (firstentry):
                  firstentry = False
                  new_line = '      if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                      '(k.eq.'+str(k+1)+')) then\n        call smatrix_'+\
                      self._scattering_process_names[i][j][k]+'(p_ext,smatrix_scattering) \n'
                  smatrix_file.write(new_line)
                else:
                  new_line = '      else if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                      '(k.eq.'+str(k+1)+')) then\n        call smatrix_'+\
                      self._scattering_process_names[i][j][k]+'(p_ext,smatrix_scattering) \n'
                  smatrix_file.write(new_line)
          if (self._total_scattering_nprocesses > 0):
            smatrix_file.write('      endif')

        # If none of the flags are present then just write the line from the template to the file
        else:
	        smatrix_file.write(line)

      smatrix_template.close()
      smatrix_file.close()



    #-----------------------------------------------------------------------#
    def Write_makefile(self):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  This routine creates the makefile for compiling all the matrix       #
    #  elements generated by madgraph.                                      #
    #                                                                       #
    #-----------------------------------------------------------------------#

      # Flags that are needed to create the matrix element makefile and smatrix.f files        
      __flag = '#MAKEFILE_MADDM'

      # Creates the make file in the matrix_elements folder.
      # includes the object files for all the individual matrix elements
      makefile = open('Projects/'+self._projectname+'/matrix_elements/makefile','w')
      makefile_template = open('Projects/'+self._projectname+'/matrix_elements/makefile_template','r')

      makefile_lines = makefile_template.readlines()
      for line in makefile_lines:
        if __flag in line:
          new_line = 'objs = smatrix.o'
          # Add all the annihilation matrix.o files to the makefile
          for i in range(len(self._dm_particles)):
            for j in range(i,len(self._dm_particles)):
              for k in range(self._ann_nprocesses[i][j]):
                new_line = new_line+' matrix_'+self._ann_process_names[i][j][k]+'.o'

          # Add all the DM -> DM matrix.o files to the makefile
          for i in range(len(self._dm_particles)):
            for j in range(i,len(self._dm_particles)):
              for k in range(self._dm2dm_nprocesses[i][j]):
                new_line = new_line+' matrix_'+self._dm2dm_process_names[i][j][k]+'.o'

          # Add all the DM/SM scattering matrix.o files to the makefile
          for i in range(len(self._dm_particles)):
            for j in range(len(self._dm_particles)):
              for k in range(self._scattering_nprocesses[i][j]):
                new_line = new_line+' matrix_'+self._scattering_process_names[i][j][k]+'.o'

          new_line = new_line + '\n'
          makefile.write(new_line)
        else:
          makefile.write(line)
       
      makefile.close()



    #-----------------------------------------------------------------------#
    def CalculateRelicAbundance(self):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  Runs the maddm.x code created by CreateNumericalSession()            #
    #                                                                       #
    #-----------------------------------------------------------------------#

      cwd = os.getcwd()
      cmd = self._projectpath+'/maddm.x > omega' #'Projects/'+self._projectname+'/maddm.x > omega'
      exec_dir = self._projectpath
      exec_cmd = os.path.join(cwd, cmd)
      subprocess.call(exec_cmd, cwd = exec_dir, shell=True)
      
      #result = open('Projects/'+self._projectname+'/omega', 'r')
      result = open(self._projectpath+'/omega', 'r')
      omegah2 = result.read()
      result.close()
      return omegah2
      

    #-----------------------------------------------------------------------#
    def ChangeParameter(self, flag, value):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  Changes the parameter in the param card                              #
    #  flag is a string name of the parameter. Value is the numerical value #
    #  that the parameter should be set to.                                 #
    #                                                                       #
    #-----------------------------------------------------------------------#

       try:
        
        #Change the CM energy
        foundflag = False
        for line in fileinput.FileInput(self._paramcard,\
                inplace=1):
            line_list = line.split()
            if flag in line_list:
            	#Make sure that the flag is not a part of another variable name
            	foundflag = True
                oldvalue = line_list[1]
                newline=line.replace(oldvalue,str(value))
                print newline,
            else: 
                print line,

        fileinput.close()
        if not foundflag:
        	print "Error: Parameter not found!"

       except OSError:
        print "Error: Could not open the Param Card file!"
        exit

    #-----------------------------------------------------------------------#
    def init_param_scan(self,name_of_script):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  Starts the param_scan.py program. The function prompts the user to   #
    #  edit the param_scan.py script and set up the parameter scan.         #
    #   Editing the param_scan.py file requires basic knowledge of Python   #
    #-----------------------------------------------------------------------#

		#Edit the parameter scan script
		shutil.copy2('param_scan_default.py', self._projectpath+'/'+name_of_script)
		subprocess.call('vi '+self._projectpath+'/%s' % name_of_script, shell= True)

		curr_dir = os.getcwd()
        	os.chdir(self._projectpath)

        	cmd = os.path.join('./',name_of_script)
        	subprocess.call(cmd, shell=True)
		os.chdir(curr_dir)
		
