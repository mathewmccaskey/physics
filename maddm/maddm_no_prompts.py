#! /usr/bin/env python

from init import *
from darkmatter import *

#Start the user interface
#[new_proj, model_name, project_name] = initialize_maddm_session(True)

# Create an instance of the base darkmatter class and initialize it
# with the user input
# if the project is already existent, do not generate diagrams etc.
dm = darkmatter()
new_project = True
dm.init_from_model('rsxSM', 'rsxSM', new_project)


#If it is a new project...
if new_project:

	ans=raw_input('Warning: The code will erase the existing Project folder. Would you like to continue?[n] (y/n):')
	if ans !='y':
		exit()

	# Find the DM candidate and coannihilaton particles
	dm.FindDMCandidate(prompts = False)
	dm.FindCoannParticles(prompts = False, coann_eps = 0.2)
	dm.GetProjectName()

	print "------ Generating Diagrams ------",      
	dm.GenerateDiagrams()    
	print " Done!"
	print "------ Creating the Numerical Session ------",
	dm.CreateNumericalSession(prompts = False)
	print " Done!"
	
	print "Diagnostics:"
	print dm._projectname
	print dm._paramcard
	print dm._projectpath
	
	
	print "------ Calculating Relic Abundance ------ ",
	oh2 = dm.CalculateRelicAbundance()
	print " Done!"
	print "\n     Omega h^2 = "+str(oh2)

