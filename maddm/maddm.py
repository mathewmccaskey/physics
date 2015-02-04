#! /usr/bin/env python

from init import *
from darkmatter import *

#Start the user interface
[new_proj, model_name, project_name] = initialize_MadDM_session(True)

# Create an instance of the base darkmatter class and initialize it
# with the user input
# if the project is already existent, do not generate diagrams etc.
dm = darkmatter()
dm.init_from_model(model_name, project_name, new_proj)


#If it is a new project...
if new_proj == True:
	# Find the DM candidate and coannihilaton particles
	dm.FindDMCandidate()
	dm.FindCoannParticles()
	dm.GetProjectName()

	print "------ Generating Diagrams ------\n",      
	dm.GenerateDiagrams()    
	print " Done!"
	print "------ Creating the Numerical Session ------\n",
	dm.CreateNumericalSession()
	print " Done!"
	
	print "Diagnostics:"
	print dm._projectname
	print dm._paramcard
	print dm._projectpath
	
	
	print "------ Calculating Relic Abundance ------ ",
	oh2 = dm.CalculateRelicAbundance()
	print " Done!"
	print "\n     Omega h^2 = "+str(oh2)

#If a project already exists skip to the optional parameter scan.
legit_answer = False 
while not legit_answer:
	do_param_scan = raw_input('Would you like to perform a parameter scan?[n] (y/n):')
	if do_param_scan == 'y' or do_param_scan == 'Y':
		legit_answer = True
		param_scan_name = raw_input('Enter the name of your parameter scan: ')
		dm.init_param_scan(param_scan_name)
						
	elif do_param_scan == 'n' or do_param_scan == 'N' or do_param_scan =='':
		legit_answer = True
		print "Finished! Exiting!"
		exit(0)
	else:
		print "Not a legitimate input. Try again."