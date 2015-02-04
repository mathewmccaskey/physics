#! /usr/bin/env python

import numpy as np
import fileinput
import sys
import os
import pickle

#Set the proper paths 
maddm_path = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]
pieces = maddm_path.split('/')
maddm_path = maddm_path.replace(pieces[-1], '')
sys.path.append(maddm_path)

from init import *
from darkmatter import *

#Read in the DM information
with open('dm_object.pik', 'r') as f:
	dm = pickle.load(f)

#Change the directory to MadDM root folder so that code can run properly.
os.chdir(maddm_path)

#Print out some basic information about the dark matter particle.
print "--------------------------------------------"
print "Model: "+dm._modelname
print "Project Name: "+dm._projectname
print "Parameter Card: "+dm._paramcard
print "DM particles: "+dm._dm_particles[0].get('name')
print "DM spin: "+str(dm._dm_particles[0].get('spin'))
print "--------------------------------------------"

#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#CHANGE THIS PART FOR YOUR PARAMETER SCAN
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#Add as many for loops as you wish to scan over parameters
#NOTE:  The default script varies over two parameters (param1 and param2), and outputs the
#		results in the format "param1    param2    omegah^2" on the screen and in the output file.
#		Please make sure to change the number of "for" loops, calls to ChangeParameter() and the
#		output commands to account for the number of parameters you wish to vary. The critical
#		places where changes are usually required are marked with "<-----Change here ----!!!"


#Define the arrays of values your parameters should take.
#If you are using np.arrange, the format is np.arange(init_val, final_val, step_size)
param1_values = np.arange(0.001, 0.501, 0.01) #"<-----Change here ----!!!"
param2_values = np.arange(0.001, 0.501, 0.01) #"<-----Change here ----!!!"

outputfile = open("output_file.txt", "w") #"<-----Change here ----!!!"

for param1 in param1_values:
	for param2 in param2_values: #"<-----Change here ----!!!" (for number of loops)
		#Change the parameter with name 'parameter1' in the param_card.dat to the value param1
		dm.ChangeParameter('parameter1', param1) #"<-----Change here ----!!!"
		#Change the parameter with name 'parameter2' in the param_card.dat to the value param2
		dm.ChangeParameter('parameter2', param2) #"<-----Change here ----!!!"
		#Calculate relic density
		omega = dm.CalculateRelicAbundance()
		
		#Output the results on the screen and in the output file
		print param1, param2, omega, #"<-----Change here ----!!!"
		outputstring = str(param1)+" "+str(param2)+" "+str(omega) #"<-----Change here ----!!!"
		outputfile.write(outputstring)
		
outputfile.close()
#---------------------------------------------------------------------------
#-------------------------------------------------------------------------
