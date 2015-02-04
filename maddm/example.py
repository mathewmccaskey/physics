#! /usr/bin/env python
from init import *
from darkmatter import *

#Create the relic density object. 
dm=darkmatter()
#Initialize it from the rsxSM model in the MadGraph model folder, 
#and store all the results in the Projects/rsxSM subfolder. 
dm.init_from_model('rsxSM', 'rsxSM', new_proj = True)

# Determine the dark matter candidate...
dm.FindDMCandidate('x1')

#...and all the coannihilation partners.
dm.FindCoannParticles('x2')

#Get the project name with the set of DM particles and see 
#if it already exists.
dm.GetProjectName()

#Generate all 2-2 diagrams.        
dm.GenerateDiagrams()    

#Print some dark matter properties in the mean time.
print "------ Testing the darkmatter object properties ------"
print "DM name: "+dm._dm_particles[0].get('name')
print "DM spin: "+str(dm._dm_particles[0].get('spin'))
print "DM mass var: "+dm._dm_particles[0].get('mass')
print "Mass: "+ str(dm.GetMass(dm._dm_particles[0].get('pdg_code')))+"\n"
print "Project: "+dm._projectname

#Output the FORTRAN version of the matrix elements 
#and compile the numerical code.
dm.CreateNumericalSession()

#Calculate relic density.
omega = dm.CalculateRelicAbundance()
print "----------------------------------------------"
print "Relic Density: "+str(omega),
print "----------------------------------------------"


