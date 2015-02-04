import os
import sys
import glob
import shutil

#Root path
rp = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]
sys.path.append(rp)

#-----------------------------------------------------------------------#
def initialize_MadDM_session(print_banner = False):
#-----------------------------------------------------------------------#
#                                                                       #
#  This routine starts the whole program.  The desired model and        #
#  project name are input by the user.  The model is checked against    #
#  the model directory in the main madgraph folder.                     #
#                                                                       #
#-----------------------------------------------------------------------#

  if (print_banner):
    print "\n"
    print \
  "            ====================================================\n"+\
  "            |                    MadDM v1.0                     |\n"+\
  "            ====================================================\n"+\
  "                                                                               \n"+\
  "                #########                                                      \n"+\
  "             ###\\\\####//#####                                                \n"+\
  "           ######\\\\##//########                                              \n"+\
  "          ########\\\\//###########                                            \n"+\
  "         #########//\\\\############                                           \n"+\
  "        #########//##\\\\############    Created at the University of Kansas   \n"+\
  "       ########//#####\\\\###########                                          \n"+\
  "       ######################### ## ___________________________________________\n"+\
  "       ####################### 0  #  _     _               _  _____   _     _  \n"+\
  "       #############   0  ###    ## | \   / |   ___    ___|| | ___ \ | \   / | \n"+\
  "       ##############    #########  ||\\\\ //|| / __ |  / __ | ||   || ||\\\\ //|| \n"+\
  "        ##########################  ||  V  || ||__||  ||__|| ||___|| ||  V  || \n"+\
  "         ###################   ##   ||     || \_____\ \____| |_____/ ||     || \n"+\
  "          ############       ###    ___________________________________________\n"+\
  "           ##########    ######                                                 \n"+\
  "             ################                                                   \n"+\
  "                 ########                                                       \n"+\
  "                                                                                    \n" 

  # Get the model name and check if it exists in the Models folder
  model_list = os.listdir(os.path.join(rp, 'models'))

  # Loop until the model name is ok
  leg_ans = False
  while not leg_ans:
        	
	  model_name = raw_input('Enter model name (press Enter to list the available models): ')
	  if (model_name == ''):
		  # List the available DM models in the Models directory of MadGraph

		  for mdl in model_list:
			  mdl_dir = os.path.join('models', mdl)
			  if os.path.isdir(os.path.join(rp, mdl_dir)):
				  print mdl+"   "
		
	  elif model_name in model_list:
			leg_ans = True
	  else:
		  print "The model you entered is not available! Please try again"

  # Loop until the project name is set
  leg_ans = False
  while not leg_ans:
	
    project_name = raw_input("Enter project name ("+model_name+"): ")
    if (project_name == ''):
      project_name = model_name

    project_list = glob.glob('Projects/'+project_name+'_*')
    if (project_list == []):
      leg_ans = True
      new_project = True
    else:
      answer = raw_input("Project directory already exists. Overwrite?[n] (y/n):")
      if ((answer == 'y') or (answer == 'Y')):
        leg_ans = True
        new_project = True
        for project in project_list:
          shutil.rmtree(project)
      if ((answer == 'n') or (answer == 'N') or (answer == '')):
        leg_ans = True
        new_project = False

  names = [new_project, model_name,project_name]
  return names
