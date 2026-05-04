# GENRAY_CQL3D_scripts
Set of scripts to run and analyze GENRAY/CQL3D
This framework is meant to run on engaging

If you've not run GENRAY/CQL3D before, your .bashrc will need to be modified.
You can either ask someone more knowedgeable than I exactly what's needed or just copy what's in mine (grantr/.bashrc)

Simulation runs should be stored in your scratch folder
For examples all my runs are located in /home/grantr/scratch/genray_batch/{machine}\_shots/{machine}_{shotnum}.{time}/{specific simulation run}'
If you have a different path, you will likely need to make changes in targetDirectory.py to make sure it returns all of its values properly.
At the very least, you'll need to change the username.

You'll also need to copy genr_sam.pbs and cql.pbs from /home/grantr/codes to wherever you want them to live
You'll probably need to change some permissions on those files
These are used to run the codes.

To set up a simulation run, set the targetDir as desired, then use setupInputFiles.py
If the target directory doesn't exist yet, this script will also make it.
To run the codes, cd into the target dir and either call
  sbatch genr_sam.pbs
  sbatch cql.pbs (these are copied here by setupInputFiles.py)
or call
  ./runGENThenCQL.sh . (don't forget the period!)
  
Congratulations, it should just work, but probably something will go wrong.
If nothing went wrong, then you can use the various scripts to analyze the outputs.


IMPORTANT: I have always used rho_pol = sqrt(psi_N) as my radial variable in GENRAY/CQL3D.
Many of the scripts assume this is the case, either in what they list as the x label or in the assumption that (radial variable)^2 = psi_N
If you choose a different radial variable, be sure to look through each script very carefully to see where things will need to be changed.
