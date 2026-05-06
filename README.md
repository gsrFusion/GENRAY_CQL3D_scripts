# GENRAY_CQL3D_scripts
Set of scripts to run and analyze GENRAY/CQL3D
This framework is meant to run on engaging

If you've not run GENRAY/CQL3D before, there will likely need to be some module setup required to be able to run GENRAY/CQL3D.
I suggest asking your favorite neighborhood computationalist

These scripts require that simulations runs are stored in a folder structure of 
  {machine}\_shots/{machine}\_{shotNum}.{time}/{machine}\_{shotNum}.{time}\_descriptionOfParticularRun
I suggest locating these folders in scratch.
For examples all my runs are located in /home/grantr/scratch/genray_batch/{machine}\_shots/{machine}_{shotnum}.{time}/{specific simulation run}

To choose which simulation the scripts are looking at, modify targetDirectory.py 

I suggest copying the example case in testCase/DIIID_DIIID_180403.04400_n2.8Npara_1MW_frameworkTest to
  DIIID_shots/DIIID_180403.04400/DIIID_180403.04400_n2.8Npara_1MW_frameworkTest and trying to run it with sbatch genr_sam.pbs, then sbatch cql.pbs
  You will need to modify these .pbs files to point to your GENRAY and CQL3D executables.

From there, hopefully the scripts just work. If they don't you can send me an email at grantr@psfc.mit.edu

See the howToAddaNewShot.txt file for a simple guide on how to set up a new case

IMPORTANT: I have always used rho_pol = sqrt(psi_N) as my radial variable in GENRAY/CQL3D.
Many of the scripts assume this is the case, either in what they list as the x label or in the assumption that (radial variable)^2 = psi_N
If you choose a different radial variable, be sure to look through each script very carefully to see where things will need to be changed.
