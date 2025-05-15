#!/bin/bash

#####################################################
##        Important parameters of your job         ##
##             are specified here 		   ##
#####################################################

#SBATCH --time=4:00:00				## total computing time
#SBATCH --nodes=1				## number of nodes 
#SBATCH --ntasks-per-node=1			## number of tasks per node
#SBATCH --cpus-per-task=32			## number of CPUs per task
#SBATCH --mem=128GB				## memory per node
#SBATCH --partition=secondary			## queue
#SBATCH --output=di_runtime_test_cen2.out		## file that will receive output from execution
#SBATCH --error=di_runtime_test_cen2.err		## file that will receive any error messages
#SBATCH --job-name=csr_runtime_test		## job name
#SBATCH --mail-user=vikramr2@illinois.edu
#SBATCH --mail-type=BEGIN,END

########## Run your executable ######################

cd ../build
cmake .. && make
time ./reccs ../data/cen.tsv -v -c ../data/cen_0.01.tsv
