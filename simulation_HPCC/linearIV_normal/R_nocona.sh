#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=linearIV_normal
#SBATCH --output=./report/report_%j.txt
#SBATCB –-error=./report/err_%j.txt
#SBATCH --partition nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=10:00:00


module load gcc/10.1.0
module load r/4.0.2



#Allow R to perform some automatic parallelization.
#	MKL_NUM_THREADS - The maximum number of threads you want R to spawn on your behalf.
#	$NSLOTS - This will be replaced by the number of slots you request in yout parallel environment.
#		Example:  -pe sm 36 -> $NSLOTS=36.
#export MKL_NUM_THREADS=$NSLOTS
#export MKL_NUM_THREADS=$NSLOTS
#Run the example R script using the Rscript application.
#Allow R to perform some automatic parallelization.
#	MKL_NUM_THREADS - The maximum number of threads you want R to spawn on your behalf.
#	$NSLOTS - This will be replaced by the number of slots you request in yout parallel environment.
#		Example:  -pe sm 36 -> $NSLOTS=36.
#export MKL_NUM_THREADS=$NSLOTS

#Run the example R script using the Rscript application.
#for t in 0.001 0.004 0.007  0.01 0.02 0.03 0.04 0.05 0.1 0.15 0.2 0.25 0.3
#for n in 200 500 1000 2000 5000 10000

RandomSeed=$1
#100 300 1000 
for n in 100 300 500 800 1000
do
	for s in 2 5
	do
		Rscript ./R_main.R $n $s $RandomSeed;
	done	
  
done
 
 
 


