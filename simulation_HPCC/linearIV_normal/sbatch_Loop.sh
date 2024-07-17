# !/bin/bash

for Index_loop in {1..500};
do
	sbatch R_nocona.sh $Index_loop;
	echo 'Submitted times = ' $Index_loop;
done