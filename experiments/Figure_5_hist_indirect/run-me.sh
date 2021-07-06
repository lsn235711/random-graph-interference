#!/bin/bash

# load the module
for i in {1..40}
do
   sbatch --export=seed=$i script.sh
done



  
