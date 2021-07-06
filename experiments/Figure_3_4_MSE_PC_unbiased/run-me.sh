#!/bin/bash

for i in {1..20}
do
   sbatch --export=seed=$i script.sh
done



  
