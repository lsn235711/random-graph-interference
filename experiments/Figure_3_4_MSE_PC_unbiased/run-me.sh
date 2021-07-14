#!/bin/bash

for i in {1..200}
do
   sbatch --export=seed=$i script.sh
done



  
