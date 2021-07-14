# Random Graph Asymptotics for Treatment Effect Estimation under Network Interference
This repository contains the code to reproduce the numerical results in Li and Wager (2021).

## Folder
* experiments/: contains the code to reproduce the numerical results in our paper. Each subfolder contains the code of one figure or a group of figures. For example in folder "Figure_2_hist_direct/", the ".pdf" file is Figure 2 in the paper, and the ".R" file can be run to make the plot. 

## Replicating the experiments
* For Figure 2,6,7-9, the ".R" files in the corresponding folders can be directly executed on a laptop.
* For Figure 3,4,5, one needs to run the "run-me.sh" file on a computing cluster, and then run the "\*output\*.R" file to make the plot. The "script.sh" file is specifically designed for the [Stanford Sherlock](https://www.sherlock.stanford.edu/) cluster. One may need to modify the files correspondingly. 

## References
Shuangning Li and Stefan Wager. <b>Random Graph Asymptotics for Treatment Effect Estimation under Network Interference</b>. 2021. [[arXiv](https://arxiv.org/abs/2007.13302)]
