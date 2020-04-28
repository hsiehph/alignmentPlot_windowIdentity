# alignmentPlot_windowIdentity
The pipeline consists of three steps:
  1) using MAFFT to produce pairwise sequence alignments;
  2) running a custom PYTHON script to smooth identity along the alignment using three choices of windows/steps;
  3) plotting the results using R (ggplot2)
  
Note that 
1) it is advised to have input sequences with reasonable lengths (<1Mb ideally) so the alignment can be effectively done;
2) the MAFFT command in the Snakemake file may need to be modified to get a desired alignemnt.
