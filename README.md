How Well Can Multivariate and Univariate GWAS Distinguish Between True and Spurious Pleiotropy?
================

Fernandes, Samuel B.; Zhang, Kevin, S.; Jamann, Tiffany M.; Lipka, Alexander E. Frontiers in Genetics (2021). https://www.frontiersin.org/article/10.3389/fgene.2020.602526     

The following files contain:

pipeline_ST.R: Pipeline to simulate and analyze ST scenarios  
pipeline_LD.R: Pipeline to simulate and analyze LD (Spurious Pleiotropy) scenarios  
pipeline_P.R: Pipeline to simulate and analyze P scenarios  

run_simulation_ST.R, run_simulation_LD.R, run_simulation_P.R: Functions to call simplePHENOTYPES in each genetic architecture  

gemma.R: a wrapper to run the software GEMMA from R.  
run_gemma_ST, run_gemma_LD, run_gemma_P: Function to call gemma.R in each genetic architecture  

graphs_combined.R: scripts used to generate all figures. Figures for different distances or FDR were generated with the same script by changing the file with the dataset.  

data: Marker data and kinship information for each dataset used in the simulation  

It is necessary to create a folder named 'simulation' containing subfolders named 'ST', 'LD', and 'P' to run the pipeline script.  

Contact: <samuelf@illinois.edu> or <fernandessb101@gmail.com>
