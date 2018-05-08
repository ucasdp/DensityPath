# DensityPath
The novel algorithm, DensityPath can accurately and efficiently reconstruct the underlying cell developmental trajectories for large-scale scRNAseq data. 

The code is runned as follows: 
1.	The main part of DensityPath is performed in R and the dependent packages has been listed in DensityPath.R script. All the dependencies will be automatically installed using the commands in the set_packages.R script.
2.	The first dimension reduction part through PCA can be performed in after_pca.R file.
3.  The code of dimension reduction part elastic embedding is available in the website of  Miguel  A ́.Carreira-Perpin ̃a ́n and can be downloaded in http://faculty.ucmerced.edu/mcarreira-perpinan/research/software/code-EE_SNE_tSNE.tar.gz, the usage of EE can refer to the examples. In our study, the script demo_COIL.m is adjusted to perform the dimension reduction as follows:
  A.	The variable Y is the processed data after PCA via the after_pca() function and the useless variable Y0, M, t0 and t can be removed. Each row of Y represents a cell and each column means one gene.
  B.	The variable Y needs to be normalized as Y0, which means:
                                               Y = Y-min(Y(:));
                                               Y = Y./max(Y(:));
  C.	The calculation and normalization of the weights Wp and Wn are as in the code and the parameters do not change.
  D.	To achieve the full convergence, we cancel the limitation of the running time and maximum iteration, which means we deleted the variable opts.runtime and opts.maxit.
  E.	The variable l in the code is the parameter lambda in the objective function of elastic embedding. It was set to be 10 in our study while it’s user-specified parameter actually.
4.	The data after dimension reduction elastic embedding now is the coordination of cells in 2-d space and it can be applied in the following DensityPath analysis. 
  A.	The function densitypathimage() will output a 2×2 figure of DensityPath, which are the scatter plot of single cell points on the reduced-dimension 2-d space of gene expression after dimensionality reduction, density landscape, high density clusters and the cell state-transition path constructed by DensityPath on the 2-d density function heatmap plot, respectively.
  B.	The function dens_bran_pseu() can figure out the cells in representative cell states (RCSs), the pseudo-branch of the points’ structure. With the prior knowledge of start cell, the function can still output the pseudo-time of each cell along the trajectory, the paths from the start RCS which is closest to the start cell to the end RCSs, and pseudo-order of each path from start RCS to one end RCS.
  
