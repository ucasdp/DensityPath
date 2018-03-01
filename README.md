# DensityPath
The novel algorithm, DensityPath can accurately and efficiently reconstruct the underlying cell developmental trajectories for large-scale scRNAseq data. 

The code is runned as follows: 
(1) source the function: "after_pca", "densitypathimage" and "dens_bran_pseu"; 
(2) load gene expression data as OriginalData, each row of data represents one cell and each column represents one gene; 
(3) use function "after_pca(OriginalData)" to reduce dimensionality for data using PCA, and save the CSV file of the data after PCA; 
(4) load the CSV file of the data after PCA in Matlab, perform elastic embedding(EE) for reducing dimensionality to 2D and save the 2D data after EE as X; 
(5) load X in R, call the function "densitypathimage(X)" to obtain the trajectory using DensityPath algorithm; 
(6) call the function " dens_bran_pseu (X)" to get the assignment of branches and pseudotime of each cell.


Output: 
(1) the function "densitypathimage()" will output a 2×2 figure of DensityPath, which are the scatter plot of single cell points on the reduced-dimension 2-d space of gene expression after dimensionality reduction, density landscape, high density clusters and the cell state-transition path constructed by DensityPath on the 2-d density function heatmap plot, respectively;
(2) the function "dens_bran_pseu()" will return the list variable, which contains "RCSs", "pseudotime", "pseudobran", "MSTtree", "allpath". "RCSs" represents the cells contained in each high density cluster of representative cell states (RCSs). "pseudotime" represents the pesudotime calculated using DensityPath algorithm for each cell. "pseudobran" represents the cells contained in each branch assigned by DensityPath algorithm. "MSTtree" represents the minimum spanning tree of the path. "allpath" represents the labels and 2D coordinates of the peak points in each density cluster, which are on the path from the start density cluster to the end of each branch.


Example: Example data is a simulated dataset (s8971.mat, the dataset utilized here was from the material of [38]) which is charactered by 475 cells and 48 genes.
