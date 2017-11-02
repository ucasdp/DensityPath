# DensityPath
A novel algorithm, DensityPath, which accurately and efficiently reconstructs the underlying cell developmental trajectories for large-scale scRNAseq data.

DensityPath, that can accurately and efficiently reconstruct the underlying cell developmental trajectories for large-scale scRNAseq data. DensityPath algorithm not only extract the separate high density representative cell states(RCSs) to handling the heterogeneous scRNAseq data accurately based on the powerful level-set cluster method, but also constructs cell state-transition path by finding the shortest paths (geodesic) of the representative cell states on the surface of density landscape.

The main steps for DensityPath are as follows:
1.Reduce the dimensionality of scRNAseq data. 
2.Estimate the density function (landscape) and the level sets.
3.Select high density clusters as the RCSs.
4.Construct the cell state-transition path.
5.Map the single cells onto the cell state-transition path.

The code is running as follows:
(1).Read the data files.
(2).Set the parameters of density clustering k and h.
(3).Set the path of graphic output.
(4).Select the points that needs to be mapped to the density path.
(5).Call the function: DensityPath(XX,k,h,SI_fig_dir,X), where XX represents the data in (1) ,k and h are parameters of density clustering in (2), SI_fig_dir represents the path of graphic output in (3), x represents the points that need to be mapped to the density path in (4).


The DensityPath function will output the figure of density path and mapping path, then return the list variable, "densitypath", which contains "densityKDE", "KDEdensitypeaks", "dis", "minspantreepath" and "minadjpaths". 
"DensityKDE" represents the estimated density of the sample points.
"KDEdensitypeaks" represents the two-dimensional coordinates of the density peaks.
"Dis" represents the geodesic distance between the density peaks.
"Minspantreepath" represents the paths of the minimum spanning tree on the three-dimensional density surface.
"Minadjpaths" represents the the mapping paths for the points.


Example data is a simulated dataset (testdata.csv) which contains 2480 data points on the two-dimensional space by sampling independently from 7 classes of bivariate normal distributions. The data has 3 columns, the first two columns are the two-dimensional coordinates of the points, and the third column represents the categories that the points belong to.

