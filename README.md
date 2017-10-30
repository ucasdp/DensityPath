# DensityPath
A novel algorithm, DensityPath, which accurately and efficiently reconstructs the underlying cell developmental trajectories for large-scale scRNAseq data.

Cell fates are determined by transition-states which occur during complex biological processes such as proliferation and differentiation. The advance in single-cell RNA sequencing (scRNAseq) provides the snapshots of single cell transcriptomes, offering the essential opportunity to study such complex biological processes. DensityPath, that can accurately and efficiently reconstruct the underlying cell developmental trajectories for large-scale scRNAseq data. DensityPath algorithm not only extract the separate high density representative cell states(RCSs) to handling the heterogeneous scRNAseq data accurately based on the powerful level-set cluster method, but also constructs cell state-transition path by finding the shortest paths (geodesic) of the representative cell states on the surface of density landscape.

The main steps for DensityPath are as follows:
1.Reduce the dimensionality of scRNAseq data. 
2.Estimate the density function (landscape) and the level sets.
3.Select high density clusters as the RCSs.
4.Construct the cell state-transition path.
