# Spatial Perturb-Seq

This repository contains a vignette and R code showing some functions for a single-cell Spatial Perturb-Seq dataset, consisting of cell x (genes + barcodes + XY) matrix.

## Vignette

[View the HTML Vignette](https://kimberle9.github.io/spatialperturbseq/vignette.html)

## Description

This project aims to add to Seurat's [1] functionality , for the analysis and visualization of single-cell perturb-seq data. These functions manipulate Seurat objects and require Seurat v5. The main components are an R Markdown file, an HTML file, and an R script. The functions are described as follows.

## Installation

You can install the latest version of **spatialperturbseq** from GitHub:
```R
# Install devtools if not already installed
install.packages("devtools")

# Install spatialperturbseq
devtools::install_github("kimberle9/spatialperturbseq")
```

## Features

#### 1. add_image_slice

This function is for Stereo-seq and Stereopy-generated Seurat objects, which do not come with an image component by default. This function adds an image slice to the Seurat object, allowing it to be compatible with Seurat [1], BANKSY [2], and other packages.

#### 2. annotate_pos_neg*

Given a Seurat object and a list of features (eg. gRNA barcodes), creates a metadata column in the Seurat object stating "Pos" or "Neg" for each cell, ie classifies cells containing or not containing gRNA/barcodes. 


#### 3. barplot_pos*

After Pos/Neg is annotated, we want to see if positive cells are enriched in any identity classes of cell types. This function takes an ident input and plots the percentage of barcode-positive cells for each identity class.


#### 4. pos_neg_DEG*

After Pos/Neg is annotated, calculates DEGs and shows top upregulated genes for each population. 

#### 5. spatial_pos_neg

Plots all cells on XY axes, with negative cells in black and positive cells in red.

#### 6. nndist_pos_neg

Nearest neighbour distances aids in assessing clustering or spatial distribution of positive cells. This function plots a frequency histogram of Nearest Neighbor distances for Pos points in a spatial region, highlighting distribution of distances. A left-skewed graph suggests clustering of Pos points while a right-skewed graph suggests a more dispersed pattern. A uniform graph could indicate random distribution without clustering or dispersion, while multimodal peaks could indicate multiple clusters. 

#### 7. annotate_gRNA*

Given a Seurat object and a list of features (eg. gRNA barcodes), creates a metadata column in the Seurat object stating the gRNA/barcode expressed in each cell. If the cell contains no barcodes, it will be labelled "Neg". If the cell contains multiple barcodes, it will be labelled "Multiple". 

#### 8. piechart_barcode_uniqueness*

Ideally, most gRNA/barcode-containing cells contain only 1 unique gRNA/barcode, unless a combination of multiple perturbations is desired. Given a Seurat object and a list of features (eg. gRNA barcodes), creates a piechart showing the proportion of cells containing 1, 2, 3 etc barcodes. Remove0 is an optional parameter (default TRUE) that allows inclusion of Neg cells when set to FALSE. This plot gives a sense of how spread out the perturbations are among the cells. This does not require the Seurat object to have spatial coordinates.

#### 9. barplot_barcode_count*

Bar plot showing the number of cells for each perturbation, after gRNA is annotated.

#### 10. plot_perturbation_effect*

To visualize impact of each perturbation on gene expression, this function creates a bar plot of number of DEGs and average log2FC of those DEGs. Input parameters include the list of perturbations and the control (ideally perturbation in a safe harbour region like "sgrna-msafe"), logfc.threshold (default 0.5), test.use (default = "wilcox"). 

#### 11. dotplot_upregulated*

Dotplot of upregulated genes for each perturbation. Usually, one of the upregulated genes would be the barcode itself, therefore the parameter include_sgrna is FALSE by default. Set include_sgrna to TRUE for a sanity check. 

#### 12. plot_barcode_spread

Given a Seurat object with spatial coordinates, and annotateWhichgRNA function ran, plots the spatial distribution of the pertubations. Default pt.size is 1.  

#### 13. get_euclidean_neighbours

Given a Seurat object with spatial coordinates, and a list of CellIDs, uses Rccp [2] to identify the nearest neighbors of specific cells based on euclidean distance. The k_geom parameter (default = 15) is the number of neighbours to compute per cell. For example, if 2 cells are given, and the k_geom parameter is 15, the maximum number of neighbours returned in 30 (2 x 15). However, there could be less than 30 cells returned as the 2 cells might share neighbours. This allows analysis (eg. DEGs) to be run on a perturbed cell's neighbours, allowing the study of non-cell autonomous function of genes. This requires the Seurat object to have spatial coordinates.

#### 14. get_KNN

Given a Seurat object with spatial coordinates, and a list of CellIDs, uses BANKSY's [2] computeNeighbors function to return a list of CellIDs of nearest neighbours of input cells. Similar to the get_euclidean_neighbours function, k_geom parameter is 15 by default. 


#### 15. get_radial_neighbours

This function is used to identify neighboring cells that are within a specified radial distance from a given set of target cells in a spatial scRNA-seq dataset. It calculates the Euclidean distances between the spatial coordinates of each target cell in cell_list and all other cells in the dataset. The function then filters these distances to include only those cells that are within the specified radius (default 200 units), excluding the target cell itself. Each target cell could have a different number of neighbours, as the number of neighbours is not fixed. 


#### 16. positive_neighbor_histogram

Given a Seurat object with spatial coordinates, and list of cells (ie gRNA/barcode-positive cells), plots a frequency histogram showing the frequency of cells having 0, 1, 2 etc positive neighbours. The k_geom parameter is the number of neighbours to compute per cell. The smaller the k_geom, the more likely a cell is to have 0 positive neighbours. This function visualizes the spread of perturbations spatially. Ideally, gRNA/barcode positive cells do not have gRNA/barcode positive neighbours, to ensure that the effects of perturbations are not confounded by interactions between adjacent cells. In other words, most cells should have 0 positive neighbours. This requires the Seurat object to have spatial coordinates.


#### 17. module_score*

Given a Seurat object with a module score computed with Seurat's AddModuleScore function, the name of the module, and a list of Idents to plot, plots a VlnPlot showing significance between all the combinations of pairs. Eg. If 2 idents are given, computes significance once. If 3 idents are given, computes significance 3 times. If 4 idents are given, computes significance 6 times. Increment parameter (default 0.25) adjusts the distance between the significance bar. This does not require the Seurat object to have spatial coordinates.


#### 18. liana_table* 

This interrogates the LIANA database [4] to identify ligand-receptor interactions in a Seurat object. The output is a table of ligand, receptors, p-values and scores. This does not require spatial coordinates and can be used with any single-cell perturb-seq data. With spatial transcriptomics data, this function can be used in conjunction with getBanksyNeighbours function to identify LR interactions between a cell and its local neighbours, as cell-cell communication usually occurs between cells and nearby neighbours. 


#### 19. plot_ligand_receptor 

This function aids visualization of ligand and receptor interactions between cells and their local neighbours. It returns a plot of individual cells (source cells) and their neighbours (target cells) on a spatial map, with expression of ligand in red (source cells) and receptor in blue (target cells). Together with the seurat object with spatial coordinates, the list of source cells and target cells must be provided, as well as a ligand and receptor. 


#### 20. power_estimation* 

This function takes in a seurat object, and the cells for comparison (perturbation and control), and estimates the power for detecting differentially expressed (DE) genes between the two groups of cells using the scPower package [5]. Priors are calculated based on estimation of negative binomial parameters for each gene using the DESeq2 package, followed by fitting a gamma distribution on the dataset, to capture mean and dispersion. This is used to quantify expression probability of genes. This function returns a detection power estimate (number between 0-1). Overall detection power = expression probability of a gene x DE power of a gene. A power of 0.5 means there is a 50% chance of detecting a true effect or difference if it exists. A power estimate of >0.8 is generally accepted. Values <0.8 suggest the need to increase the sample size and number of cells, given the parameters of the dataset.  

* Does not require the Seurat object to have spatial coordinates, and can be used for non-spatial perturb-seq experiments. 

## Citation

If you use **spatialperturbseq** in your research, please cite it as: Shen, K., Seow, W. Y., Keng, C. T., Shern, D. L., Guo, K., Meliani, A., … Chew, W. L. (2024). Spatial Perturb-Seq: Single-cell functional genomics within intact tissue architecture. doi:10.1101/2024.12.19.628843


## Files

- `vignette.Rmd`: The R Markdown file containing the vignette.
- `vignette.html`: The HTML output of the vignette.
- `vignette.R`: The R script with the code used in the vignette.
- `README.md`: This README file.


## References

1. Stuart, T., Butler, A., Hoffman, P., Hafemeister, C., Papalexi, E., Mauck III, W. M., ... & Satija, R. (2019). Comprehensive integration of single-cell data. Cell, 177(7), 1888-1902.e21. https://doi.org/10.1016/j.cell.2019.05.031
2. Singhal, V., Chou, N., Lee, J., Yue, Y., Liu, J., Chock, W. K., Lin, L., Chang, Y. C., Teo, E. M. L., Aow, J., Lee, H. K., Chen, K. H., & Prabhakar, S. (2024). BANKSY unifies cell typing and tissue domain segmentation for scalable spatial omics data analysis. Nature genetics, 56(3), 431–441. https://doi.org/10.1038/s41588-024-01664-3
3. Eddelbuettel D, François R (2011). “Rcpp: Seamless R and C++ Integration.” Journal of Statistical Software, 40(8), 1–18. doi:10.18637/jss.v040.i08
4. Dimitrov, D., Türei, D., Garrido-Rodriguez M., Burmedi P.L., Nagai, J.S., Boys, C., Flores, R.O.R., Kim, H., Szalai, B., Costa, I.G., Valdeolivas, A., Dugourd, A. and Saez-Rodriguez, J. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022). https://doi.org/10.1038/s41467-022-30755-0
5. Schmid, K. T., Höllbacher, B., Cruceanu, C., Böttcher, A., Lickert, H., Binder, E. B., Theis, F. J., & Heinig, M. (2021). scPower accelerates and optimizes the design of multi-sample single cell transcriptomic studies. Nature communications, 12(1), 6625. https://doi.org/10.1038/s41467-021-26779-7

