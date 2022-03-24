# DSNN_MatrixImputation
DSNN algorithm stands for Doubly Sparse in DCT domain with Nuclear Norm minimization. It is developed for imputing missing values in microarray gene expression data functions-
1) data_preprocessing.m - for log-transforming data if the data has high range of values.
2) data_reverseprocess.m - converting back log-transformed data to original values
3) DSNN- main algorithm for imputing missing values in gene expression 


%% NOTE: please install spgl solver and keep it in your path in matlab for this algorithm to work

Please cite the following if you use our work in your research-
Farswan A, Gupta A, Gupta R, and Kaur G, "Imputation of gene expression data in blood cancer and its significance in inferring biological pathways", Frontiers in Oncology, vol. 9, article no. 1442, pp. 1-14, January 08, 2020.
