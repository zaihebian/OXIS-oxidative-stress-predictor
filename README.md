# OXIS-oxidative-stress-predictor
Title: Quantitative Estimation of Intracellular Oxidative Stress in Human Tissues 

For usage:
test.ipynb is the code for estimating the oxidative stress for test samples.
inputs:
1) Gene expression matrix, such as "GSE143155_expr.csv" in the test code. You can replace this file by your own data, in which rows correspond to samples and columns correspond to gene names.
2) Sample information, such as "GSE143155_group.csv" in the test code.  You can replace this file by your own data, noting that the samples should be in the same order with the columns of the gene expression matrix.

Please put these two files in the same fold with the files set_coefs.csv and module_genes_input.csv, and make sure all these files are in the current working directory.

Then you run the code in test.ipynb, a file named "results.csv" will be created with the predicted oxidative stress in the column named OxiStress.

Plus:
OxisModel.py is the source code for training the model  
