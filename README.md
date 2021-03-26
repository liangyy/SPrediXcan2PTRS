# SPrediXcan2PTRS

In this repository, we provide a light-weight script to convert S-PrediXcan results to PTRS weights.

# Input

1. S-PrediXcan results.
2. Genotype reference panel (match the GWAS/prediction model population!)
3. Hyperparameters (more details in below): 
    - `alpha`: The relative amount of l1 and l2 penalty. 
    - `lambda` sequence: `lambda` controls the amount of weights to put on the regularization terms. No need to specify the whole sequence but need to provide `nlambda` (the length of `lambda` sequence), the script will determine `lambda_max`. 
    - `offset` sequence: `offset` controls the amount of weights to be put on the diagonal of predicted expression correlation matrix. Need to provide the whole `offset` sequence, typically, use something like `[0.01, 0.1, 0.2, 0.3]`. 
    
# Output

A gene x lambda x offset matrix containing the gene weights at each setting.

# Model details

TODO
