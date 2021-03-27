# SPrediXcan2PTRS

In this repository, we provide a light-weight script to convert S-PrediXcan results to PTRS weights. 
The idea behind the theme is similar to `lassosum` in the context of GWAS: [Mak et al (2017) Polygenic scores via penalized regression on summary statistics. Genetic Epidemiology 41(6) 469-480](https://onlinelibrary.wiley.com/doi/full/10.1002/gepi.22050)

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

We fit the following model:

```
argmin_x || y - G x ||_2^2 + lambda_0 * ( alpha * ||x||_1 + (1 - alpha) / 2 ||x||_2^2 ) 
```
, where `y` is the phenotpe and `G` is the predicted expression matrix. And we assume that `y` and `G` are standardized.

This problem is equivalent to:

```
argmin_x x' R x - 2 b' x + lambda * ( alpha * ||x||_1 + (1 - alpha) / 2 ||x||_2^2 ) 
```
, where `R` is the predicted expression correlation matrix, `b = z_predixcan / sqrt(N)` with `N` being the GWAS sample size. Note that `lambda =/= lambda_0` and they are matched up to a factor: `lambda = lambda_0 / N`. 

In practice, we need to calculate an approximation of `R` from a reference panel. To stablize the model fitting, we add an `offset` term (see more details in `lassosum` paper).
Eseentially, we calculate `R_tilde` from a LD reference panel and then we use `(1 - offset) R_tilde + offset I` to replace `R`. 

# Dependencies

Use `conda env create -f environment.yml` to install the dependencies. 
Also, need to pre-append `transethnic_prs` [link](https://github.com/liangyy/transethnic_prs) to `PYTHONPATH`.
