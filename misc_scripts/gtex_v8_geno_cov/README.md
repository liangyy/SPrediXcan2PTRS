# Generate genotype covariance for GTEx V8 EN models

```
# on CRI
bash submit_gtex_v8_en_geno_cov.sh
```

```
# on washington
TISSUE=Whole_Blood
screen -dmS geno_cov_$TISSUE bash gtex_v8_en_geno_cov.screen $TISSUE
```

