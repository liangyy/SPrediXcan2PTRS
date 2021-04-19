# Generate genotype covariance for Gtex V8 EN models

```
# on CRI
bash submit_gtex_v8_en_geno_cov.sh
```

```
# on washington
TISSUE=Whole_Blood
screen -dmS bash gtex_v8_en_geno_cov.screen $TISSUE
```