This should be done with a modified version of the scrip hsp.py in picrust2


Run separetely due to the need for conda activation

```
Rscript --vanilla create_folds.R growth_temperature
## ./piscrust.sh needs conda with picrust2
./picrust2.sh # this will just run for any physiology (no need for pairwise processing)
Rscript --vanilla validation.R growth_temperature
```
