This should be done with a modified version of the scrip hsp.py in picrust2


Run separetely due to the need for conda activation

## Create folds

```bash
physiologies=(
    'growth temperature'
    'coding genes'
    'genome size'
    'optimal ph'
    'width'
    'mutation rate per site per year'
    'length'
    'mutation rate per site per generation'
)

for i in "${physiologies[@]}"
do
	Rscript --vanilla create_folds.R "$i" &
done
```

## Run picrust2

```bash
## conda activate picrust2
./picrust2.sh
```

## Get validation metrics and plots 

```bash
physiologies=(
    'growth_temperature'
    'coding_genes'
    'genome_size'
    'optimal_ph'
    'width'
    'mutation_rate_per_site_per_year'
    'length'
    'mutation_rate_per_site_per_generation'
)

for i in "${physiologies[@]}"
do
	Rscript --vanilla validation.R "$i" &
done
```
