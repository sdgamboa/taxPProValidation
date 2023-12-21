#! /bin/bash

for i in $(ls train_folds | sed -e "s/\.csv$//")
do 
    echo -e "Running propagation for $i"
    hsp.py -n -t LTP_all_08_2023.newick --observed_trait_table train_folds/$i.csv -m scp -o predicted/$i.tsv --processes 1
done

