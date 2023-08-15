#! /usr/bin/bash

physiologies=(
    'acetate producing'
    'aerophilicity'
    'growth temperature'
)

for i in "${physiologies[@]}"
do
    echo "$i"
    /usr/bin/Rscript generate_data.R "$i"
done
