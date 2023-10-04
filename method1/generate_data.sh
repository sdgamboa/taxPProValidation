#! /usr/bin/bash

physiologies=(
    'aerophilicity'
    'growth temperature'
    'disease association'
    'biosafety level'
    'gram stain'
    'shape'
    'motility'
    'arrangement'
    'coding genes'
    'genome size'
    'extreme environment'
    'spore formation'
    'animal pathogen'
    'plant pathogenicity'
    'optimal ph'
    'width'
    'COGEM pathogenicity rating'
    'mutation rate per site per year'
    'antimicrobial sensitivity'
    'pathogenicity human'
    'hemolysis'
    'length'
    'biofilm forming'
    'growth medium'
    'mutation rate per site per generation'
    'spore shape'
    'health associated'
    'sphingolipid producing'
    'acetate producing'
    'butyrate producing'
    'lactate producing'
    'hydrogen gas producing'
)

for i in "${physiologies[@]}"
do
    echo "generateing data for $i"
    /usr/bin/Rscript generate_data.R "$i"
done

