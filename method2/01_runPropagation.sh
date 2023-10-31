#! /usr/bin/bash

physiologies=(
#    'aerophilicity'
    'habitat'
#    'antimicrobial resitance'
#    'growth temperature'
#    'disease association'
#    'biosafety level'
#    'gram stain'
#    'shape'
#    'motility'
#    'arrangement'
#    'coding genes'
#    'genome size'
#    'extreme environment'
#    'spore formation'
#    'animal pathogen'
#    'plant pathogenicity'
#    'optimal ph'
#    'width'
#    'COGEM pathogenicity rating'
#    'mutation rate per site per year'
#    'antimicrobial sensitivity'
#    'pathogenicity human'
#    'hemolysis'
#    'length'
#    'biofilm forming'
#    'growth medium'
#    'mutation rate per site per generation'
#    'spore shape'
#    'health associated'
#    'sphingolipid producing'
#    'acetate producing'
#    'butyrate producing'
#    'lactate producing'
#    'hydrogen gas producing'
)

myVar=$(which R)
myRes=$(echo $myVar | grep -e "waldronlab" | wc -l)

for i in "${physiologies[@]}"
do
    echo generating data for "$i"
    if [ $myRes -eq 0 ]; then
        echo "I'm not on supermicro"
        Rscript 01_runPropagation.R "$i" "$1"
    elif [ $myRes -eq 1 ]; then
        echo "I'm on supermicro"
        /usr/bin/Rscript 01_runPropagation.R "$i" "$1"
    fi
done

