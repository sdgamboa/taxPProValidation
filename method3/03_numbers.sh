#! /bin/bash

## this is meant to be run in method2/ from the terminal
physiologies=($(ls -l *csv | sed -e "s/.* \(.*\)_\(test\|propagated\)_\(all\|species\|strain\|genus\|\)_.*$/\\1/" | sort | uniq))

myVar=$(which R)
myRes=$(echo $myVar | grep -e "waldronlab" | wc -l)

for i in "${physiologies[@]}"
do
    echo generating count summary for "$i"
    if [ $myRes -eq 0 ]; then
        echo "I'm not on supermicro"
        Rscript 03_numbers.R "$i"
    elif [ $myRes -eq 1 ]; then
        echo "I'm on supermicro"
        /usr/bin/Rscript 03_numbers.R "$i"
    fi
done
