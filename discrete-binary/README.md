
# Discrete binary attributes

This script runs the code for running the 10-fold
cross-validation of the attributes/physiologies with
attribute type "binary". An example is 'sphingolipid
producing' and it must have values TRUE and FALSE. 
Both TRUE and FALSE values are necessary since the
ASR didn't perform well when only TRUE or only FALSE
values were used.


## Steps for reproducing

### 1. Modify the 01_prediction_binary.sh file

Choose physiologies from the `bugphyzz::showPhys`
function that should be cross-validated. Only add physiologies with attribute type "binary."

Add the chosen physiologies to the `physiologies`
variable in the script. No spaces are allowed, use
underscore (_) instead if the name of the physiology
contains any.

### 2. Run the sh scripts

It is recommended to run this interactively using tmux
and check the intermediate results.

Run the following one-liners in a linux terminal with
bash:

```bash
## These will be executed in new tmux sessions
tmux new -d -s binary-all './01_prediction_binary.sh all'
tmux new -d -s binary-genus './01_prediction_binary.sh genus'
tmux new -d -s binary-species './01_prediction_binary.sh species'
tmux new -d -s binary-strain './01_prediction_binary.sh strain'
```
### 3. Check output

Verify the results. Check the *.log and *.msg files
for errors that could have interrupted the process.
Some \*Fold\*.csv files should have been created.

### 4. Run the validation

Run validation script in bash:

```bash
./02_validation_binary.sh phytools-ltp
```

After running this script, some *.csv files must have
been created. These contain the results of the
10-fold cross-validation with the Matthew's
Correlation Coefficient (MCC).

### 5. Create summary

Run the following one-liner in bash:

```bash
## On supermicro
/usr/bin/Rscript --vanilla 03_summary_binary.R 

## On a machine with R installed in PATH
Rscript --vanilla 03_summary_binary.R 
```
The result is in the [discrete_binary_summary.tsv](./discrete_binary_summary.tsv).

This file is a summary of all of the MCC results
generated above.

### 6. Compress the *csv files in a tar.gz file

Compress the files in a tar.gz file with date of
creation. *csv files are ignore in the
`.gitignore` file so the files should be
compressed in a tar.gz file for uploading to GitHub.

```
DAT=$(date +%Y-%m-%d)
tar -czvf discrete-binary-$DAT.tar.gz *Fold*
```