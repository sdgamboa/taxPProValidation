
# Discrete attributtes of type "binary"

The phytools package is used for the propagation of discrete attributes.
Attributes of type "binary" are those with logical attribute values,
i.e., TRUE and FALSE. An example is 'sphingolipid
producing'.  Both TRUE and FALSE values are necessary since the
ASR didn't perform well when only TRUE or only FALSE
values were used.

The scripts here:

1. Create files with the propagation results in a 10-fold
cross-validation approach (folds 1 through 10). The "test" files contain
the real annotations that were left out of the propgation. The
"propagation" files contain the results of the propagation withouht the
test set.

2. Calculate the Matthew's correlelation coefficient using the \*Fold\*csv files.

3. Create a summary with MCC values at different ranks and other metrics.

## Steps for reproducing

### 1. Modify the 01_prediction_binary.sh file

Choose physiologies from the `bugphyzz::showPhys`
function that should be cross-validated. Only add physiologies with attribute type "binary."

Add the chosen physiologies to the `physiologies`
variable in the script. No spaces are allowed. Use
underscore (_) instead of spaces if any.

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

After running this script, some *.tsv files must have
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

Compress the *Fold*csv files in a .tar.gz file with date of
creation. *csv files are ignored in the
`.gitignore` file so the files should be
compressed in a .tar.gz file for uploading to GitHub.

```bash
DAT=$(date +%Y-%m-%d)
tar -czvf discrete-binary-$DAT.tar.gz *Fold*
```

## Columns