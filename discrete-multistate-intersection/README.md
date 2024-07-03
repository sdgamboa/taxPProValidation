
# Discrete attributes of type "multistate-intersection"

The attributes of type "multistate-intersection" are those with several
mutually exclusive attribute values of class "character". An example,
would be the attribute "aerophilicity" with attribute values "aerobic",
"anaerobic", and "facultatively anaerobic".

In some cases, the same taxon could have two mutually exclusive annotations.
This could happen in a few situtations. For example, a genus could contain
species with "anaerobic" and "facultive anaerobic" microbes. This could
also happen if a taxon is annotated differently in two sources and these
sources have the same "confidence in curation" value. In those cases,
the "Score" is divided among the annotations. For example, 
"aerobic" with score 0.7, "anaerobic" with score 0.1, and "facultatively
anerobic" with score 0.2.




Still, this cases
are rare.

The phytools package was used for the propagation of discrete attributes,
includiong those of type "multistate-intersection".

The scripts here:

1. Create files with the propagation results in a 10-fold
cross-validation approach (folds 1 through 10). The "test" files contain
the real annotations that were left out of the propagation. The
"propagation" files contain the results of the propagation without the
test set.

2. Calculate the Matthew's correlation coefficient using the \*Fold\*csv files.

3. Create a summary with MCC values at different ranks and other metrics.

## Steps for reproducing

### 1. Modify the 01_prediction_binary.sh file

Choose physiologies from the `bugphyzz::showPhys`
function that should be cross-validated. Only add physiologies with attribute type "binary."

Add the chosen physiologies to the `physiologies`
variable in the script. No spaces are allowed. Use
underscore (_) instead of spaces, if any.

### 2. Run the sh scripts

It is recommended to run this interactively using tmux and check the
intermediate results.

Run the following one-liners in a linux terminal with bash:

```bash
## These will be executed in new tmx sessions
tmux new -d -s multistate-all './01_prediction_multistate.sh all'
tmux new -d -s multistate-genus './01_prediction_multistate.sh genus'
tmux new -d -s multistate-species './01_prediction_multistate.sh species'
tmux new -d -s multistate-strain './01_prediction_multistate.sh strain'
```

### 3. Check output

Verify the results. Check the *.log and *.msg files
for errors that could have interrupted the process.
Some \*Fold\*.csv files should have been created.

### 4. Run the validation

Run validation script in bash:

```bash
./02_validation_multistate.sh phytools-ltp
```
After running this script, some *.tsv files must have
been created. These contain the results of the
10-fold cross-validation with the Matthew's
Correlation Coefficient (MCC) values.


### 5. Create summary

Run the following one-liner in bash:

```bash
Rscript --vanilla 03_summary_multistate.R
/usr/bin/Rscript --vanilla 03_summary_multistate.R
```
The result is in the [discrete_multistate_summary.tsv](./discrete_multistate_summary.tsv) file.

This file is a summary of all of the MCC results
generated above.

### 6. Compress the *csv files in a tar.gz file

Compress the \*Fold\*csv files in a .tar.gz file with date of
creation. \*csv files are ignored in the
`.gitignore` file, so the files should be
compressed in a .tar.gz file for uploading to GitHub.

```
DAT=$(date +%Y-%m-%d)
tar -czvf discrete-multistate-$DAT.tar.gz *Fold*
```













