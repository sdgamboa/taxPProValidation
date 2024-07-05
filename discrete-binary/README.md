
# Discrete attributes of type "binary"

The attributes of type "binary" are those with logical
attribute values, i.e., TRUE or FALSE. An example is the "sphingolipid
producing" attribute, which could have the values
"sphingolipid producing--TRUE" or "sphingolipid producing--FALSE".

Both TRUE and FALSE values are necessary since the
ASR methods didn't perform well when only all TRUE or all FALSE
values were used.

The phytools package is used for the propagation of discrete attributes.
Attributes of type "binary" are those with logical attribute values,
i.e., TRUE and FALSE. An example is 'sphingolipid
producing'.



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
Correlation Coefficient (MCC) values.

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

Compress the \*Fold\*csv files in a .tar.gz file with date of
creation. \*csv files are ignored in the
`.gitignore` file, so the files should be
compressed in a .tar.gz file for uploading to GitHub.

```bash
DAT=$(date +%Y-%m-%d)
tar -czvf discrete-binary-$DAT.tar.gz *Fold*
```

## File column description

### Columns in the "\*test\*Fold\*.csv" files

These columns are very similar to the output of the
`bugphyzz::physiologies` function.

+ NCBI_ID. NCBI taxonomy ID with rank prefix (e.g., g__561).
+ Attribute. Physiology/Attribute name fused with the attribute value.
For example, "animal pathogen--TRUE" or "animal pathogen--FALSE".
+ Taxon name. Scientific name of the microbe.
+ Rank. From strain to Superkingdom.
+ Attribute_source. Source of the annotation.
+ Confidence_in_curation. Confidence of the source.
+ Evidence. Evidence codes from bugphyzz.
+ Frequency. Frequency codes from bugphyzz.
+ Score. Conversion of Frequency to Scores. Taken from bugphyzz.
+ Attribute_group.
+ Attribute_type. "binary".
+ taxid. NCBI taxonomy ID.

### Columns in the "\*predicted\*Fold\*csv" files

+ tip_label. Label as annotated in the phylogenetic tree used
for propagation. Some labels are NCBI_IDs with the rank prefix, e.g.,
"g__561" (Escherichia) while others are the original tips from
the LTP (Living Tree Project) tree.
+ Attribute. Physiology/Attribute name fused with the attribute value.
For example, "animal pathogen--TRUE" or "animal pathogen--FALSE".
+ Score. Posterior probability results from the ASR.
Values between 0 and 1. This will be used for MCC.
+ Accession. This is the genome accession of the species and strains
in the phylogenetic tree (at the tips). It does not apply to
the genus.
+ taxid. The NCBI taxonomy ID.
+ Taxon_name. Scientific name of the tip.
+ Rank. From strain to superkingdom.
+ \*rank\*_taxid columns. These columns contain the taxonomy ids
at different taxonomic ranks for each tip of the tree.

### Columns in \*mcc.tsv files

+ Physiology. Name of the physiology. Same as Attribute in the case
of binary attributes.
+ Rank. Taxonomic rank. all == genus, species, and strain
+ Method. ASR method (phytools) and tree (LTP).
+ Attribute. Name of the attribute.
+ Mean. Mean MCC (considering all 10 folds).
+ SD. The standard deviation of the MCC (considering all 10 folds).

### Columns in \*counts.tsv files

+ Attribute. Name of the physiology.
+ ltp_bp_phys. Number of taxa matching in both the ltp tree and bugphyzz.
+ ltp. Number of tips in the LTP tree.
+ nsti_mean. Mean nearest sequence taxonomic index (check castor definition.)
+ nsti_sd. Standard deviation of NSTI.


### Columns in the discrete_binary_summary.tsv file

The columns are the same as the ones described above. This file
is a combination of all of the \*counts.tsv and \*mcc.tsv files
into a single one.
