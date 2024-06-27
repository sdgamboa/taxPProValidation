# taxPProValidation

Version: v1.0.4

Code for generating the main table:

```
/usr/bin/Rscript --vanilla create_summary_table.R
```

This must be run after running the scripts in the discrete-binary,
discrete-multistate-intersection, and numeric directories.

Output is in [validation_summary.tsv](./validation_summary.tsv)

## Version log

### v1.0.4
Update the code based on the bugphyzzExports script.
Compress "Folds" into a tar.gz file for numeric, binary, and multistate.

### v1.0.3
Remove previous attempts at performing ASR, such as using castor for
discrete attributes and using data.tree and taxonomy.

### v1.0.2
The validation was run by setting the values accepted as TRUE to 0.8 for
binary and multistate attributes. The code and summary tables were
updated accordingly.
