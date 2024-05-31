# taxPProValidation

Version: v1.0.2

Code for generating the main table:

```
/usr/bin/Rscript --vanilla create_summary_table.R
```

This must be run after running the scripts in the discrete-binary,
discrete-multistate-intersection, and numeric directories.

Output is in [validation_summary.tsv](./validation_summary.tsv)


## Version log

### v1.0.2
The validation was run by setting the values accepted as TRUE to 0.8 for
binary and multistate attributes. The code and summary tables were
updated accordingly.
