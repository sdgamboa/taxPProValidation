
One-liners for running propagation (in new tmux sessions)

```
tmux new -d -s multistate-all './01_prediction_multistate.sh all'
tmux new -d -s multistate-genus './01_prediction_multistate.sh genus'
tmux new -d -s multistate-species './01_prediction_multistate.sh species'
tmux new -d -s multistate-strain './01_prediction_multistate.sh strain'
```
Run validation:

```
./02_validation_multistate.sh phytools-ltp
```

Create summary:

```
/usr/bin/Rscript --vanilla 03_summary_multistate.R
```

Output is in the [discrete_multistate_summary.tsv](./discrete_multistate_summary.tsv) file.

```
DAT=$(date +%Y-%m-%d)
tar -czvf binary-multistate-$DAT.tar.gz *Fold*
```
