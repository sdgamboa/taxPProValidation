
One-liners for running propagation:

```
## These will be executed in new tmux sessions
tmux new -d -s binary-all './01_prediction_binary.sh all'
tmux new -d -s binary-genus './01_prediction_binary.sh genus'
tmux new -d -s binary-species './01_prediction_binary.sh species'
tmux new -d -s binary-strain './01_prediction_binary.sh strain'
```

Run validation script

```
./02_validation_binary.sh phytools-ltp
```

Create summary:

```
/usr/bin/Rscript --vanilla 03_summary_binary.R 
```

Output is in the `discrete_binary_summary.tsv` file.
