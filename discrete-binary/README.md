
One-liners for running propagation:

```
## These will be executed in new tmux sessions
tmux new -d -s binary-all './01_prediction_binary.sh all'
tmux new -d -s binary-genus './01_prediction_binary.sh genus'
tmux new -d -s binary-species './01_prediction_binary.sh species'
tmux new -d -s binary-strain './01_prediction_binary.sh strain'
```



