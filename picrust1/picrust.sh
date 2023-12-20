#! /bin/bash
# 
# myVar=$(which R)
# myRes=$(echo $myVar | grep -e "waldronlab" | wc -l)

# echo -e "$myVar"

# echo -e ">>>>> Creating folds..."
# if [ $myRes -eq 0 ]; then
#     echo ">>>>> I'm not on supermicro"
#     Rscript --vanilla create_folds.R 
# elif [ $myRes -eq 1 ]; then
#     echo ">>>>> I'm on supermicro"
#     /usr/bin/Rscript --vanilla create_folds.R
# fi
# 
# 
# Rscript --vanilla create_folds.R

## Step 1 - Get prunde and reference trees, and trait table in the right format.
for i in $(ls -1 train_folds | sed -e 's/\.csv$//')
do
    echo -e "Formatting input files for $i"
    formattedDirName=formatted/$i
    format_tree_and_trait_table.py --input_table_delimiter=tab -t LTP_all_08_2023.newick -i train_folds/$i.csv -o $formattedDirName
    sed -i -e "s/'//g" $formattedDirName/reference_tree.newick
    sed -i -e "s/'//g" $formattedDirName/pruned_tree.newick

    echo -e "Running ASR for $i"
    asrDirName=asr/$i
    ancestral_state_reconstruction.py -i $formattedDirName/trait_table.tab -t $formattedDirName/pruned_tree.newick -o $asrDirName/asr_counts.tab -c $asrDirName/asr_ci.tab

    echo -e "Predicting traits for $i"
    predDirName=predicted/$i
    predict_traits.py -a -i $formattedDirName/trait_table.tab -t $formattedDirName/reference_tree.newick -r $asrDirName/asr_counts.tab -c $asrDirName/asr_ci.tab -o $predDirName/predict_traits.tab
done

## Perform ancestral state reconstruction 


#ancestral_state_reconstruction.py -i formatted/trait_table.tab -t formatted/pruned_tree.newick -o asr_counts.tab -c asr_ci.tab

## Create predictions
#predict_traits.py -a -i formatted/trait_table.tab -t formatted/reference_tree.newick -r asr_counts.tab -c asr_ci.tab -o predict_traits.tab
