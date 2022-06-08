#!/bin/bash 
EARLY_STOPPING=5

while (( $# )); do
  case "$1" in
    -p|--prefix)
        PREFIX=$2
                ;;
    -t|--templates)
        TEMPLATES=$2
                ;;
    -m|--model)
        MODEL=$2
                ;;
    -e|--early_stopping)
        EARLY_STOPPING=$2
                ;;
    -*)   echo "Bad option" >&2; exit 1 ;;
    --)   shift; args+=( "$@" ); set -- ;;
    *)    args+=( "$1" ) ;;
  esac

    shift
done

source activate my-rdkit-env

echo starting

date
echo "Command mpirun -n 16 python calculate_applicability.py --templates ${PREFIX}retro.templates.bkms.json.gz --train-smiles ${PREFIX}train.input.smiles.npy --valid-smiles ${PREFIX}valid.input.smiles.npy --test-smiles ${PREFIX}test.input.smiles.npy --output-prefix $PREFIX" 

mpirun -n 16 python calculate_applicabilty.py --templates ${PREFIX}retro.templates.${TEMPLATES}.json.gz --train-smiles ${PREFIX}train.input.smiles.npy --valid-smiles ${PREFIX}valid.input.smiles.npy --test-smiles ${PREFIX}test.input.smiles.npy --output-prefix $PREFIX
echo ending
date

for HIGHWAY in 0 1 3;
do
        for LAYERS in 1 2 3;
        do
                for HIDDEN in 1024 2048 4096;
                do
                        for REGULARIZER in "none";
                        do
        LEARNING_RATE=0.0003
        NAME="${LAYERS}x_${HIDDEN}_h${HIGHWAY}_lr${LEARNING_RATE}_r${REGULARIZER}_es${EARLY_STOPPING}"
      command="python train.py --templates ${PREFIX}retro.templates."$TEMPLATES".json.gz --train-smiles ${PREFIX}train.input.smiles.npy --valid-smiles ${PREFIX}valid.input.smiles.npy --train-labels ${PREFIX}train.labels.classes.npy --valid-labels ${PREFIX}valid.labels.classes.npy --model-name "$MODEL"-baseline"$NAME" --nproc 16 --num-hidden $LAYERS --hidden-size $HIDDEN --num-highway $HIGHWAY --learning-rate $LEARNING_RATE --regularizer $REGULARIZER --early-stopping $EARLY_STOPPING"
      echo "starting baseline model training"

      echo "Command: $command" 
      date
      $command > log.$NAME
      echo "ending"
      date

      command="python train_appl.py --train-smiles ${PREFIX}train.input.smiles.npy --valid-smiles ${PREFIX}valid.input.smiles.npy --train-labels ${PREFIX}train.appl_matrix.npz --valid-labels ${PREFIX}valid.appl_matrix.npz --model-name "$MODEL"-appl-"$NAME" --nproc 16 --num-hidden $LAYERS --hidden-size $HIDDEN --num-highway $HIGHWAY --learning-rate $LEARNING_RATE --early-stopping $EARLY_STOPPING"

      echo "Command: $command" 
      echo "starting pretraining"
      date
      $command > log_appl.$NAME

      echo "ending"
      date

      echo "starting retraining"
      date

      command="python train.py --templates ${PREFIX}retro.templates."$TEMPLATES".json.gz --train-smiles ${PREFIX}train.input.smiles.npy --valid-smiles ${PREFIX}valid.input.smiles.npy --train-labels ${PREFIX}train.labels.classes.npy --valid-labels ${PREFIX}valid.labels.classes.npy --nproc 16 --model-name "$MODEL"-pretrained-"$NAME" --pretrain-weights training/"$MODEL"-appl-"$NAME"-weights.hdf5 --num-hidden $LAYERS --hidden-size $HIDDEN --num-highway $HIGHWAY --learning-rate $LEARNING_RATE --regularizer $REGULARIZER --early-stopping $EARLY_STOPPING"


      echo "Command: $command" 

      $command > log.${NAME}_pretrain_${MODEL}

      echo ending
      date
        done
        done
    done
  done
