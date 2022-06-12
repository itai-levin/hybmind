#!/bin/bash 
#SBATCH -N 1 
#SBARCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --time=0-01:00:00
#SBATCH --partition=sched_mit_ccoley
#SBATCH --nodelist=node1238 
#SBATCH --mem=50000

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

for LAYERS in 1;
do
  for HIDDEN in 4096;
  do
      for REGULARIZER in "none"
      do
        LEARNING_RATE=0.0003
        NAME="${LAYERS}x_${HIDDEN}_h${HIGHWAY}_lr${LEARNING_RATE}_r${REGULARIZER}_es${EARLY_STOPPING}_noval"
      	HIGHWAY=0
      command="python train.py --templates ${PREFIX}retro.templates."$TEMPLATES".json.gz --train-smiles ${PREFIX}train,valid,test.input.smiles.npy --no-validation --train-labels ${PREFIX}train,valid,test.labels.classes.npy --model-name "$MODEL"-baseline-"$NAME" --nproc 16 --num-hidden $LAYERS --hidden-size $HIDDEN --num-highway $HIGHWAY --learning-rate $LEARNING_RATE --regularizer $REGULARIZER --epochs 9"
      echo "starting baseline model training"

      echo "Command: $command" 
      date
      $command > log.$NAME
      echo "ending"
      date

     command="python train_appl.py --train-smiles ${PREFIX}train,valid,test.input.smiles.npy --train-labels ${PREFIX}train,valid,test.appl_matrix.npz --model-name "$MODEL"-appl-"$NAME" --nproc 16 --num-hidden $LAYERS --hidden-size $HIDDEN --num-highway $HIGHWAY --learning-rate $LEARNING_RATE  --no-validation --epochs 25"


      echo "Command: $command" 
      echo "starting pretraining"
      date
      $command > log_appl.$NAME

      echo "ending"
      date

      echo "starting retraining"
      date

      command="python train.py --templates ${PREFIX}retro.templates."$TEMPLATES".json.gz --train-smiles ${PREFIX}train,valid,test.input.smiles.npy --train-labels ${PREFIX}train,valid,test.labels.classes.npy --nproc 16 --model-name "$MODEL"-pretrained-"$NAME" --pretrain-weights training/"$MODEL"-appl-"$NAME"-weights.hdf5 --num-hidden $LAYERS --hidden-size $HIDDEN --num-highway $HIGHWAY --learning-rate $LEARNING_RATE --regularizer $REGULARIZER --no-validation --epochs 10"


      echo "Command: $command" 

      command > log.${NAME}_pretrain_${MODEL}

      echo ending
      date 
      done
    done
  done

