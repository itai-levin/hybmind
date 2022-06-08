#!/bin/bash

while getopts i:t:o:c: option
do
case "${option}"
in
i) SMILESFILE=${OPTARG};; #name of input file containing Reaction IDs in column 1 and SMILES in column 2
t) TASKNUM=${OPTARG};; #Number of nodes to split the mappings across
o) OUTFILE=${OPTARG};; #Prefix for the output files
c) COLUMN=${OPTARG};; #Column number for unmapped smiles
esac
done

LINENUM=$(wc -l $SMILESFILE| cut -f4 -d' ')
BINSIZE=$(expr $LINENUM / $TASKNUM + 1)

echo $LINENUM
echo $BINSIZE
echo $TASKNUM
echo $OUTFILE
echo $COLUMN
echo $SMILESFILE

for (( LOWER=1; LOWER<$LINENUM; LOWER+=$BINSIZE )); do
        UPPER=$(expr $LOWER + $BINSIZE - 1)
        bash map_smiles.sh -i $SMILESFILE -l $LOWER -u $UPPER -o "$OUTFILE-$LOWER-$UPPER.txt" -c $COLUMN #&
done
