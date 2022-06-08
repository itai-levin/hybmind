#!/bin/bash
while getopts d:p:o:x: option
do
case "${option}"
in
d) DIRPATH=${OPTARG};;
p) PREFIX=${OPTARG};;
o) OUTFILE=${OPTARG};;
x) OUTDIR=${OPTARG};;
esac
done

mkdir $OUTDIR/output/

cp $DIRPATH/$PREFIX*txt $OUTDIR/output/

echo "saving to $OUTDIR/$PREFIX$OUTFILE.txt"
cat $OUTDIR/output/$PREFIX*txt > $OUTDIR/$PREFIX$OUTFILE-nobackslash.txt
sed 's/?/\\/g' $OUTDIR/$PREFIX$OUTFILE-nobackslash.txt > $OUTDIR/$PREFIX$OUTFILE.txt

rm $OUTDIR/output/*
rm $OUTDIR/*nobackslash*
rmdir $OUTDIR/output
