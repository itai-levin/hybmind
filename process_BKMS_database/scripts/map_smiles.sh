#!/bin/bash 

COLUMN=2

while getopts i:l:u:o:c: option
do
case "${option}"
in
i) SMILESFILE=${OPTARG};;
l) LOWER=${OPTARG};;
u) UPPER=${OPTARG};;
o) OUTFILE=${OPTARG};;
c) COLUMN=${OPTARG};;
esac
done

CURRDIR=$(pwd)
[ ! -d mapped ] && mkdir mapped

for i in $(awk -v l=$LOWER -v u=$UPPER -F"\t" '{if (NR >= l && NR <= u) print NR}' $SMILESFILE); do
    smiles=\"$(awk -F"\t" -v var1=$i -v var2=$COLUMN 'NR==var1 {print $var2}' $SMILESFILE)\";
    echo "SMILES: $smiles"
    ids=$(awk -F"\t" -v var=$i 'NR==var {print $1}' $SMILESFILE)
    echo "ID: $id"
    java -jar ReactionDecoder.jar -Q SMI -q $smiles -j AAM -f TEXT -p mapped/$ids 

    rm mapped/${ids}_ECBLAST_smiles_AAM.rxn
    mapsmiles=$(cat mapped/${ids}_ECBLAST_smiles_AAM.txt |head -n 4 | tail -n 1 | sed 's/\\/?/g')
    awk -F"\t" -v var1=$i -v var2=$mapsmiles 'NR==var1 {print $0,"\t", var2}' $SMILESFILE >> $OUTFILE
    rm mapped/${ids}_ECBLAST_smiles_AAM.txt
    done
