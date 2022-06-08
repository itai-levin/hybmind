#!/bin/bash
  
#IDEALLY LINENUM WOULD BE SET AUTOMATICALLY BUT NOT SURE HOW TO GET THE NUMBER OF LINES FROM THE JSON
LINENUM=188

TASKNUM=10

BINSIZE=$(expr $LINENUM / $TASKNUM + 1)
TASKID=0

for (( LOWER=0; LOWER<$LINENUM ; LOWER+=$BINSIZE )); do
        UPPER=$(expr $LOWER + $BINSIZE)
        echo TASKID $TASKID, LOWER $LOWER, UPPER $UPPER
        python apply_all_reaxys_to_bkms_subset_rdchiral.py --start $LOWER --end $UPPER --task-id $TASKID --reactions ../data/reactions_with_over_conjugated_phosphate_fixed.json.gz --templates ../data/reaxys_templates_wcharge-wostereo_wocharge-wstereo_wocharge-wostereo.json.gz --prefix ../data/nochargeorstereo &
        TASKID=$(expr $TASKID + 1)
done