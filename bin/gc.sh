#!/bin/bash
echo -e "id\tlen\tbp\tgc\tfgc"
for LIB in "$@"; do 
    ID=`echo $LIB | cut -d_ -f3`; 
    LEN=`echo $LIB | cut -d. -f3`; 
    GC=`sed -n '2~4p' "$LIB" | perl -ne 'chomp; $t+=length($_); $gc+= $_=~tr/gcGC//; END{printf("%s\t%s\t%.4f\n", $t, $gc, $gc/$t)}'`; 
    printf "%s\t%s\t%s\t%s\t%s\n" $ID $LEN $GC; 
done;
