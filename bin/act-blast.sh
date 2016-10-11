#!/usr/bin/env bash

# Author: Thomas Hackl - thackl@lim4.de
# Date: 2016-05-12

## subs
usage(){
cat <<EOF
Usage:
  act-blast.sh A.fa B.fa "blast options .."';
EOF
exit 0;
}

check_bin(){
    hash $1 || { echo "$1 required in PATH" >&2; exit 1;}
}

## prep
[[ $# -eq 0 ]] && usage;

check_bin seq-join;
check_bin blastn;


## main
A_PRE=`basename $1 .fa`;
A_FA=/tmp/$A_PRE-single.fa
B_PRE=`basename $2 .fa`;
B_FA=/tmp/$B_PRE-single.fa

/home/thackl/projects/coding/code/seq-scripts/bin/seq-join < $1 > $A_FA;
/home/thackl/projects/coding/code/seq-scripts/bin/seq-join < $2 > $B_FA;

blastn -subject $A_FA -query $B_FA $3 \
       -outfmt '6 score pident qstart qend qseqid sstart send sseqid' |
    tr '\t' ' ' >"$A_PRE"_"$B_PRE"-act-blast.tsv
