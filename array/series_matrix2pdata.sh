#!/usr/bin/env bash
# series_matrix to pdata

FILE="$1"
sed -n '/!Sample_title/{:loop;p;/!Sample_data_row_count/q;N;s/.*\n//;b loop};' "$FILE" > temp
# row 2 columns
awk -F "\t" '{ for (i=1; i<=NF; i++) RtoC[i]= (RtoC[i]!=""? RtoC[i] FS $i: $i) }    END{ for (i in RtoC) print RtoC[i] }' temp | sed -E 's/\!//g' > pdata.txt
# Remove temporary file
rm temp
