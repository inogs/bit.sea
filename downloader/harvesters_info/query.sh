#! /bin/bash

FLAGS="-P pager=off " #-t 

QUERY="SELECT tbl_float.id_float, tbl_float.wmo, tbl_float.id_type, tbl_type.type, tbl_float.nome_fs, tbl_float.status FROM tbl_float INNER JOIN tbl_type ON tbl_float.id_type = tbl_type.id_type WHERE (((tbl_float.id_type)=64 Or (tbl_float.id_type)=62 Or (tbl_float.id_type)=59 Or (tbl_float.id_type)=58 Or (tbl_float.id_type)=57 Or (tbl_float.id_type)=52 Or (tbl_float.id_type)=49 Or (tbl_float.id_type)=44))  ORDER BY tbl_float.wmo;"

psql -h oceano.inogs.it -p 5432 -U ogspub $FLAGS -d float -c "$QUERY" > wmo.txt
## just cutting last two lines ###
NL=`awk 'END{print NR}' wmo.txt`
head -$(( NL -2 )) wmo.txt > junk.txt
mv junk.txt wmo.txt
#################################
