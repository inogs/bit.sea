#! /bin/bash


# id_type    type
# 1    Apex
# 2    Provor CTS2
# 3    Provor CTS3
# 4    Nemo
# 7    Arvor I
# 8    Apex APF8
# 23    Apex SBE
# 24    Apex APF9A_2207a
# 25    Arvor-a3
# 26    Arvor I - 2
# 28    Arvor C
# 29    Arvor-L
# 31    Arvor
# 41    Provor-BIO
# 43    Apex - rudics
# 44    Provor-BIO
# 45    Apex SBE
# 49    Provor-DO
# 51    NOVA
# 52    Provor BIO
# 53    Provor
# 57    Provor-DO-I
# 58    Provor NUT
# 59    Provor-DO
# 61    NOVA
# 62    DOVA
# 63    APEX I
# 64    Provor NUT
# 66    Provor CTS2
# 67    Arvor C
# 68    NEMO
# 69    APEX
# 70    Apex APF9A_2207B&C
# 71    Apex APF9A_2361.2
# 74    Arvor deep
# 75    Arvor
# 76    Arvor I - v2015
# 77    Arvor_N
# 78    Provor_III
# 79    Arvor C sbd
# 80    DOVA
# 81    Arvor I - 2600
# 82    Arvor I - ICE
# 88    Arvor-DO
# 126    Provor CTS3
# 127    Provor III NC
# 129    Provor BIO CTS4
# 131    Arvor deep
# 132    Arvor I DO
# 133    Provor NUT CTS4
# 135    Provor CTS5
# 136    Provor NUT CTS4
# 137    Arvor-L

for TYPE in 44 49 52 57 58 59 62 64 80 88 129 132 133 135 136 137; do

curl "http://maos.inogs.it/api/api.php?key=ogs112211&sql=SELECT%20tbl_float.id_float,%20tbl_float.wmo,%20tbl_float.id_type,%20tbl_type.type,%20tbl_float.nome_fs,%20tbl_float.status%20FROM%20tbl_float%20INNER%20JOIN%20tbl_type%20ON%20tbl_float.id_type%20=%20tbl_type.id_type%20WHERE%20((tbl_float.id_type)=${TYPE})%20ORDER%20BY%20tbl_float.wmo;"

done | sed -e "s/\]\[/,/g" 

# multiple calls of curl generate a stdout with some "][", which needed to be replaced by a ","
