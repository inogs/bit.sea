#! /bin/bash

QUERY="SELECT DISTINCT tbl_float.id_float, tbl_float.wmo, tbl_float.id_type, tbl_type.type, tbl_float.nome_lov, tbl_float.status FROM ((tbl_float LEFT JOIN tbl_float_sensori_web ON tbl_float.id_float = tbl_float_sensori_web.id_float) LEFT JOIN tbl_type ON tbl_float.id_type = tbl_type.id_type) LEFT JOIN tbl_sensori_web ON tbl_float_sensori_web.id_sensori_web = tbl_sensori_web.id_sensori_web WHERE (((tbl_sensori_web.nome_sensore)='DO' Or (tbl_sensori_web.nome_sensore)='BIO' Or (tbl_sensori_web.nome_sensore)='NUT')) ORDER BY tbl_float.id_type,tbl_float.id_float"

#echo ${QUERY// /%20}
curl "http://maos.inogs.it/api/api.php?key=ogs112211&sql=${QUERY// /%20}"
