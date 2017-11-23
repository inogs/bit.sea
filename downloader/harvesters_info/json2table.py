import sys, json
fid=open("a.json",'r') 
A=json.load(fid)  # list of dicts
fid.close()

nFloats= len(A)
FLOAT_dtype=[('id',np.int),('wmo','S20'),('id_type',np.int),('type','S20'),('nome_fs','S20'),('status','S1')]



f=open("new_wmo.txt","w")
f.write("id_float |   wmo   | id_type |     type      |  nome_fs   | status\n")
f.write("----------+---------+---------+---------------+------------+--------\n")



for iFloat in range(nFloats):
    D = A[iFloat]
    id_float=int(D['id_float'])
    wmo     =int(D['wmo'])
    id_type =int(D['id_type'])
    type=    str(D['type'])
    nomefs=  str(D['nome_fs'])
    status=  str(D['status'])
    
    
    row="%10.f\t|\t%s\t|\t%d\t|\t%s\t|\t%s\t|\t%s\n"  % (id_float ,wmo ,id_type, type, nomefs,status)
    #print row
    f.write(row)
f.close()