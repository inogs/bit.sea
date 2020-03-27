
LayerLIST = [
    [0,10],
    [10,30],
    [30,60],
    [60,100],
    [100,150],
    [150,300],
    [300,600]
]

DICTlayer = {}
dep5m = {}
DICTlayers = {
    'chl' : [],
    'nit' : [],
}
for ll in LayerLIST:
    layername = '%s' %(ll[0]) + '-' + '%s' %(ll[1])
    DICTlayer[layername] = ll
    dep5m[layername] = range(ll[0],ll[1]+1,5)
    DICTlayers['nit'].append(layername)
    if ll[1]<200:
        DICTlayers['chl'].append(layername)



