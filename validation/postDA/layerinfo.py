
LayerLIST = [
    [0,10],
    [10,30],
    [30,60],
    [60,100],
    [100,150],
    [150,300],
    [300,600]
]

DICTlayerQ = {}
dep5mQ = {}
DICTlayersQ = {
    'chl' : [],
    'nit' : [],
}
for ll in LayerLIST:
    layername = '%s' %(ll[0]) + '-' + '%s' %(ll[1])
    DICTlayerQ[layername] = ll
    dep5mQ[layername] = range(ll[0],ll[1]+1,5)
    DICTlayersQ['nit'].append(layername)
    if ll[1]<200:
        DICTlayersQ['chl'].append(layername)



