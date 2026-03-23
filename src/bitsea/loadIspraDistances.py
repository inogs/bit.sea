import numpy as np

StatNameDist = {}
Distances = {}

DIRINP = '/pico/scratch/userexternal/ateruzzi/DISTANZE_ISPRA/' 


filename = 'DistanzeARPAint.csv'
inputf = open(DIRINP + filename)
StatNameDist['int'] = []
Distances['int'] = []
for il,line in enumerate(inputf):
    if il>0:
        linelist = line.strip().split(';')
        StatNameDist['int'].append(linelist[0])
        Distances['int'].append(float(linelist[1]))

inputf.close()


filename = 'DistanzeARPAest.csv'
inputf = open(DIRINP + filename)
StatNameDist['est'] = []
Distances['est'] = []
for il,line in enumerate(inputf):
    if il>0:
        linelist = line.strip().split(';')
        StatNameDist['est'].append(linelist[0])
        Distances['est'].append(float(linelist[1]))

inputf.close()

statExcludeInt = StatNameDist['int']

# Exclude stations close to the coast
limdist = 500 #m
print('Threshold for the distance from the coast: ' + str(limdist) + 'm')
excludedist = []
for iid in range(len(Distances['est'])):
    dd = Distances['est'][iid]
    if dd<limdist: excludedist.append(StatNameDist['est'][iid])


statExcludeDist = excludedist


