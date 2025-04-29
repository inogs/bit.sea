from bitsea.basins.region import Polygon, Rectangle
from bitsea.basins.basin import SimplePolygonalBasin, ComposedBasin

# Boxes of 0.5 degrees or composed boxes or polygons for rivers with discharges higher than 100 m3s-1 in 2019

Ebro05 = Rectangle(0.6666667,  1.1666667,	40.4791679,	40.9791679)
Ebro05 = SimplePolygonalBasin('Ebro', Ebro05, 'Ebro box 0.5 degrees')

Po05 = Rectangle(12.2500000,	12.7500000,	44.6041679,	45.1041679)
Po05 = SimplePolygonalBasin('Po', Po05, 'Po box 0.5 degrees')

PoAdigeBrenta05 = Rectangle(12.2500000,	12.7500000,	44.6041679,	45.2708321)
PoAdigeBrenta05 = SimplePolygonalBasin('Po+Adige+Brenta', PoAdigeBrenta05, 'Po-Adige-Brenta composed box')

GrandRhone05 = Rectangle(4.5833335,	5.0833335,	43.0625000,	43.5625000)
GrandRhone05 = SimplePolygonalBasin('GrandRhone', GrandRhone05, 'Grand Rhone box 0.5 degrees')

PetitRhone05 = Rectangle(4.1250000,	4.6250000,	43.1875000,	43.6875000)
PetitRhone05 = SimplePolygonalBasin('PetitRhone', PetitRhone05, 'Petit Rhone box 0.5 degrees')

SemanVjoseShkumbini05 = Rectangle(19.1250000,	19.6250000,	40.5625000,	41.0625000)
SemanVjoseShkumbini05 = SimplePolygonalBasin('Seman+Vjose+Shkumbini', SemanVjoseShkumbini05, 'Seman-Vjose-Shkumbini box 0.5 degrees')

BunaBojanaMati05  = Rectangle(19.1250000,	19.6250000,	41.5625000,	42.0625000)
BunaBojanaMati05  = SimplePolygonalBasin('Buna_Bojana+Mati', BunaBojanaMati05 , 'Buna/Bojana and Mati box 0.5 degrees')

Piave05  = Rectangle(12.4583330,	12.9583330,	45.2708321,	45.7708321)
Piave05  = SimplePolygonalBasin('Piave', Piave05 , 'Piave box 0.5 degrees')

SocaIsonzo05  = Rectangle(13.2916670, 13.7916670,	45.4791679,	45.9791679)
SocaIsonzo05  = SimplePolygonalBasin('Soca_Isonzo', SocaIsonzo05 , 'Soca/Isonzo box 0.5 degrees')

Arno05  = Rectangle(10.0000000,	10.5000000,	43.4375000,	43.9375000)
Arno05  = SimplePolygonalBasin('Arno', Arno05 , 'Piave box 0.5 degrees')

#Neretva05  = Rectangle(17.1250000,	17.6250000,	42.7708321,	43.2708321)
#Neretva05  = SimplePolygonalBasin('Neretva05', Neretva05 , 'Neretva box 0.5 degrees')

Neretva05red  = Rectangle(17.1250000,	17.6250000,	42.983,	43.2708321)
Neretva05red  = SimplePolygonalBasin('Neretva', Neretva05red, 'Neretva box 0.5 degrees reduced')

Tevere05  = Rectangle(11.9583330,	12.4583330,	41.4791679,	41.9791679)
Tevere05  = SimplePolygonalBasin('Tevere', Tevere05 , 'Tevere box 0.5 degrees')

Volturno05  = Rectangle(13.6250000,	14.1250000,	40.7708321,	41.2708321)
Volturno05  = SimplePolygonalBasin('Volturno', Volturno05, 'Volturno box 0.5 degrees')

Meric05 = Rectangle(25.7500000,	26.2500000,	40.4791679,	40.9791679)
Meric05 = SimplePolygonalBasin('Meric_Evros_Maritsa', Meric05, 'Meric/Evros/Maritsa box 0.5 degrees')

Axios05 = Rectangle(22.4583340,	22.9583340,	40.2291679,	40.7291679)
Axios05 = SimplePolygonalBasin('Axios_Vadar', Axios05, 'Axios/Vadar box 0.5 degrees')

Acheloos05 = Rectangle(20.8333340,	21.3333340,	38.1041679,	38.6041679)
Acheloos05 = SimplePolygonalBasin('Acheloos', Acheloos05, 'Acheloos box 0.5 degrees')

Gediz05 = Rectangle(26.5000000,	27.0000000,	38.3541679,	38.8541679)
Gediz05 = SimplePolygonalBasin('Gediz', Gediz05, 'Gediz box 0.5 degrees')

#Buyuk05 = Rectangle(26.9166660,	27.4166660,	37.2708321,	37.7708321)
#Buyuk05 = SimplePolygonalBasin('Buyuk-Menderes05', Buyuk05, 'Buyuk box 0.5 degrees')

Buyuk05red = Rectangle(26.916666,	27.416666,	37.2708321,	37.683)
Buyuk05red = SimplePolygonalBasin('BuyukMenderes', Buyuk05red, 'Buyuk box 0.5 degrees reduced')

Ceyhan05 = Rectangle(35.2916679,	35.7916679,	36.3125000,	36.8125000)
Ceyhan05 = SimplePolygonalBasin('Ceyhan', Ceyhan05, 'Ceyhan box 0.5 degrees')

Goksu05 = Rectangle(33.7916679,	34.2916679,	36.0208321,	36.5208321)
Goksu05 = SimplePolygonalBasin('Goksu', Goksu05, 'Goksu box 0.5 degrees')

Medjerda05 = Rectangle(9.9583330,	10.4583330,	36.7708321,	37.2708321)
Medjerda05 = SimplePolygonalBasin('Medjerda', Medjerda05, 'Goksu box 0.5 degrees')

#PoAdigeBrenta05large = Polygon(
#    [12.25,12.75,12.75,12.833333,12.833333,12.75,12.75,12.25],
#    [44.6041679,44.6041679 ,44.6875 ,44.6875 ,45.1875 ,45.1875 ,45.2708321 ,45.2708321]
#)
#PoAdigeBrenta05large = SimplePolygonalBasin(
#    'PoAdigeBrenta05large',
#    PoAdigeBrenta05large,
#    'Po-Adige-Brenta composed large box'
#)

P = ComposedBasin(
    'rivers',
    [Ebro05, Po05, PoAdigeBrenta05, GrandRhone05, PetitRhone05, SemanVjoseShkumbini05, BunaBojanaMati05, Piave05, SocaIsonzo05,Arno05, Neretva05red,Tevere05,Volturno05, Meric05, Axios05, Acheloos05 , Gediz05,  Buyuk05red, Ceyhan05, Goksu05,Medjerda05],
    'River boxes'
)
