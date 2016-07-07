from basins.region import Polygon
from basins.basin import SimplePolygonalBasin, ComposedBasin

alb = Polygon([-5.5,-1.0,-1.0,-5.5],
              [32.0,32.0,40.0,40.0])
alb = SimplePolygonalBasin('alb', alb, 'Alboran Sea')

swm1 = Polygon([-1.0, 5.00,  5.00,-1.00 ],
               [32.0,32.00, 39.50,39.50 ])
swm1 = SimplePolygonalBasin('swm1', swm1,'South Western Mediterranean west')

swm2 = Polygon([5, 9.25,  9.25, 5 ],
               [32.0,32.00, 39.50,39.5 ])
swm2 = SimplePolygonalBasin('swm2', swm2,'South Western Mediterranean east')

nwm = Polygon([-1.0,9.25, 9.25,-1],
        [39.5,39.5, 46,46])
nwm = SimplePolygonalBasin('nwm', nwm,'North Western Mediterranean ')

tyr1 = Polygon([09.25,15.00,10,09.25],
               [41.25,41.25,46,46.00])
tyr1 = SimplePolygonalBasin('tyr1', tyr1,'Northern Tyrrhenian')


tyr2 = Polygon([ 9.25, 15.00, 15.00, 15.60,  16.10, 16.50, 15.00, 9.25],
               [36.75, 36.75, 38.00, 38.20,  38.20, 39.50, 41.25, 41.25] )
tyr2 = SimplePolygonalBasin('tyr2', tyr2,'Southern Tyrrhenian')

ion1 = Polygon([ 9.25, 15.00, 15.00, 10.70,  9.25],
               [32.00, 32.00, 36.75, 36.75, 35.00])
ion1 = SimplePolygonalBasin('ion1', ion1,'Western Ionian')

ion2 = Polygon([15.00, 21.85, 21.85, 15.00],
               [30.00, 30.00, 36.75, 36.75])
ion2 = SimplePolygonalBasin('ion2', ion2,'Eastern Ionian')

ion3 = Polygon([15.00, 21.85, 21.85, 18.50, 17.00, 16.10, 16.50, 16.10, 15.60, 15.00],
               [36.75, 36.75, 40.00, 40.00, 41.00, 40.00, 39.50, 38.20, 38.20, 38.00])
ion3 = SimplePolygonalBasin('ion3', ion3,'Northern Ionian')

adr1 = Polygon([10,13.0,20.0, 15],
               [46,42.5,42.5, 46])
adr1 = SimplePolygonalBasin('adr1', adr1,'Northern Adriatic')

adr2 = Polygon([ 14.0,20.0,20,18.5,18.0,16.6,13.0],
               [ 42.5,42.5,40,40.0,40.5,41.0,42.5])
adr2 = SimplePolygonalBasin('adr2', adr2,'Southern Adriatic')

lev1 = Polygon([21.85, 26.25, 26.25, 24.90, 24.00, 21.85 ],
               [30.00, 30.00, 35.10, 35.10, 35.30, 35.30 ]) 
lev1 = SimplePolygonalBasin('lev1', lev1,'Western Levantine')

lev2 = Polygon([26.25, 33.00, 33.00, 28.00, 28.00, 26.30, 26.25 ],
               [33.60, 33.60, 38.00, 38.00, 35.30, 35.30, 35.28 ]) 
lev2 = SimplePolygonalBasin('lev2', lev2,'Northern Levantine')

lev3 = Polygon([26.25, 26.25, 33.00, 33.00],
               [30.00, 33.60, 33.60, 30.00]) 
lev3 = SimplePolygonalBasin('lev3', lev3,'Southern Levantine')

lev4 = Polygon([33.00, 37.00, 37.00, 33.00],
               [30.00, 30.00, 38.00, 38.00]) 
lev4 = SimplePolygonalBasin('lev4', lev4,'Eastern Levantine')

aeg = Polygon([21.85, 24.00, 24.90, 26.25, 26.25, 26.30, 28.00, 28.00, 21.85],
        [35.30, 35.30, 35.10, 35.10, 35.28, 35.30, 35.30, 42.00, 42.00])
aeg = SimplePolygonalBasin('aeg', aeg, 'Aegean Sea')


atl  = Polygon([ -5.5,-5.5,-9.0,-9.0], 
               [ 32.0,40.0,40.0,32.0]) 
atl  = SimplePolygonalBasin('atl', atl,'Atlantic buffer')

med = ComposedBasin('med',[alb ,swm1,swm2,nwm,tyr1,tyr2,adr1,adr2,aeg,ion1,ion2,ion3,lev1,lev2,lev3,lev4],'Mediterranean Sea')

wes = ComposedBasin('wes', [alb ,swm1,swm2,nwm,tyr1,tyr2], 'West Mediterranean Sea')
eas = ComposedBasin('eas', [ion1,ion2,ion3,lev1,lev2,lev3,lev4], 'East Mediterranean Sea')
mnm = ComposedBasin('mnm', [alb ,swm1,swm2,nwm,tyr1,tyr2,ion1,ion2,ion3,lev1,lev2,lev3,lev4],"Med without marginal seas")

Pred = ComposedBasin('Pr',[alb ,swm1,swm2,nwm,tyr1,tyr2,adr1,adr2,aeg,ion1,ion2,ion3,lev1,lev2,lev3,lev4])
P    = ComposedBasin('P',[alb ,swm1,swm2,nwm,tyr1,tyr2,adr1,adr2,aeg,ion1,ion2,ion3,lev1,lev2,lev3,lev4,med])

#for NRT3

lev = ComposedBasin('lev',[lev1, lev2, lev3, lev4],'Levantine Sea')
ion = ComposedBasin('ion',[ion1, ion2, ion3],'Ionian Sea')
adr = ComposedBasin('adr',[adr1, adr2],'Adriatic Sea')
tyr = ComposedBasin('tyr',[tyr1, tyr2],'Tyrrenian Sea')
swm = ComposedBasin('swm',[swm1, swm2], 'South Western Mediterraneaan Sea')

NRT3 = ComposedBasin('NRT3',[alb, swm, nwm, tyr, adr, ion, lev],'Gruped Subbasin for Near Real Time medeaf page')


