from basins.region import Polygon
from basins.basin import SimplePolygonalBasin, ComposedBasin

alb = Polygon([-5.5,-1.0,-1.0,-5.5],
              [32.0,32.0,40.0,40.0])
alb = SimplePolygonalBasin('alb', alb, 'Alboran Sea')

swm = Polygon([-1.0,  9.25,   9.25, -1.0 ],
              [32.0, 32.00,  39.50, 39.5 ])
swm = SimplePolygonalBasin('swm', swm, 'South West Med')

nwm = Polygon([-1.0,9.25, 9.25,15.00,13.0, 10,-1],
              [39.5,39.5,41.25,41.25,42.5 ,46,46])
nwm = SimplePolygonalBasin('nwm', nwm, 'North West Med')

tyr = Polygon([ 9.25,15.00,16.5,16.5,15.85,15.85,15.46875,14.6,12.5,12.5, 9.25],
              [41.25,41.25,40.0,38.6,38.20,38.30,38.21875,38.0,37.8,36.8,36.80])
tyr = SimplePolygonalBasin('tyr', tyr, 'Tyrrhenian Sea')

adn = Polygon([10, 13.0, 20.0, 15],
              [46, 42.5, 42.5, 46])
adn = SimplePolygonalBasin('adn', adn)

ads = Polygon([14.0, 20.0, 20, 18.5, 18.0, 16.6, 13.0],
              [42.5, 42.5, 40, 40.0, 40.5, 41.0, 42.5])
ads = SimplePolygonalBasin('ads', ads)

aeg = Polygon([22.80, 23, 23.00, 28.00, 28,22, 22.8 ],
              [37.15, 37, 35.25, 35.25, 43,43, 37.15])
aeg = SimplePolygonalBasin('aeg', aeg, 'Aegean Sea')

ion = Polygon([16.5,16.5,15.85,15.85,15.46875,14.6,12.5,12.5,9.25,9.25,21.63,21.63,23.0,20,18.5,18.0,16.6],
              [40.0,38.6,38.20,38.30,38.21875,38.0,37.8,36.8,36.80,30.0,30.00,36.82,38.7,40,40.0,40.5,41.0])
ion = SimplePolygonalBasin('ion', ion, 'Ionian Sea')

lev = Polygon([21.63, 21.63, 21.90, 22.80, 23.00, 23.00, 28.00, 28.00, 36.30, 36.30], 
              [30.00, 36.82, 37.15, 37.15, 37.00, 35.25, 35.25, 38.00, 38.00, 30.00])
lev = SimplePolygonalBasin('lev', lev, 'Levantine Sea')

atl = Polygon([ -5.5, -5.5, -9.0, -9.0],
              [ 32.0, 40.0, 40.0, 32.0])
atl = SimplePolygonalBasin('atl', atl)


med = Polygon([ -5.5, 36.0, 36.0,  -5.5],
              [ 29.0, 29.0, 46.0,  46.0])
med = SimplePolygonalBasin('med',med,'Mediterranean Sea')
med = ComposedBasin('med',[alb ,swm,nwm,tyr,adn,ads,aeg,ion,lev],'Mediterranean Sea')

wes = ComposedBasin('wes', [alb, swm, nwm, tyr], 'West Mediterranean Sea')
eas = ComposedBasin('eas', [ion, lev], 'East Mediterranean Sea')
mnm = ComposedBasin('mnm', [alb, swm, nwm, tyr, ion, lev],"Med without marginal seas")
adr = ComposedBasin('adr',[adn,ads], 'Adriatic')

Pred = ComposedBasin('Pr',[alb ,swm,nwm,tyr,adr,aeg,ion,lev])
P    = ComposedBasin('P', [alb ,swm,nwm,tyr,adr,aeg,ion,lev,med])

