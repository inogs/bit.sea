from pathlib import Path

from bitsea.basins.basin import SimplePolygonalBasin, ComposedBasin
from bitsea.basins.region import Polygon

DATADIR = Path(__file__).resolve().parent / 'WKTcoords_12nm'

filename = DATADIR / 'GolfoTrieste.txt'
got = Polygon.from_WKT_file(filename)
got = SimplePolygonalBasin('got', got, 'Gulf of Trieste 12nm')

filename = DATADIR / 'LagunaGradoMarano.txt'
lgm = Polygon.from_WKT_file(filename)
lgm = SimplePolygonalBasin('lgm', lgm, 'Grado-Marano lagoon 12nm')

filename = DATADIR / 'VenetoEst.txt'
vne = Polygon.from_WKT_file(filename)
vne = SimplePolygonalBasin('vne', vne, 'Eastern Veneto 12 nm')

filename = DATADIR / 'LagunaVenezia.txt'
lve = Polygon.from_WKT_file(filename)
lve = SimplePolygonalBasin('lve', lve, 'Venice lagoon 12nm')

filename = DATADIR / 'Po.txt'
poa = Polygon.from_WKT_file(filename)
poa = SimplePolygonalBasin('poa', poa, 'Po area 12nm')

filename = DATADIR / 'EmiliaRomagna.txt'
emr = Polygon.from_WKT_file(filename)
emr = SimplePolygonalBasin('emr', emr, 'Emilia-Romagna 12nm')

filename = DATADIR / 'Marche.txt'
mar = Polygon.from_WKT_file(filename)
mar = SimplePolygonalBasin('mar', mar, 'Marche 12 nm')

filename = DATADIR / 'Abruzzo.txt'
abr = Polygon.from_WKT_file(filename)
abr = SimplePolygonalBasin('abr', abr, 'Abruzzo 12nm')

filename = DATADIR / 'Molise.txt'
mol = Polygon.from_WKT_file(filename)
mol = SimplePolygonalBasin('mol', mol, 'Molise 12nm')

# ... add other polygons


# Composed basins

fvg = ComposedBasin('fvg', [got, lgm], 'Friuli-Venezia Giulia 12 nm')
ven = ComposedBasin('ven', [vne, lve, poa], 'Veneto 12 nm')

NAd = ComposedBasin('NAd', [got, lgm, vne, lve, poa, emr, mar], 'NAdri 12 nm')

NAd_coastal_basins = (got, lgm, vne, lve, poa, emr, mar, fvg, ven, NAd)
