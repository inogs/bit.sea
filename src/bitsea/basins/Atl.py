from bitsea.basins.region import Polygon, Rectangle
from bitsea.basins.basin import SimplePolygonalBasin, ComposedBasin

Gib = Rectangle(-8.0,  -6.0,	35.0,	36.5)
Gib = SimplePolygonalBasin('GibBox', Gib, 'GIB 2 x 1.5 degrees')


P = ComposedBasin(
    'box',
    [Gib],
    'Gib Boxes'
)
