from basins.region import Polygon, Rectangle
from basins.basin import SimplePolygonalBasin, ComposedBasin

Ebro = Rectangle(0.41666669, 1.41666669, 40.22916794, 41.22916794)
Ebro = SimplePolygonalBasin('Ebro', Ebro, 'Ebro box')

Po = Rectangle(11.91666698, 12.91666698, 44.27083206, 45.27083206)
Po = SimplePolygonalBasin('Po', Po, 'Po box')

Rhone = Rectangle(4.33333349, 5.33333349, 42.81250000, 43.81250000)
Rhone = SimplePolygonalBasin('Rhone', Rhone, 'Rhone box')

Buna = Rectangle(18.87500000, 19.87500000, 41.31250000, 42.31250000)
Buna = SimplePolygonalBasin('Buna', Buna, 'Buna box')

Neretva = Rectangle(16.87500000, 17.87500000, 42.52083206, 43.52083206)
Neretva = SimplePolygonalBasin('Neretva', Neretva, 'Neretva box')

Meric = Rectangle(25.50000000, 26.50000000, 40.22916794, 41.22916794)
Meric = SimplePolygonalBasin('Meric', Meric, 'Meric/Evros/Maritsa box')

Manavgat = Rectangle(30.95833397, 31.95833397, 36.22916794, 37.22916794)
Manavgat = SimplePolygonalBasin('Manavgat', Manavgat, 'Manavgat box')

GoskuSeyhanCeyhan = Polygon(
    [33.54166794, 34.54166794, 34.54166794, 35.04166794, 35.04166794, 36.04166794, 36.04166794, 33.54166794],
    [35.77083206, 35.77083206, 36.22916794, 36.22916794, 36.06250000, 36.06250000, 37.22916794, 37.22916794]
)
GoskuSeyhanCeyhan = SimplePolygonalBasin(
    'GoskuSeyhanCeyhan',
    GoskuSeyhanCeyhan,
    'Gosku-Seyhan-Ceyhan box'
)

P = ComposedBasin(
    'rivers',
    [Ebro, Po, Rhone, Buna, Neretva, Meric, Manavgat, GoskuSeyhanCeyhan],
    'River boxes'
)
