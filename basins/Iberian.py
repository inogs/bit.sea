from basins.region import Polygon
from basins.basin import SimplePolygonalBasin, ComposedBasin, SimpleBasin
from basins.region import Rectangle

reg1 = Rectangle(4,8,41,43)
reg2 = Rectangle(-3,-0.5,35,37.5)
reg3 = Rectangle(3,6,38,40.5)
reg4 = Rectangle(0.5,2.0,36.5,38.0)

Reg1 = SimpleBasin('Reg1',reg1,'Region1')
Reg2 = SimpleBasin('Reg2',reg2,'Region2')
Reg3 = SimpleBasin('Reg3',reg3,'Region3')
Reg4 = SimpleBasin('Reg4',reg4,'Region4')


P    = ComposedBasin('P',[Reg1, Reg2, Reg3, Reg4])
