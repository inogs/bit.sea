from bitsea.basins.basin import ComposedBasin
from bitsea.basins.basin import SimplePolygonalBasin
from bitsea.basins.region import Polygon
from bitsea.basins.region import Rectangle

alb = Rectangle(lonmin=-5.5, lonmax=-1.0, latmin=32.0, latmax=40.0)
alb = SimplePolygonalBasin("alb", alb, "Alboran Sea")

swm1 = Rectangle(lonmin=-1.0, lonmax=5.0, latmin=32.0, latmax=39.5)
swm1 = SimplePolygonalBasin("swm1", swm1, "South Western Mediterranean west")

swm2 = Rectangle(lonmin=5.0, lonmax=9.25, latmin=32.0, latmax=39.5)
swm2 = SimplePolygonalBasin("swm2", swm2, "South Western Mediterranean east")

nwm = Rectangle(lonmin=-0.7, lonmax=9.25, latmin=39.5, latmax=44.7)
nwm = SimplePolygonalBasin("nwm", nwm, "North Western Mediterranean ")

tyr1 = Rectangle(lonmin=9.25, lonmax=15.0, latmin=41.25, latmax=46.0)
tyr1 = SimplePolygonalBasin("tyr1", tyr1, "Northern Tyrrhenian")

tyr2 = Polygon(
    [9.25, 15.00, 15.00, 15.60, 16.10, 16.50, 15.00, 9.25],
    [36.75, 36.75, 38.00, 38.20, 38.20, 39.50, 41.25, 41.25],
)
tyr2 = SimplePolygonalBasin("tyr2", tyr2, "Southern Tyrrhenian")

ion1 = Polygon(
    [9.25, 15.00, 15.00, 10.70, 9.25], [32.00, 32.00, 36.75, 36.75, 35.00]
)
ion1 = SimplePolygonalBasin("ion1", ion1, "Western Ionian")

ion2 = Rectangle(lonmin=15.0, lonmax=21.85, latmin=30.0, latmax=36.75)
ion2 = SimplePolygonalBasin("ion2", ion2, "Eastern Ionian")

ion3 = Polygon(
    [15.00, 21.85, 21.85, 18.50, 17.00, 16.10, 16.50, 16.10, 15.60, 15.00],
    [36.75, 36.75, 40.00, 40.00, 41.00, 40.00, 39.50, 38.20, 38.20, 38.00],
)
ion3 = SimplePolygonalBasin("ion3", ion3, "Northern Ionian")

adr1 = Polygon([10, 13.0, 20.0, 15], [46, 42.5, 42.5, 46])
adr1 = SimplePolygonalBasin("adr1", adr1, "Northern Adriatic")

adr2 = Polygon(
    [14.0, 20.0, 20, 18.5, 18.0, 16.6, 13.0],
    [42.5, 42.5, 40, 40.0, 40.5, 41.0, 42.5],
)
adr2 = SimplePolygonalBasin("adr2", adr2, "Southern Adriatic")

lev1 = Polygon(
    [21.85, 26.25, 26.25, 24.90, 24.00, 21.85],
    [30.00, 30.00, 35.10, 35.10, 35.30, 35.30],
)
lev1 = SimplePolygonalBasin("lev1", lev1, "Western Levantine")

lev2 = Polygon(
    [26.25, 33.00, 33.00, 28.00, 28.00, 26.30, 26.25],
    [33.60, 33.60, 38.00, 38.00, 35.30, 35.30, 35.28],
)
lev2 = SimplePolygonalBasin("lev2", lev2, "Northern Levantine")

lev3 = Rectangle(lonmin=26.25, lonmax=33.00, latmin=30.0, latmax=33.6)
lev3 = SimplePolygonalBasin("lev3", lev3, "Southern Levantine")

lev4 = Rectangle(lonmin=33.00, lonmax=37.00, latmin=30.0, latmax=38.0)
lev4 = SimplePolygonalBasin("lev4", lev4, "Eastern Levantine")

aeg = Polygon(
    [21.85, 24.00, 24.90, 26.25, 26.25, 26.30, 28.00, 28.00, 21.85],
    [35.30, 35.30, 35.10, 35.10, 35.28, 35.30, 35.30, 42.00, 42.00],
)
aeg = SimplePolygonalBasin("aeg", aeg, "Aegean Sea")


atl = Rectangle(lonmin=-9.0, lonmax=-6.0, latmin=30.0, latmax=43.0)
atl = SimplePolygonalBasin("atl", atl, "Atlantic buffer")

med = ComposedBasin(
    "med",
    [
        alb,
        swm1,
        swm2,
        nwm,
        tyr1,
        tyr2,
        adr1,
        adr2,
        aeg,
        ion1,
        ion2,
        ion3,
        lev1,
        lev2,
        lev3,
        lev4,
    ],
    "Mediterranean Sea",
)

wes = ComposedBasin(
    "wes", [alb, swm1, swm2, nwm, tyr1, tyr2], "West Mediterranean Sea"
)
eas = ComposedBasin(
    "eas", [ion1, ion2, ion3, lev1, lev2, lev3, lev4], "East Mediterranean Sea"
)
eas2 = ComposedBasin(
    "eas2",
    [ion1, ion2, ion3, lev1, lev2, lev3, lev4, adr1, adr2, aeg],
    "East Mediterranean Sea with marginal seas",
)
mnm = ComposedBasin(
    "mnm",
    [
        alb,
        swm1,
        swm2,
        nwm,
        tyr1,
        tyr2,
        ion1,
        ion2,
        ion3,
        lev1,
        lev2,
        lev3,
        lev4,
    ],
    "Med without marginal seas",
)

# Mediterranean Sea divided in 3 areas: West, Central, East:

eas3 = ComposedBasin(
    "eas3", [lev1, lev2, lev3, lev4, aeg], "Eastern Mediterranean Sea"
)
wes3 = ComposedBasin(
    "wes3", [alb, nwm, tyr1, tyr2, swm1, swm2], "Western Mediterranean Sea"
)
mid3 = ComposedBasin(
    "mid3", [adr1, adr2, ion1, ion2, ion3], "Central Mediterranean Sea"
)
med3 = ComposedBasin("MED3", [eas3, wes3, mid3, med])
# med3 = ComposedBasin('MED3', [eas3,wes3,mid3])

# Mediterranean Sea divided in 4 subregions: --> for density-nitrate rel. analysis
eas4 = ComposedBasin(
    "East", [lev1, lev2, lev3, lev4, aeg], "Eastern Mediterranean Sea"
)
wes4 = ComposedBasin(
    "West", [alb, nwm, swm1, swm2], "Western Mediterranean Sea"
)
midw4 = ComposedBasin("midW", [tyr1, tyr2, ion1], "Mid West Mediterranean Sea")
mide4 = ComposedBasin(
    "midE", [adr1, adr2, ion2, ion3], "Mid East Mediterranean Sea"
)
med4 = ComposedBasin("MED4", [eas4, wes4, midw4, mide4])

Pred = ComposedBasin(
    "Pr",
    [
        alb,
        swm1,
        swm2,
        nwm,
        tyr1,
        tyr2,
        adr1,
        adr2,
        aeg,
        ion1,
        ion2,
        ion3,
        lev1,
        lev2,
        lev3,
        lev4,
        atl,
    ],
)
P = ComposedBasin(
    "P",
    [
        alb,
        swm1,
        swm2,
        nwm,
        tyr1,
        tyr2,
        adr1,
        adr2,
        aeg,
        ion1,
        ion2,
        ion3,
        lev1,
        lev2,
        lev3,
        lev4,
        med,
        atl,
    ],
)

# for NRT3

lev = ComposedBasin("lev", [lev1, lev2, lev3, lev4], "Levantine Sea")
ion = ComposedBasin("ion", [ion1, ion2, ion3], "Ionian Sea")
adr = ComposedBasin("adr", [adr1, adr2], "Adriatic Sea")
tyr = ComposedBasin("tyr", [tyr1, tyr2], "Tyrrenian Sea")
swm = ComposedBasin("swm", [swm1, swm2], "South Western Mediterraneaan Sea")

NRT3 = ComposedBasin(
    "NRT3",
    [alb, swm, nwm, tyr, adr, ion, lev],
    "Gruped Subbasin for Near Real Time medeaf page",
)
MVR = ComposedBasin(
    "MVR ",
    [alb, swm, nwm, tyr, adr, ion, lev, med],
    "Grouped Subbasin for Monthly Validation Report",
)
