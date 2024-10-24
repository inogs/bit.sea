# bit.sea
![ogstm picture](https://github.com/inogs/ogstm/blob/master/DOC/PICTURES/PPN_MED_OGSTM_BFM.png)
## What is bit.sea ?
bit.sea is a python tool for scientific oceanographic applications,
to manage the [preprocessing](https://github.com/inogs/ogstm_preproc) and [postprocessing](https://github.com/inogs/ogstm/ogstm_postproc) of [ogstm](https://github.com/inogs/ogstm).




The tool is employed in the context of operational oceanography within [Marine Copernicus](https://data.marine.copernicus.eu/product/MEDSEA_ANALYSISFORECAST_BGC_006_014/description)
to generate quality checked observational datasets used in production and to manage the [validation framework](https://catalogue.marine.copernicus.eu/documents/QUID/CMEMS-MED-QUID-006-014.pdf).

## Installation:

```bash
git clone git@github.com:inogs/bit.sea.git
cd bit.sea
pip install bit.sea .
```

## Usage

```bash
import bit.sea
from bitsea.commons import genUserDateList as DL
from bitsea.commons import timerequestors, TimeList
Min15 = DL.getTimeList("20180301-00:00:00","20200310-00:00:00", minutes=15)
TL     = TimeList(Min15)
Hourly_req=requestors.Hourly_req(2018,3,5,12,delta_hours=2)
indexes,weights = TL.select(Hourly_req)
```


The software responsible is [@gbolzon](https://www.github.com/gbolzon).