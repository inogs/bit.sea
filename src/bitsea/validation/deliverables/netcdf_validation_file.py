import netCDF4 as NC
from pathlib import Path

def write(outfile:Path,
          kwargs) -> Path:
    nSub,nCoast = kwargs["MODEL_MEAN"].shape
    with NC.Dataset(outfile, "w") as ncOUT:
        ncOUT.createDimension("nsub", nSub)
        ncOUT.createDimension("ncoast", nCoast)
        s = ""
        for coast in kwargs["COASTNESS_LIST"]:
            s = s + coast + ","
        setattr(ncOUT, "coastlist", s[:-1])
        s = ""
        for sub in kwargs["basinlist"]:
            s = s + sub + ","
        setattr(ncOUT,'basinlist',s[:-1])
        
        setattr(ncOUT,"modelfile",kwargs["modelfile"])
        setattr(ncOUT,"satfile",kwargs["satfile"])
        
        
        ncvar = ncOUT.createVariable("VALID_POINTS", "i", ("nsub", "ncoast"))
        ncvar[:] = kwargs["VALID_POINTS"]
        
        Metric_list=[
            "BGC_CLASS4_CHL_RMS_SURF_BASIN",
            "BGC_CLASS4_CHL_BIAS_SURF_BASIN",
            "BGC_CLASS4_CHL_CORR_SURF_BASIN",
            "MODEL_MEAN",
            "SAT___MEAN",
            "MODEL__STD",
            "SAT____STD",
            "BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG",
            "BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG",
        ]
        for metric in Metric_list:
            ncvar= ncOUT.createVariable(metric,"f", ("nsub", "ncoast"))
            ncvar[:] = kwargs[metric]
    return outfile

