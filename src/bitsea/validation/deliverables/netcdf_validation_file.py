import netCDF4 as NC
from pathlib import Path
import xarray as xr
from bitsea.commons.Timelist import TimeList, TimeInterval
import numpy as np

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

def read(inputfile:Path)-> xr.Dataset:
    return xr.open_dataset(inputfile)

class dir_reader():
    def __init__(self, TI: TimeInterval, dirname:Path, var:str, coastness:str):
        print(f"reading {var}")
        TL = TimeList.fromfilenames(TI, dirname, "valid.*", filtervar=var, prefix="valid.", dateformat="%Y%m%d")
        TIMES=TL.Timelist
        nFrames = TL.nTimes
        First = read(TL.filelist[0])
        nSub, nCoast=First.VALID_POINTS.shape
        iCoast=First.coastlist.rsplit(",").index(coastness)
        MODEL_MEAN = np.zeros((nFrames,nSub), np.float32)
        SAT___MEAN = np.zeros((nFrames,nSub), np.float32)
        MODEL__STD = np.zeros((nFrames,nSub), np.float32)
        SAT____STD = np.zeros((nFrames,nSub), np.float32)

        for iFrame, filename in enumerate(TL.filelist):
            H=read(filename)
            MODEL_MEAN[iFrame,:] = H.MODEL_MEAN[:,iCoast]
            SAT___MEAN[iFrame,:] = H.SAT___MEAN[:,iCoast]
            MODEL__STD[iFrame,:] = H.MODEL__STD[:,iCoast]
            SAT____STD[iFrame,:] = H.SAT____STD[:,iCoast]

        self.TIMES = TIMES
        self.MODEL_MEAN = MODEL_MEAN
        self.SAT___MEAN  = SAT___MEAN
        self.MODEL__STD  = MODEL__STD
        self.SAT____STD  = SAT____STD

if __name__ == "__main__":
    A = read("~/Downloads/valid.20220106.P4l.nc")
    B = read("~/Downloads/valid.20220113.P4l.nc")