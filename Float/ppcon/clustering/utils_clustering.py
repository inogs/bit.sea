import os
import pandas as pd

def make_ds(training_folder,  INDIR):
    if training_folder == "SUPERFLOAT":
        from make_ds_clustering import make_pandas_df

        if not os.path.exists( INDIR  ):
            os.mkdir( INDIR  )

        if not os.path.exists( INDIR   + "/clustering/"):
            os.mkdir(          INDIR   + "/clustering/")

        if not os.path.exists( INDIR   + "/clustering/ds_sf_clustering.csv"):
            print("making ds...")
            make_pandas_df(    INDIR   + '/SUPERFLOAT/Float_Index.txt' , INDIR   + "/clustering/ds_sf_clustering.csv", INDIR)
            print("superfloat clustering complete dataset created")

    return
