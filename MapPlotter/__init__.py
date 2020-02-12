#!/usr/bin/env python
'''
Map Plotter is a toolkit that provides a framework for
2D plots of NetCDF files

	MAPPLOTTER class

	Plot NETCDF data using CARTOPY. Therefore Cartopy must be installed
	on your system.

	Example usage:

		# Define class instance
		plotter = mp.MapPlotter(projection='PlateCarree')
		params  = plotter.defaultParams() # Create basic parameters dictionary

		# To plot already loaded fields
		plotter.plot(lon,lat,data,params=params)

		# To plot data from NetCDF data
		plotter.plot_from_file(filename,varname,lonname,latname,iTime=0,iDepth=0,params=params)

		# To plot data from NetCDF data with meshmask
		plotter.plot_from_file_and_mask(filename,varname,maskfile,iTime=0,iDepth=0,masklon="glamt",masklat="gphit",params=params)

		# To see the data
		plotter.save('outfile.png',dpi=300)
		plotter.show() # wrapper of plt.show()


	Usage: map_plotter [-h] -f FILE -v VAR [-m MASK] [--lon LON] [--lat LAT]
	                   [-t TIME] [-d DEPTH] [-c CONF] -o OUT [--dpi DPI]
	
	Plot NetCDF data in beautiful 2D maps.
	
	optional arguments:
	  -h, --help               show this help message and exit
	  -f FILE, --file FILE     NetCDF file path
	  -v VAR, --var VAR        Variable to plot
	  -m MASK, --mask MASK     Mask file
	  --lon LON                Longitude variable name (default: glamt)
	  --lat LAT                Latitude variable name (default: gphit)
	  -t TIME, --time TIME     Time index for NetCDF (default: 0)
	  -d DEPTH, --depth DEPTH  Depth index for NetCDF (default: 0)
	  -c CONF, --conf CONF     Configuration file path
	  -o OUT, --outfile OUT    Output file name
	  --dpi DPI                Output file DPI (default: 300)
'''
__version__ = "0.1"

from .MapPlotter import MapPlotter