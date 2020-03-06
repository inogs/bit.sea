#!/usr/bin/env python
#
# Map Plotter Class
#
# Plot NETCDF Data using CARTOPY
#
# Arnau Miro, OGS (2019)
from __future__ import print_function

import io, requests, json
import numpy as np, matplotlib, matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import netCDF4 as NC

class MapPlotter():
	'''
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

	'''
	def __init__(self,projection='PlateCarree',**kwargs):
		'''
		Class constructor. Define a MapPlotter instance.

		Inputs:
			> Projection: Type of projection according to cartopy
						  as a string 
				https://scitools.org.uk/cartopy/docs/latest/crs/projections.html

			Additional arguments for the projection can be added
			as kwargs to this function.

		Output:
			> MapPlotter class instance
		'''
		self._projection = getattr(ccrs,projection)(**kwargs)
		# Class variables
		self._fig = None
		self._ax  = None
		self._gl  = None
		self._plt = None

	@property
	def projection(self):
		return self._projection
	@projection.setter
	def projection(self,value):
		self._projection = getattr(ccrs,value)

	@staticmethod
	def loadNC(fname,varname,fill_value=1e20,mask_value=np.nan):
		'''
		Load data from NETCDF file. This function is identical to
		bit.sea commons.netcdf4 readfile. It will also filter any
		fill value with a mask value.

		Inputs:
			> fname:      NetCDF file name
			> varname:    Variable name inside NetCDF file
			> fill_value: NetCDF fill value (default: 1e20)
			> mask_value: Value to replace the fill value (default: NaN)

		Output:
			> variable as a numpy array.
		'''
		D = NC.Dataset(fname,'r')
		v = np.array(D.variables[varname])
		D.close()
		v[v >= fill_value] = mask_value
		return v

	@staticmethod
	def defaultParams():
		'''
		Return default plot params dictionary.

		Outputs:
			> params dictionary
		'''
		params = {
			# Figure and axes params
			'fig'         : None,
			'ax'          : None,
			# Figure params
			'size'        : (8,6),
			'dpi'         : 100,
			'style'       : 'ggplot',
			# Axes params
			'xlim'        : [-180,180],
			'ylim'        : [-90,90],
			'max_div'     : 4,
			'top_label'   : False,
			'bottom_label': True,
			'right_label' : False,
			'left_label'  : True,
			'grdstyle'    : {},
			'grdargs'     : {'draw_labels':True,'linewidth':0},
			'features'    : ['coastline','continents','rivers','image'],
			'res'         : '50m',
			'img'         : None,
			# Title and labels
			'title'       : [], # [title,kwargs]
			'xlabel'      : [],
			'ylabel'      : [],
			# Plot params
			'alpha'       : 1.,
			# Colormap and colorbar
			'cmap'        : 'coolwarm',
			'ncol'        : 256, 
			'draw_cbar'   : True,
			'orientation' : 'horizontal',
			'bounds'      : [-1e30, 1e30], # [min,max]
			'numticks'    : 10,
			'tick_format' : '%.2f',
			'tick_font'   : None,
			'label'       : {'label':'','weight':None,'style':None}
		}
		return params

	@staticmethod
	def loadParams(filename):
		'''
		Load params dictionary from file
		'''
		file   = open(filename,'r')
		params = json.load(file)
		file.close() 
		return params

	@staticmethod
	def saveParams(filename,params):
		'''
		Save params dictionary to file
		'''
		file = open(filename,'w')
		json.dump(params,file)
		file.close()

	@staticmethod
	def show():
		'''
		Wrapper to matplotlib.pyplot.show()
		'''
		plt.show()

	def save(self,filename,dpi=300):
		'''
		Save figure to disk.

		Inputs:
			> filename: Output file name
			> dpi:      Pixel depth (default: 300)
		'''
		if self._fig and self._ax and self._plot:
			self._fig.savefig(filename,dpi=dpi,bbox_inches='tight')

	def createFigure(self,sz=(8,6),dpi=100,style=None):
		'''
		Create a new figure.

		Inputs:
			> sz:    Figure size (default: (8,6))
			> dpi:   Pixel depth (default: 100)
			> style: Name of a matplotlib style sheet
				https://matplotlib.org/3.1.1/gallery/style_sheets/style_sheets_reference.html

		Outputs:
			> Figure object
		'''
		fig = None
		if style:
			plt.style.use(style)
			fig = plt.figure(figsize=sz,dpi=dpi)
		else:
			fig = plt.figure(figsize=sz,dpi=dpi,facecolor='w',edgecolor='k')
		return fig

	def createAxes(self):
		'''
		Create new axes.

		Outputs:
			> Axes object
		'''
		return plt.axes(projection=self._projection)

	def createGridlines(self,xlim=[-180,180],ylim=[-90,90],top=False,bottom=True,left=True,right=False,max_div=4,style={},gridlines_kwargs={'draw_labels':True,'linewidth':0}):
		'''
		Create gridlines for the current axes.

		Inputs:
			> xlim:             Minimum and maximum extend for the x axis (default: [-180,180])
			> ylim:             Minimum and maximum extend for the y axis (default: [-90,90])
			> top:              Display labels at the top (default: False)
			> right:            Display labels at the right (default: False)
			> max_div:          Maximum number of divisions on the axes (defaut: False).
			> style:            Style dictionary for the axis labels.
			> gridlines_kwargs: Extra gridlines arguments dictionary.

		Outputs:
			> Gridlines object
		'''
		gl = None
		if self._ax:
			# Axes limits
			self._ax.set_xlim(xlim)
			self._ax.set_ylim(ylim)
			# Grid lines
			gl = self._ax.gridlines(crs=self._projection,**gridlines_kwargs)
			gl.xlocator       = matplotlib.ticker.FixedLocator(np.arange(xlim[0],xlim[1],max_div))
			gl.ylocator       = matplotlib.ticker.FixedLocator(np.arange(ylim[0],ylim[1],max_div))
			gl.xformatter     = LONGITUDE_FORMATTER
			gl.yformatter     = LATITUDE_FORMATTER
			gl.xlabels_top    = top
			gl.xlabels_bottom = bottom
			gl.xlabel_style   = style
			gl.ylabels_left   = left
			gl.ylabels_right  = right
			gl.ylabel_style   = style
		return gl

	def setTitle(self,title,**kwargs):
		'''
		Set the title on the current axes.

		Inputs:
			> title:  Title string.

		Any kwargs to format the title are also accepted by this function.
		'''
		if self._ax:
			self._ax.set_title(title,**kwargs)

	def setXAxis(self,label,**kwargs):
		'''
		Set the x axis label on the current axes.

		Inputs:
			> label:  Label string.

		Any kwargs to format the label are also accepted by this function.
		'''
		if self._ax:
			self._ax.set_xlabel(label,**kwargs)

	def setYAxis(self,label,**kwargs):
		'''
		Set the y axis label on the current axes.

		Inputs:
			> label:  Label string.

		Any kwargs to format the label are also accepted by this function.
		'''
		if self._ax:
			self._ax.set_ylabel(label,**kwargs)

	def drawCoastline(self,res='50m',color=(0,0,0,1),linewidth=.75):
		'''
		Draws the coastline using Natural Earth Features.

		Inputs:
			> res:       Resolution (default: 50m)
			> color:     Color of the coastline (default (0,0,0,1))
			> linewidth: Width of the coastline (default .75)
		'''
		if self._ax:
			self._ax.add_feature(
				cfeature.NaturalEarthFeature('physical', 'coastline', res, 
					edgecolor=color, facecolor='none', linewidth=linewidth)
				)

	def drawContinentBorders(self,res='50m',color=(0,0,0,.6),linewidth=.75):
		'''
		Draws the continent borders using Natural Earth Features.

		Inputs:
			> res:       Resolution (default: 50m)
			> color:     Color of the border line (default (0,0,0,.6))
			> linewidth: Width of the border line (default .75)
		'''
		if self._ax:
			self._ax.add_feature(
				cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', res, 
					edgecolor=color, facecolor='none', linewidth=linewidth)
				)

	def drawRivers(self,res='50m',color=(0,0,.545,.75),linewidth=.5):
		'''
		Draws the rivers using Natural Earth Features.

		Inputs:
			> res:       Resolution (default: 50m)
			> color:     Color of the rivers (default (0,0,.545,.75))
			> linewidth: Width of the rivers (default .5)
		'''
		if self._ax:
			self._ax.add_feature(
				cfeature.NaturalEarthFeature('physical', 'rivers_lake_centerlines', res, 
					edgecolor=color, facecolor='none', linewidth=linewidth)
				)

	def drawBackground(self,img=None):
		'''
		Draw a background image. If no image is provided it uses cartopy's stock image.

		Inputs:
			> img: URL or full path to a high resolution image.
		'''
		if self._ax:
			if not img == None or img == '':
				# Detect if we are dealing with a URL or a path
				if 'https://' in img or 'http://' in img:
					self._ax.imshow(plt.imread(io.BytesIO(requests.get(img).content)), origin='upper', 
						transform=self._projection, extent=[-180, 180, -90, 90])
				else:
					self._ax.imshow(plt.imread(img), origin='upper', 
						transform=self._projection, extent=[-180, 180, -90, 90])
			else:
				self._ax.stock_img()

	def setColormap(self,cmap='coolwarm',ncol=256):
		'''
		Set the colormap.

		Inputs:
			> cmap: Colormap name (default: 'coolwarm')
				https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html
			> ncol: Number of colors on the colormap (default: 256)

		Outputs:
			> Colormap object
		'''
		return plt.get_cmap(cmap,ncol)

	def setColorbar(self,orientation='horizontal',extend='neither',numticks=10,tick_format='%.2f',tick_font=None,
		label={}):
		'''
		Set the colorbar.

		Inputs:
			> orientation: Orientation of the colorbar (default: 'horizontal')
			> extend:      Extend of the colorbar (default: 'neither')
			> numticks:    Number of ticks on the colorbar (default: 10)
			> tick_format: Format of the colorbar ticks (default: "%.2f")
			> tick_font:   Tick font size
			> label:       Label dictionary for the colorbar

		Outputs:
			> Colorbar object
		'''
		cbar = None
		if self._fig and self._plot:
			cbar = self._fig.colorbar(self._plot,orientation=orientation,extend=extend)
			cbar.set_label(**label)
			if tick_font: cbar.ax.tick_params(labelsize=tick_font)
			cbar.locator   = matplotlib.ticker.LinearLocator(numticks=numticks)
			cbar.formatter = matplotlib.ticker.FormatStrFormatter(tick_format)
			cbar.update_ticks()
		return cbar

	def plot_empty(self,params=None):
		'''
		Plot the map with all the settings and without data.

		Outputs:
			> Figure object
		'''
		if params == None:
			params = self.defaultParams()
		# Do we have figure and axes?
		if not params['fig'] and not params['ax']:
			self._fig = self.createFigure(sz=params['size'],dpi=params['dpi'],style=params['style'])
			self._ax  = self.createAxes()
		if params['fig'] and not params['ax']:
			self._fig = params['fig']
			self._ax  = self.createAxes()
		if params['fig'] and params['ax']:
			self._fig = params['fig']
			self._ax  = params['ax']

		# Axes grid lines
		self._gl = self.createGridlines(xlim=params['xlim'],ylim=params['ylim'],top=params['top_label'],
			bottom=params['bottom_label'],right=params['right_label'],left=params['left_label'],
			max_div=params['max_div'],style=params['grdstyle'],gridlines_kwargs=params['grdargs'])

		# Title and axis labels
		if not params['title'] == []:
			if len(params['title']) == 1:
				self.setTitle(params['title'][0])
			else:
				self.setTitle(params['title'][0],**params['title'][1])

		if not params['xlabel'] == []:
			if len(params['xlabel']) == 1:
				self.setXAxis(params['xlabel'][0])
			else:
				self.setXAxis(params['xlabel'][0],**params['xlabel'][1])

		if not params['ylabel'] == []:
			if len(params['ylabel']) == 1:
				self.setYAxis(params['ylabel'][0])
			else:
				self.setYAxis(params['ylabel'][0],**params['ylabel'][1])

		# Add features
		if 'coastline' in params['features']:
			self.drawCoastline(res=params['res'])
		if 'continents' in params['features']:
			self.drawContinentBorders(res=params['res'])
		if 'rivers' in params['features']:
			self.drawRivers(res=params['res'])
		if 'image' in params['features']:
			self.drawBackground(img=params['img'])

		return self._fig

	def plot(self,lon,lat,data,params=None):
		'''
		Main plotting function. Plots given the longitude, latitude and data.
		An optional params dictionary can be inputted to control the plot.

		Inputs:
			> lon:    Longitude vector or matrix
			> lat:    Latitude vector or matrix
			> data:   Data matrix
			> params: Optional parameter dictionary

		Outputs:
			> Figure object
		'''
		self.plot_empty(params=params)

		# Set maximum and minimum
		z_min, z_max = np.nanmin(data), np.nanmax(data)
		
		cbar_min = params['bounds'][0] if params['bounds'][0] >= -1e20 else z_min
		cbar_max = params['bounds'][1] if params['bounds'][1] <= 1e20  else z_max

		# Set extend
		extend = 'neither'
		if (cbar_min > z_min): extend = 'min'
		if (cbar_max < z_max): extend = 'max'
		if (cbar_min > z_min and cbar_max < z_max): extend = 'both'

		# Plot
		self._plot = self._ax.pcolormesh(lon,lat,data,
					cmap=self.setColormap(cmap=params['cmap'],ncol=params['ncol']),
					norm=matplotlib.colors.Normalize(cbar_min,cbar_max),
					alpha=params['alpha'],
					transform=self._projection
					 		 )

		# Colorbar
		if params['draw_cbar']:
			self.setColorbar(orientation=params['orientation'],
							 extend=extend,
							 numticks=params['numticks'],
							 tick_format=params['tick_format'],
							 tick_font=params['tick_font'],
							 label=params['label']
							)

		return self._fig

	def plot_from_file(self,filename,varname,lonname,latname,iTime=0,iDepth=0,params=None):
		'''
		Plot function. Plots data given a NetCDF file and the names of the variables
		as well as the current depth and time index.

		Inputs:
			> filename: NetCDF file path
			> varname:  Name of the NetCDF variable to plot
			> lonname:  Name of the longitude dimension
			> latname:  Name of the latitude dimension
			> iTime:    Time index for NetCDF (default: 0)
			> iDepth:   Depth index for NetCDF (default: 0)
			> params: Optional parameter dictionary

		Outputs:
			> Figure object
		'''
		# Load longitude and latitude
		lon  = self.loadNC(filename,lonname)
		lat  = self.loadNC(filename,latname)
		# Load data
		data = self.loadNC(filename,varname)
		if len(data.shape) == 4:
			data = data[iTime,iDepth,:,:]
		else:
			data = data[iDepth,:,:]
		# Plot
		return self.plot(lon,lat,data,params=params)

	def plot_from_file_and_mask(self,filename,varname,maskfile,iTime=0,iDepth=0,
		masklon='glamt',masklat='gphit',params=None):
		'''
		Plot function. Plots data given a NetCDF file, a mask file, the names of 
		the variables as well as the current depth and time index.

		Inputs:
			> filename: NetCDF file path
			> varname:  Name of the NetCDF variable to plot
			> maskfile: Path to the mask file
			> iTime:    Time index for NetCDF (default: 0)
			> iDepth:   Depth index for NetCDF (default: 0)
			> masklon:  Name of the longitude dimension (default: 'glamt')
			> masklat:  Name of the latitude dimension (default: 'gphit')
			> params: Optional parameter dictionary

		Outputs:
			> Figure object		
		'''
		# Load longitude and latitude
		lon  = self.loadNC(maskfile,masklon)
		if len(lon.shape) == 4:
			lon = lon[0,0,:,:]
		else:
			lon = lon[0,:,:]
		lat  = self.loadNC(maskfile,masklat)
		if len(lat.shape) == 4:
			lat = lat[0,0,:,:]
		else:
			lat = lat[0,:,:]
		# Load data
		data = self.loadNC(filename,varname)
		if len(data.shape) == 4:
			data = data[iTime,iDepth,:,:]
		else:
			data = data[iDepth,:,:]
		# Plot
		return self.plot(lon,lat,data,params=params)