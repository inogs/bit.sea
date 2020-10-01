import matplotlib.pyplot as pl
import numpy as np
from layer_integral import coastline

class figure_generator():
    def __init__(self, S):
        self.submask=S
    def _gen_structure(self, IDrun,season,subbasin_name,var_list, TEXTxlabel,xmin,xmax,xinc=None):
        '''
        Generates structure figure
        Arguments : 
        * submask * 2d logical array
        * IDrun * string 
        * season * string
        * subbasin_name * string
        
        Returns : 
         * fig *  a figure object
         * axes * a list of axis objects 
        '''
        fig, axs = pl.subplots(2,6, facecolor='w', edgecolor='k')
        axs = axs.ravel()

        x_shift = -0.025
        for ax in axs: # now shift all the axes to left
            BBox = ax.get_position()
            rect = [BBox.xmin+x_shift, BBox.ymin, BBox.width, BBox.height]
            ax.set_position(rect)

        bool2d=self.submask.mask_at_level(0)
        smaskplot = np.ones_like(bool2d,dtype=np.float32)
        smaskplot[~bool2d] = np.nan
        lon_min = self.submask.xlevels.min()
        lon_max = self.submask.xlevels.max()
        lat_min = self.submask.ylevels.min()
        lat_max = self.submask.ylevels.max()
        axs[6].imshow(smaskplot, extent=[lon_min, lon_max, lat_max, lat_min])
        axs[6].invert_yaxis()
        axs[6].set_xlim([-6, 36])
        axs[6].set_ylim([30, 46])
        c_lon, c_lat=coastline.get()

        axs[6].plot(c_lon,c_lat, color='#000000',linewidth=0.5) # can generate Warning because of nans
        axs[6].xaxis.set_tick_params(labelsize=6)
        axs[6].yaxis.set_tick_params(labelsize=6)
        axs[6].set_xticks([ 0,10,20,30])
        axs[6].set_yticks([32,36,40,44])

        BBox = axs[6].get_position() # shift only the map a bit to left
        rect = [BBox.xmin+x_shift, BBox.ymin, BBox.width, BBox.height]
        axs[6].set_position(rect)

        # alcune impostazioni #
        fig.set_dpi(200)
        hsize=9
        vsize=6
        fig.set_size_inches(hsize,vsize)
        #fig.subplots_adjust(hspace = 0.1, wspace=0.1)

#        fig.text(0.12,0.15,IDrun,fontsize=15)
        fig.text(0.12,0.40,season,fontsize=10)
        fig.text(0.12,0.35,subbasin_name,fontsize=18)
    
        for i,ax in enumerate(axs[:6]):
            ax.set_xlim(xmin[i],xmax[i])
            ax.set_ylim(-200,0)
            ax.grid(color='k',linestyle='--')

            start, end = ax.get_xlim()
            if xinc is not None:
                ax.xaxis.set_ticks(np.arange(start, end, xinc[i]))
            ax.xaxis.set_tick_params(labelsize=8)
            ax.yaxis.set_tick_params(labelsize=0)
            ax.set_title(var_list[i])
            ax.set_xlabel(TEXTxlabel[i])

        for i,ax in enumerate(axs[7:]):
            ax.set_xlim(xmin[i+1],xmax[i+1])
            ax.set_ylim(-2000,-200) # QUI CI VORREBBE UN ylim(-lastlev,-200) per ciascun subbasin
            ax.grid(color='k',linestyle='--')

            start, end = ax.get_xlim()
            if xinc is not None:
                ax.xaxis.set_ticks(np.arange(start, end, xinc[i+1]))
            #ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
            ax.xaxis.set_tick_params(labelsize=8)
            ax.set_yticks([-2000,-1500,-1000,-500])
            ax.yaxis.set_tick_params(labelsize=0)
        
        axs[0].set_yticks([-200,-150,-100,-50,0])
        axs[0].yaxis.set_tick_params(labelsize=8)
        axs[7].set_yticks([-2000,-1500,-1000,-500])
        axs[7].yaxis.set_tick_params(labelsize=8)
        
        return fig, axs

    def gen_structure_1(self,IDrun,season,subbasin_name):
        var_list=['CHL','NO$_3$','PO$_4$','O$_2$','SiO$_2$','NH$_4$']
        TEXTxlabel=['mgChl/m'+u'\u00B3','mmolN/m'+u'\u00B3','mmolP/m'+u'\u00B3','mmolO$_2$'+'/m'+u'\u00B3','mmolSi/m'+u'\u00B3','mmolNH$_4$/m'+u'\u00B3']
        xmin=[0, 0,0,  180,0,0]
        xmax=[1,10,0.6,280,9,1]
        return self._gen_structure(IDrun, season, subbasin_name, var_list, TEXTxlabel, xmin, xmax)

    def gen_structure_2(self,IDrun,season,subbasin_name):
        var_list=['CHL','PAR','Pc','PP','Bc'] 
        TEXTxlabel=['mgChl/m'+u'\u00B3','xPAR','mgC/m'+u'\u00B3','mgC'+'/m'+u'\u00B3'+'/day','mgC/m'+u'\u00B3']
        xmin=[0, 0,0,  -5, 0]
        xmax=[1,10,25, 70,15]
        return self._gen_structure(IDrun, season, subbasin_name, var_list, TEXTxlabel, xmin, xmax)


    def gen_structure_3(self,IDrun,season,subbasin_name):
        var_list=['PCO$_2$','DIC','ALK','pH_T','R2c'] 
        TEXTxlabel=['ppm',r'$\mu$'+'mol/kg',r'$\mu$'+'mol/kg','','mgC/m'+u'\u00B3']
        xmin=[300,2050,2350,7.9, 0]
        xmax=[480,2400,2750,8.2,50]
        xinc=[ 40, 100, 100,0.1,10]
        return self._gen_structure(IDrun, season, subbasin_name, var_list, TEXTxlabel, xmin, xmax, xinc)
    
        
def add_legend(ax):
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,labelspacing=-0.25, handletextpad=0,borderpad=0.1)
        leg = ax.get_legend()
        ltext  = leg.get_texts()
        pl.setp(ltext,fontsize=8)


def profile_plotter(depth,values,color,ax_top,ax_bottom,label):
    '''
    Plots a model profile over two axis objects.

    Arguments :
    * depth     * array of depths
    * values    * numpy array
    * color     * string or rgba object
    * ax_top    * axis object
    * ax_bottom * axis object
    * label     * string
    '''

    ax_top.plot(values,depth,color=color,label=label)
    if ax_bottom is not None:
        ax_bottom.plot(values,depth,color=color)


def clim_profile_plotter(depth,mean,std,ax_top,ax_bottom):
    '''
    Plots a model profile over two axis objects.

    Arguments :
    * depth   * numpy array of depths
    * mean    * numpy array
    * std     * numpy array
    * ax_top    * axis object
    * ax_bottom * axis object
    '''
    bool_top = depth>= -200
    bool_bot = ~bool_top
    ax_top.plot(mean,depth,'ro')
    ax_top.plot(mean+std,depth,'r-.')
    ax_top.plot(mean-std,depth,'r-.')
    if ax_bottom is not None:
        ax_bottom.plot(mean[bool_bot],depth[bool_bot],'ro')
        ax_bottom.plot(mean[bool_bot]+std[bool_bot],depth[bool_bot],'r-.')
        ax_bottom.plot(mean[bool_bot]-std[bool_bot],depth[bool_bot],'r-.')



if __name__=="__main__":
    from basins import V2 as basV2
    from commons.mask import Mask
    from commons.submask import SubMask

    maskfile="/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/V1/meshmask_872.nc"
    #maskfile="/Users/gbolzon/Documents/workspace/ogs_bounday_conditions/masks/meshmask_872.nc"
    TheMask= Mask(maskfile)
    S = SubMask(basV2.lev1, maskobject=TheMask)
    
    pl.close('all')
    F = figure_generator(S)
    fig,axes = F.gen_structure_1('pre-OP','winter','lev1')
    
    x=np.arange(72)
    depth   = np.array([12, 60, 80, 150, 200, 250, 500, 1500])

    # SYNTHETIC PROFILES #####################
    Chl = 0.05 + 0.1 * np.exp(-((x-2)/4)**2) 
    N3n = 2    + 0.3 * np.exp(-((x-5)/8)**2) 
    N1p = 0.05 + 0.3 * np.exp(-((x-20)/4)**2)
    O2o = 200  + 0.5 * np.exp(-((x-60)/4)**2)
    Sil = 3    + 0.3 * np.exp(-((x-4)/4)**2)

    par = 3    + 0.3 * np.exp(-((x-5)/4)**2)
    pc  = 0.1  + 0.3 * np.exp(-((x-1)/4)**2)
    pp  = 200  + 0.5 * np.exp(-((x-60)/4)**2)
    bc  = 3    + 0.3 * np.exp(-((x-4)/4)**2) 
    pco = 350  + 0.1 * np.exp(-((x-2)/4)**2)
    dic = 2100 + 3 * np.exp(-((x-5)/4)**2) 
    alk = 2200 + 2 * np.exp(-((x-20)/4)**2) 
    pht = 7.95 + 0.5 * np.exp(-((x-1)/4)**2) 
    r2c = 0.03 + 0.3 * np.exp(-((x-4)/4)**2) 
    #----------------------------------------
    
    profile_plotter(-TheMask.zlevels, Chl, '0.8', axes[0], None,'2014')
    profile_plotter(-TheMask.zlevels, N3n, '0.8', axes[1], axes[6],'2014')
    profile_plotter(-TheMask.zlevels, N1p, '0.8', axes[2], axes[7],'2014')
    profile_plotter(-TheMask.zlevels, O2o, '0.8', axes[3], axes[8],'2014')
    profile_plotter(-TheMask.zlevels, Sil, '0.8', axes[4], axes[9],'2014')
    profile_plotter(-TheMask.zlevels, 2*Chl,'0.4', axes[0], None,'2015')
    profile_plotter(-TheMask.zlevels, 2*N3n,'0.4', axes[1], axes[6],'2015')
    profile_plotter(-TheMask.zlevels, 2*N1p,'0.4', axes[2], axes[7],'2015')
    profile_plotter(-TheMask.zlevels, 1.2*O2o,'0.4', axes[3], axes[8],'2015')
    profile_plotter(-TheMask.zlevels, 2*Sil,'0.4', axes[4], axes[9],'2015')


    MeanClim=np.array([0.04, 0.08, 0.3, 0.35, 0.4, 0.4, 0.5, 0.5])
    std=np.ones_like(MeanClim)*0.05
    clim_profile_plotter(-depth, MeanClim, std , axes[2], axes[7])

    add_legend(axes)
    #fig.show()
    fig.savefig('test1.png', dpi=150)
    

    fig,axes = F.gen_structure_2('pre-OP','winter','lev1')
    profile_plotter(-TheMask.zlevels, Chl, '0.8', axes[0], None,'2014')
    profile_plotter(-TheMask.zlevels, par, '0.8', axes[1], axes[6],'2014')
    profile_plotter(-TheMask.zlevels, pc , '0.8', axes[2], axes[7],'2014')
    profile_plotter(-TheMask.zlevels, pp , '0.8', axes[3], axes[8],'2014')
    profile_plotter(-TheMask.zlevels, bc , '0.8', axes[4], axes[9],'2014')
    profile_plotter(-TheMask.zlevels, 2*Chl, '0.4', axes[0], None,'2015')
    profile_plotter(-TheMask.zlevels, 2*par, '0.4', axes[1], axes[6],'2015')
    profile_plotter(-TheMask.zlevels, 2*pc,  '0.4', axes[2], axes[7],'2015')
    profile_plotter(-TheMask.zlevels, 1.2*pp,'0.4', axes[3], axes[8],'2015')
    profile_plotter(-TheMask.zlevels, 2*bc,  '0.4', axes[4], axes[9],'2015')

    MeanClim1=np.array([0.04, 0.08, 0.3, 0.35, 0.4, 0.4, 0.5, 0.5])
    MeanClim2=np.array([0.04, 0.08, 0.3, 0.35, 0.4, 0.4, 0.5, 0.5])
    std1=np.ones_like(MeanClim1)*0.05
    std2=np.ones_like(MeanClim2)*0.5
    clim_profile_plotter(-depth, MeanClim1, std1 , axes[2], axes[7])
    clim_profile_plotter(-depth, MeanClim2, std2 , axes[3], axes[8])
    add_legend(axes)
    fig.savefig('test2.png', dpi=150)

    fig,axes = F.gen_structure_3('pre-OP','winter','lev1')    
    profile_plotter(-TheMask.zlevels, pco, '0.8', axes[0], None,'2014')
    profile_plotter(-TheMask.zlevels, dic, '0.8', axes[1], axes[6],'2014')
    profile_plotter(-TheMask.zlevels, alk, '0.8', axes[2], axes[7],'2014')
    profile_plotter(-TheMask.zlevels, pht, '0.8', axes[3], axes[8],'2014')
    profile_plotter(-TheMask.zlevels, r2c, '0.8', axes[4], axes[9],'2014')
    profile_plotter(-TheMask.zlevels, 2*pco,'0.4', axes[0], None,'2015')
    profile_plotter(-TheMask.zlevels, 2*dic,'0.4', axes[1], axes[6],'2015')
    profile_plotter(-TheMask.zlevels, 2*alk,'0.4', axes[2], axes[7],'2015')
    profile_plotter(-TheMask.zlevels, 1.2*pht,'0.4', axes[3], axes[8],'2015')
    profile_plotter(-TheMask.zlevels, 2*r2c,'0.4', axes[4], axes[9],'2015')

    MeanClim=np.array([0.04, 0.08, 0.3, 0.035, 0.04, 0.04, 0.05, 0.05])
    std=np.ones_like(MeanClim)*0.05
    clim_profile_plotter(-depth, MeanClim, std , axes[4], axes[9])
    add_legend(axes)
    fig.savefig('test3.png', dpi=150)   
    
