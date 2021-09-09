import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matchup.taylorDiagram import TaylorDiagram
from matchup.targetDiagram import TargetDiagram
from matplotlib.cm import ScalarMappable as ScalarMappable

class matchup(object):
    def __init__(self, Model, Ref):

        mask = ~np.isnan(Model) & ~np.isnan(Ref)

        self.Model = Model[mask]
        self.Ref   = Ref[mask]

    def diff(self):
        return self.Model - self.Ref

    def number(self):        
        return len(self.Model)

    def variances(self):
        return self.Ref.std()**2, self.Model.std()**2

    def medians(self):
        '''
        return median(Reference), median(Model)
        '''
        return np.median(self.Ref), np.median(self.Model)

    def covariance(self,output_matrix=False):
        array = (self.Model - self.Model.mean())*(self.Ref - self.Ref.mean())
        return array.mean() if not output_matrix else array

    def correlation(self,output_matrix=False):
        return self.covariance(output_matrix=output_matrix)/(self.Model.std()* self.Ref.std())

    def bias(self):
        '''
        Bias
        ID 4.3
        Calculates the mean error.
        Result of zero does not necessarily indicate low error due to cancellation.
        '''
        return self.diff().mean();

    def maxdiff(self):

        return self.diff().max()

    def MSE(self):
        '''Mean Square Error
        ID 4.4
        Calculates a mean error (in data units squared), which is not effected by cancellation.
        Squaring the data may cause bias towards large events.
        '''
        return (self.diff()**2).mean()

    def RMSE(self):
        ''' Root mean Square Error
        ID 4.5
        MSE error (4.4) except result is returned in the same units as model,
        which is useful for interpretation.
        '''
        return np.sqrt(self.MSE())

    def MAE(self):
        ''' Mean Absolute Error
        ID 4.6
        Similar to RMSE (4.5) except absolute value is used instead. This reduces
        the bias towards large events; however, it also produces a non-smooth operator
        when used in optimisation.
        '''
        return (np.abs(self.diff())).mean()

    def AME(self):
        '''Absolute Maximum Error
        ID 4.7
        Records the maximum absolute error.
        '''
        return (np.abs(self.diff())).max()



    def NSE(self):
        '''Coefficient of determination
        Nash - Sutcliff Model Efficiency
        ID 6.1

        This method compares the performance of the model to a model that
        only uses the mean of the observed data.
        A value of 1 would indicate a perfect model, while a value of zero
        indicates performance no better than simply using the mean.
        A negative value indicates even worse performance.
        '''

        return 1.0 - self.MSE()/(self.Ref.std()**2)

    def PPMC(self):
        '''
        ID 6.3
        The Pearson Product moment correlation measures the correlation
        of the measured and modelled values.
        Negatives to this model are linear model assumptions and the fact
        it can return an ideal result for a model with constant offset.
        '''
        return 1

    def RSqr(self):
        '''
        Coefficient of determination
        ID 6.4
        Squared version of 6.3, with the same interpretation of results, except range is now (0, 1).
        '''
        return 1

    def IoA(self):
        '''
        Index Of Agreement
        ID 6.5
        This method compares the sum of squared error to the potential error.
        This method is similar to 6.4 however it is designed to be better at handling
        differences in modelled and observed means and variances.
        Squared differences may add bias to large data value events (Willmott, 1981).
        '''

        array_denom = np.abs(self.Model - self.Ref.mean()) + np.abs(self.Model + self.Ref.mean())
        denom = (array_denom**2).sum()
        return 1.0 - self.number()*self.MSE()/denom

    def PI(self):
        '''Persistence Index
        ID 6.6
        The persistence index compares the sum of squared error to the error
        that would occur if the value was forecast as the previous observed value.
        Similar to 6.1 except the performance of the model is being compared to the previous value.
        '''
        s = 0
        for i in range(2,self.number()):
            s+=(self.Ref[i] - self.Ref[i-1])**2
        norm = s/self.number()
        return 1.0 - self.MSE()/norm


    def RAE(self):
        '''Relative Absolute Error
        This compares the total error relative to what the total error
        would be if the mean was used for the model.
        A lower value indicates a better performance, while a score greater
        than one indicates the model is outperformed by using the mean as the prediction.
        '''
        norm = np.abs(self.Ref - self.Ref.mean()).mean()
        return self.MAE()/norm

    def RSR(self):
        ''' RMSE - Standard Deviation  ratio
        The traditional RMSE method weighted by the standard deviation
        of the observed values (Moriasi et al., 2007)
        '''
        return self.RMSE()/self.Ref.std()

    def densityplot(self, bins, dpi=72):
        '''
        Plots the density of Model data against Reference data.

        Args:
            - *bins* : number of bins for the histogram.
            - *dpi* (optional): the figure's DPI (default: 72).

        Returns: a matplotlib Figure object and a matplotlib Axes object
        '''
        #Set the figure
        fig, ax = plt.subplots()
        fig.set_dpi(dpi)
        #Ensure that each bin is at least 10 pixels wide.
        fig.set_size_inches((bins * 10) / float(dpi), (bins * 10) / float(dpi))
        #Compute the 2D histogram
        H, xedges, yedges = np.histogram2d(self.Model, self.Ref, bins)
        #Set the axes extent
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        #Plot the 2D histogram
        im = ax.imshow(H, interpolation='nearest', extent=extent, aspect='auto', cmap="Blues")
        #Set the color bar
        div = make_axes_locatable(ax)
        cax = div.append_axes("right", size="3%", pad=0.05)
        cbar = fig.colorbar(im, cax=cax)
        return fig, ax
    
    def densityplot2(self,modelname='Model',refname='Ref',units = 'mmol m-3',sub = 'med'):
        '''
        opectool like density plot
        
        Ref is in x axis
        Model  in y axis
        
        Args
         - *modelname* (optional) string , default ='Model'
         - *refname*   (optional) string , default ='Ref'
         - *units*     (optional) string , default ='mmol m-3'
         - *sub*       (optional) string , default ='med'
        
        
        Returns: a matplotlib Figure object and a matplotlib Axes object
        '''
        
        fig, ax = plt.subplots()
        plt.title('%s Density plot of %s and %s\nNumber of considered matchups: %s' % (sub, modelname, refname, self.number()))
        cmap = 'Spectral_r'
        axis_min = min(self.Ref.min(),self.Model.min())
        axis_max = max(self.Ref.max(),self.Model.max())
        extent = [axis_min, axis_max, axis_min, axis_max]

        hexbin = ax.hexbin(self.Ref, self.Model, bins=None, extent=extent, cmap=cmap, mincnt=1)
        data = hexbin.get_array().astype(np.int32)
        MAX = data.max()

        for nticks in range(10,2,-1):
            float_array=np.linspace(0,MAX,nticks)
            int___array = float_array.astype(np.int32)
            if np.all(float_array == int___array ):
                break

        mappable = ScalarMappable(cmap=cmap)
        mappable.set_array(data)
        #fig.colorbar(mappable, ticks = int___array, ax=ax)
        cbar = fig.colorbar(mappable, ax=ax)
        labels = cbar.ax.get_yticklabels()
        FloatNumberFlag = False
        for label in labels:
            numstr = str(label.get_text())
            if numstr.find(".") > -1:
                FloatNumberFlag = True

        if FloatNumberFlag:
            cbar.remove()
            cbar = fig.colorbar(mappable, ticks = int___array, ax=ax)

        ax.set_xlabel('%s %s' % (refname,  units))
        ax.set_ylabel('%s %s' % (modelname,units))
        ax.grid()
        return fig,ax

    def targetplot(self, dpi=72):
        '''
        Plots the Taylor diagram for this matchup.

        Args:
            - *dpi* (optional): the figure's DPI (default: 72).

        Returns: a matplotlib Figure object and a matplotlib Axes object
        '''
        #Create the figure
        fig = plt.figure()
        fig.set_dpi(dpi)
        #Create the TargetDiagram object
        dia = TargetDiagram(self.Model, self.Ref, fig=fig)
        return fig, dia.ax

    def taylorplot(self, dpi=72):
        '''
        Plots the Taylor diagram for this matchup.

        Args:
            - *dpi* (optional): the figure's DPI (default: 72).

        Returns: a matplotlib Figure object and a matplotlib Axes object
        '''
        #Create the figure
        fig = plt.figure()
        fig.set_dpi(dpi)
        #Create the TaylorDiagram object
        dia = TaylorDiagram(self.Ref.std(), fig=fig, label="Reference")
        #Calculate the correlation coefficient
        corr = self.correlation()
        #Plot the model point
        dia.add_sample(self.Model.std(), corr, marker='s')
        # Add RMS contours, and label them
        contours = dia.add_contours(colors='0.5')
        plt.clabel(contours, inline=1, fontsize=10)
        return fig, dia.ax


if __name__ == "__main__" :
    n = 5600
    x = 4.0 + 1.0*np.random.randn(n)
    y = 3.0 + 0.5*np.random.randn(n)
    M = matchup(y,x)
    fig, ax = M.densityplot2(modelname='s1',refname='s2',units='Kg/m3')
    fig.show()

 
    
