import numpy as np
from commons.layer import Layer
from dataextractor import DataExtractor

class MapBuilder(object):

    @staticmethod
    def get_layer_average(data_extractor, layer, timestep=0):
        """Returns a 2D NumPy array with the average weighted over depth.

        Args:
            - *data_extractor*: a DataExtractor object that contains the data to
              extract.
            - *layer*: a Layer object that contains the vertical boundaries for
              the computation.
            - *timestep* (optional): the index of the first dimension (time) in
              the model data. Default: 0.
        """
        if not isinstance(layer, (Layer,)):
            raise ValueError("layer must be a Layer object")
        if not isinstance(data_extractor, (DataExtractor,)):
            raise ValueError("dataextractor must be a DataExtractor object")
        #Find Z indices
        top_index = np.where(data_extractor._mask.zlevels >= layer.top)[0][0]
        bottom_index = np.where(data_extractor._mask.zlevels < layer.bottom)[0][-1]
        if top_index == bottom_index:
            #Just one layer so we return the sliced data
            output = data_extractor.get_filled_values[timestep,top_index,:,:]
            return output
        #Workaround for Python ranges
        bottom_index += 1
        #Build local mask matrix
        lmask = np.array(data_extractor._mask.mask[top_index:bottom_index,:,:], dtype=np.double)
        #Build dz matrix
        dzm = np.ones_like(lmask, dtype=np.double)
        j = 0
        for i in range(top_index, bottom_index):
            dzm[j,:,:] = data_extractor._mask.dz[i]
            j += 1
        #Get the slice of the values
        v = np.array(data_extractor.values[timestep,top_index:bottom_index,:,:])
        #Build integral matrix (2D)
        integral = (v * dzm * lmask).sum(axis=0)
        #Build height matrix (2D)
        height = (dzm * lmask).sum(axis=0)
        indexmask = [height > 0]
        #Build output matrix (2D)
        output = np.full_like(integral, data_extractor.fill_value, dtype=np.double)
        output[indexmask] = integral[indexmask] / height[indexmask]
        return output
