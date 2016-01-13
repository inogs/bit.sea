# Copyright (c) 2016 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>

import numpy as np
import matplotlib.pyplot as plt

class TargetDiagram(object):
    """
    Class for creating Target Diagrams, see:
    Jolliff, Jason K., et al.
    "Summary diagrams for coupled hydrodynamic-ecosystem model skill assessment."
    Journal of Marine Systems 76.1 (2009): 64-82.
    """

    def __init__(self, model, reference, fig=None):
        """
        targetDiagram constructor.

        Args:
            - *model*: Numpy array of model data.
            - *reference*: Numpy array of reference (observation) data.
            - *fig* (optional): matplotlib figure object. If None a new
              one will be created (default: None).
        """
        if (fig is None):
            # Create the figure
            self.fig = plt.figure()
        else:
            self.fig = fig

        # Fetch default axes for Target Diagram plot
        self.ax = self.fig.gca()
        # Set limits
        self.ax.set_xlim((-1.3, 1.3))
        self.ax.set_ylim((-1.3, 1.3))
        # Set equal aspect for X and Y
        self.ax.set_aspect('equal')
        # Display grid
        self.ax.grid(b=True)
        # Set spines (the X and Y axes of the resulting image)
        self.ax.spines['left'].set_position('center')
        self.ax.spines['right'].set_position('center')
        self.ax.spines['top'].set_position('center')
        self.ax.spines['bottom'].set_position('center')
        # Draw the 1.0 distance circle
        cir = plt.Circle((0,0), radius=1.0, fill=None, linestyle='dashed')
        self.ax.add_patch(cir)

        # Set up colors array and color index
        self._ca = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
        self._ci = 0

        # Set up shapes array and shape index
        self._sa = ['s', 'o', '^', '*']
        self._si = 0

        # Plot first point
        self.add_sample(model, reference)

    def __get_shape(self):
        return self._sa[self._si]

    def __get_color(self):
        # Cycle through colors array
        c = self._ca[self._ci]
        self._ci = (self._ci + 1) % len(self._ca)
        # If we cycled through all the colors then change shape
        if self._ci == 0:
            self._si = (self._si + 1) % len(self._sa)
        return c

    def add_sample(self, model, reference, label=None):
        """
        Adds a new point to the diagram.

        Args:
            - *model*: Numpy array of model data.
            - *reference*: Numpy array of reference (observation) data.
            - *label* (optional): The label for this sample (default: None).
        """
        # Calculate means
        m_mean = np.mean(model)
        r_mean = np.mean(reference)
        # Calculate sigmas
        sigma_r = reference.std()
        sigma_m = model.std()
        sigma_n = sigma_m // sigma_r
        # Calculate correlation factor (R)
        R = np.corrcoef(model, reference)[1,0]
        # Calculate normalized bias
        n_bias = (m_mean - r_mean) // sigma_r
        # Calculate the Normalized Unbiased Root Mean Square Error
        nurmse = np.sqrt(1 + sigma_n**2 - 2 * sigma_n * R)
        # Calculate sample coordinates
        x = nurmse * np.sign(sigma_m - sigma_r)
        y = n_bias
        # Plot the sample
        color = self.__get_color()
        marker = self.__get_shape()
        self.ax.plot(x,y, marker=marker, color=color, label=label)
