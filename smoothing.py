#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 25 20:47:21 2024

@author: pschlatt
"""

import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

#%%
%matplotlib qt

#%%

# Initial x and y arrays
x = np.linspace(0, 10, 30)
y = np.sin(0.5*x)*np.sin(x*np.random.randn(30))
# Spline interpolation
spline = UnivariateSpline(x, y, s = 6)
x_spline = np.linspace(0, 10, 1000)
y_spline = spline(x_spline)
# Plotting
fig = plt.figure()
plt.subplots_adjust(bottom=0.25)
ax = fig.subplots()
p = ax.plot(x,y)
p, = ax.plot(x_spline, y_spline, 'g')
# Defining the Slider button
# xposition, yposition, width and height
ax_slide = plt.axes([0.25, 0.1, 0.65, 0.03])
# Properties of the slider
s_factor = Slider(ax_slide, 'Smoothing factor',
                  0.1, 6, valinit=6, valstep=0.2)
# Updating the plot
def update(val):
    current_v = s_factor.val
    spline = UnivariateSpline(x, y, s = current_v)
    p.set_ydata(spline(x_spline))
    #redrawing the figure
    fig.canvas.draw()
    
# Calling the function "update" when the value of the slider is changed
s_factor.on_changed(update)
plt.show()

#%%