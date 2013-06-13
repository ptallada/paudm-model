# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 08:39:31 2013

@author: dpiscia
"""
import numpy
import matplotlib.pyplot as pyplot 

def plot_xy(x_data,y_data, xlabel = None, ylabel = None,  filename = None, title="Title"):
       '''
       plot xy with colorbar
       '''
       
       heatmap, xedges, yedges = numpy.histogram2d(x_data, y_data, bins=25)
       extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
       pyplot.clf()
       pyplot.imshow(heatmap, extent = extent, aspect='auto')
       pyplot.colorbar()
       pyplot.grid()
       pyplot.xlabel(xlabel)
       pyplot.ylabel(ylabel)
       pyplot.title(title)

        # Show or save
       if filename == None:
            pyplot.show()
       else:
            pyplot.savefig(filename)
            
def plot_hist(x_data, filename = None, xlabel = None, ylabel = None, title="Title"):
    pyplot.clf()
    pyplot.hist(x_data)
    pyplot.title(title)
    pyplot.xlabel(xlabel)
    pyplot.ylabel(ylabel)
    if filename == None:
            pyplot.show()
    else:
            pyplot.savefig(filename)
    