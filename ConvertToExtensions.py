import FittingUtilities
from astropy.io import fits as pyfits
import sys
import os
import numpy
import matplotlib.pyplot as plt
import HelperFunctions

left_trim = 5
right_trim = 0
bad_regions = {}


if __name__ == "__main__":
  fileList = []
  for arg in sys.argv[1:]:
    fileList.append(arg)
    
  for fname in fileList:
    outfilename = "%s-0.fits" %(fname.split(".fits")[0])
    header = pyfits.getheader(fname)

    try:
      orders = HelperFunctions.ReadFits(fname)
    except ValueError:
      orders = HelperFunctions.ReadFits(fname, errors=2)
    orders = orders[::-1]    #Reverse order so the bluest order is first
    
    column_list = []
    for i, order in enumerate(orders):
      
      left, right = left_trim, order.size()-right_trim
      if i in bad_regions.keys():
        region = bad_regions[i]
        left = numpy.searchsorted(order.x, region[0])
        right = numpy.searchsorted(order.x, region[1])
        if left == 0 or right == order.size():
          order.x = numpy.delete(order.x, numpy.arange(left, right))
          order.y = numpy.delete(order.y, numpy.arange(left, right))
          order.cont = numpy.delete(order.cont, numpy.arange(left, right))
          order.err = numpy.delete(order.err, numpy.arange(left, right))
          
        else:
          print "Warning! Bad region covers the middle of order %i" %i
          print "Interpolating rather than removing"
          order.y[left:right] = order.cont[left:right]
          order.err[left:right] = 9e9
      else:
        order = order[left:right]
      if order.size() < 10:
        continue
      


      order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=1.5, highreject=10)
      columns = columns = {"wavelength": order.x,
                           "flux": order.y,
                           "continuum": order.cont,
                           "error": order.err}
      column_list.append(columns)
    HelperFunctions.OutputFitsFileExtensions(column_list, fname, outfilename, mode="new")
      
      
