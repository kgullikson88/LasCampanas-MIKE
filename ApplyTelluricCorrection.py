from astropy.io import fits as pyfits
import sys
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import matplotlib.pyplot as plt
import DataStructures
import os
import FittingUtilities
import HelperFunctions
import numpy

plot = True
plotorder = 29

def ReadCorrectedFile(fname, yaxis="model"):
  orders = []
  headers = []
  hdulist = pyfits.open(fname)
  numorders = len(hdulist)
  for i in range(1, numorders):
    order = hdulist[i].data
    xypt = DataStructures.xypoint(x=order.field("wavelength"),
                                  y=order.field(yaxis),
                                  cont=order.field("continuum"),
                                  err=order.field("error"))

    orders.append(xypt)
    headers.append(hdulist[i].header)
  return orders, headers


def Correct(original, corrected, offset=None):
  #Read in the data and model
  original_orders = HelperFunctions.ReadFits(original, extensions=True, x="wavelength", y="flux", errors="error", cont="continuum")
  corrected_orders, corrected_headers = ReadCorrectedFile(corrected)
  test_orders, header = ReadCorrectedFile(corrected, yaxis="flux")

  if plot:
    order = test_orders[plotorder]
    model = corrected_orders[plotorder]
    #for order, model in zip(test_orders, corrected_orders):
    plt.plot(order.x, order.y/order.cont)
    plt.plot(model.x, model.y)
    plt.title("Correction in corrected file only")
    plt.show()

  
  print len(original_orders), len(corrected_orders)
  if offset == None:
    offset = len(original_orders) - len(corrected_orders)
  offset = 0
  for i in range(offset, len(original_orders)):
    data = original_orders[i]
    data.cont = FittingUtilities.Continuum(data.x, data.y)
    try:
      model = corrected_orders[i-offset]
      header = corrected_headers[i-offset]
      print "Order = %i\nHumidity: %g\nO2 concentration: %g\n" %(i, header['h2oval'], header['o2val'])
    except IndexError:
      model = DataStructures.xypoint(x=data.x, y=numpy.ones(data.x.size))
      print "Warning!!! Telluric Model not found for order %i" %i

    if plot and i == plotorder:
      plt.figure(1)
      plt.plot(data.x, data.y/data.cont)
      plt.plot(model.x, model.y)

    if model.size() < data.size():
      left = numpy.searchsorted(data.x, model.x[0])
      right = numpy.searchsorted(data.x, model.x[-1])
      if right < data.size():
        right += 1
      data = data[left:right]
    elif model.size() > data.size():
      sys.exit("Error! Model size (%i) is larger than data size (%i)" %(model.size(), data.size()))

    badindices = numpy.where(numpy.logical_or(data.y <= 0, model.y < 0.05))[0]
    model.y[badindices] = data.y[badindices]/data.cont[badindices]
    
    data.y /= model.y
    original_orders[i] = data.copy()
  if plot:
    plt.show()
  return original_orders





def main1():
  if len(sys.argv) > 2:
    original = sys.argv[1]
    corrected = sys.argv[2]
  
    outfilename = "%s_telluric_corrected.fits" %(original.split(".fits")[0])
    print "Outputting to %s" %outfilename

    corrected_orders = Correct(original, corrected, offset=None)

    column_list = []
    if plot:
      plt.figure(2)
    for i, data in enumerate(corrected_orders):
      if plot and i == plotorder:
        plt.plot(data.x, data.y/data.cont)
      #Set up data structures for OutputFitsFile
      columns = {"wavelength": data.x,
                 "flux": data.y,
                 "continuum": data.cont,
                 "error": data.err}
      column_list.append(columns)
    HelperFunctions.OutputFitsFileExtensions(column_list, original, outfilename, mode="new")

    if plot:
      plt.title("Corrected data")
      plt.show()

  else:
    allfiles = os.listdir("./")
    corrected_files = [f for f in allfiles if "Corrected_" in f and f.endswith("-1.fits")]
    #original_files = [f for f in allfiles if any(f in cf for cf in corrected_files)]
    #hip_files = [f for f in allfiles if (f.startswith("HIP_") or f.startswith("HR_")) and not f.endswith("-0.fits")]

    #for original, corrected in zip(original_files, corrected_files):
    for corrected in corrected_files:
      idx = corrected.index("_")
      original = corrected[idx+1:].replace("-1.fits", "-0.fits")
      print original, corrected      
      
      outfilename = "%s_telluric_corrected.fits" %(original.split(".fits")[0])
      print "Outputting to %s" %outfilename
      
      corrected_orders = Correct(original, corrected, offset=None)
      
      column_list = []
      if plot:
        plt.figure(2)
      for i, data in enumerate(corrected_orders):
        if plot and i == plotorder:
          plt.plot(data.x, data.y/data.cont)
        #Set up data structures for OutputFitsFile
        columns = {"wavelength": data.x,
                   "flux": data.y,
                   "continuum": data.cont,
                   "error": data.err}
        column_list.append(columns)
      HelperFunctions.OutputFitsFileExtensions(column_list, original, outfilename, mode="new")

      if plot:
        plt.title(original)
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Flux")
        plt.show()





if __name__ == "__main__":
  main1()
