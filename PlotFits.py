from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import sys
import HelperFunctions
import FittingUtilities
import numpy

def ReadFits(fitsfile, lower=812, upper=828):
  import DataStructures
  from astropy.io import fits
  hdulist = pyfits.open(fitsfile)
  print hdulist[0].header['object']
  data = hdulist[1].data
  x = data['wave'][0]/10.0
  y = data['flux_reduced'][0]
  err = data['err_reduced'][0]
  left = numpy.searchsorted(x, lower)
  right = numpy.searchsorted(x, upper)
  data = DataStructures.xypoint(x=x[left:right],
                                y=y[left:right],
                                err=err[left:right])
  data.cont = FittingUtilities.Continuum(data.x, data.y, fitorder=4, lowreject=1.5)
  #pylab.plot(data.x, data.y)
  #pylab.plot(data.x, data.cont)
  #pylab.show()
  return [data,]

if __name__ == "__main__":
  fileList = []
  tellurics = False
  normalize = False
  byorder = False   #Plots one order at a time
  pixelscale = False
  oneplot = False
  filetype = "Mine"
  for arg in sys.argv[1:]:
    if "tellcorr" in arg:
      tellurics = True
    elif "-norm" in arg:
      normalize = True
    elif "-order" in arg:
      byorder = True
    elif "-pix" in arg:
      pixelscale = True
      #byorder = True
    elif "-one" in arg:
      oneplot = True
    elif "-type" in arg:
      filetype = "Native"
    else:
      fileList.append(arg)

  linestyles = ['k-', 'r-', 'b-', 'g-']

  for fnum, fname in enumerate(fileList):
    ls = linestyles[fnum%len(linestyles)]
    if filetype == "Mine":
      orders = HelperFunctions.ReadFits(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
    else:
      orders = ReadFits(fname)
    print fname, len(orders)
    if not oneplot:
      plt.figure(fnum)
      plt.title(fname)
    if tellurics:
      model = HelperFunctions.ReadFits(fname, extensions=True, x="wavelength", y="model")
    for i, order in enumerate(orders):
      #order.cont = FittingUtilities.Continuum(order.x, order.y, lowreject=3, highreject=3)
      if pixelscale:
        order.x = numpy.arange(order.size())
      if tellurics:
        plt.plot(order.x, order.y/order.cont, 'k-')
        plt.plot(order.x, model[i].y, 'r-')
      else:
        if normalize:
          plt.plot(order.x, order.y/order.cont)
          plt.text(order.x.mean(), 1.1, str(i+1))
        else:
          plt.plot(order.x, order.y, ls)
          plt.plot(order.x, order.cont)
      if byorder:
        plt.title("Order %i" %i)
        plt.show()
  if not byorder:    
    plt.show()
