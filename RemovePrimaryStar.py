#!/opt/local/bin/python
import numpy
import os
import sys
import pylab
from scipy.interpolate import UnivariateSpline
import SpectralTypeRelations
from collections import defaultdict
from PlotBlackbodies import Planck
import Units
import DataStructures
import Correlate
import MakeModel
import RotBroad_Fast as RotBroad
import time
import FittingUtilities
import HelperFunctions

"""
   This function removes a primary model from data  
"""



def GetModel(data, model, vel=0.0, vsini=15*Units.cm/Units.km, resolution=20000):
  #Broaden
  model = RotBroad.Broaden(model, vsini=vsini, findcont=True)
  
  #Interpolate Model, and note the endpoints
  first = model.x[0]*(1.+vel/Units.c)
  last = model.x[-1]*(1.+vel/Units.c)
  a = []
  for i in range(1, model.x.size):
    a.append(model.x[i] - model.x[i-1])
  model_fcn = UnivariateSpline(model.x, model.y, s=0)

  data2 = []
  for order in data:
    data2.append(order.copy())

  ##############################################################
  #Begin main loop over the orders
  ##############################################################
  output_orders = []
  for i in range(len(data2)):
    order = data2[i]
    order.cont = FittingUtilities.Continuum(order.x, order.y, lowreject=2, highreject=5)

    model2 = DataStructures.xypoint(x=order.x, y=model_fcn(order.x*(1.+vel/Units.c)))
    
    #Get model continuum in this section
    model2.cont = FittingUtilities.Continuum(model2.x, model2.y, lowreject=1, fitorder=3)
    
    #Reduce resolution
    model2 = FittingUtilities.ReduceResolution(model2.copy(), resolution)

    #Fit velocity with a straight shift via cross-correlation
    ycorr = numpy.correlate(order.y/order.cont-1.0, model2.y/model2.cont-1.0, mode="full")
    xcorr = numpy.arange(ycorr.size)
    lags = xcorr - (order.x.size-1)
    distancePerLag = (order.x[-1] - order.x[0])/float(order.x.size)
    offsets = -lags*distancePerLag
    offsets = offsets[::-1]
    ycorr = ycorr[::-1]
    fit = numpy.poly1d(numpy.polyfit(offsets, ycorr, 3))
    ycorr = ycorr - fit(offsets)
    maxindex = ycorr.argmax()
    if i == 28:
      print "Offset = %g nm" %offsets[maxindex]
    if abs(offsets[maxindex]) < 0.5e-9:
      model2 = DataStructures.xypoint(x=order.x, y=model_fcn(order.x*(1.+vel/Units.c)+offsets[maxindex]))
    model2.cont = FittingUtilities.Continuum(model2.x, model2.y, lowreject=1, fitorder=3)
    model2 = FittingUtilities.ReduceResolution(model2.copy(), resolution)

    #Scale using Beer's Law
    #line_indices = numpy.where(model2.y / model2.cont < 0.96)[0]
    #if len(line_indices) > 0:
    #  scale = numpy.median(numpy.log(order.y[line_indices]/order.cont[line_indices]) / numpy.log(model2.y[line_indices] / model2.cont[line_indices]) )
    #  model2.y = model2.y**scale
    #  model2.cont = model2.cont**scale
    
    #print "Order %i: Scale = %g" %(i, scale)

    output_orders.append(model2.copy())

  #pylab.show()
  return output_orders




#if __name__ == "__main__":
def main1(datafile, p_spt="A0", vsini=15.0, vel=0.0, metal=0.0): 
  
  modeldir = os.environ["HOME"] + "/School/Research/Models/Sorted/Stellar/Vband/"
  files = os.listdir(modeldir)
  modelfiles = defaultdict(list)
  Modeler = MakeModel.Modeler()
  for fname in files:
    try:
      temperature = float(fname.split("lte")[-1].split("-")[0])*100
    except ValueError:
      try:
        temperature = float(fname.split("lte")[-1].split("+")[0])*100
      except ValueError:
        print "Skipping file %s" %fname
        continue
    modelfiles[temperature].append(modeldir+fname)

  #Read in data
  orders_original = list(HelperFunctions.ReadFits(datafile, extensions=True, x="wavelength", y="flux", errors="error"))

  #Do a rough telluric correction in the data
  tellurics = []
  for i, order in enumerate(orders_original):
    lowfreq = 1e7/(order.x[-1] + 10)
    highfreq = 1e7/(order.x[0] - 10)
    model = Modeler.MakeModel(humidity=20.0, lowfreq=lowfreq, highfreq=highfreq)
    model = FittingUtilities.ReduceResolution(model, 40000)
    model = FittingUtilities.RebinData(model, order.x)
    order.y /= model.y
    tellurics.append(model)
    orders_original[i] = order.copy()

  orders_original = tuple(orders_original[::-1])
  
  #Get the best logg and temperature for a main sequence star with the given spectral type
  MS = SpectralTypeRelations.MainSequence()
  p_temp = MS.Interpolate(MS.Temperature, p_spt)
  p_mass = MS.Interpolate(MS.Mass, p_spt)
  radius = MS.Interpolate(MS.Radius, p_spt)
  logg = numpy.log10(Units.G*p_mass*Units.Msun/(radius*Units.Rsun)**2)

  #Find the best fit temperature
  Ts = sorted(modelfiles.keys())
  bestidx = numpy.argmin((numpy.array(Ts) - p_temp)**2)
  filenames = []
  for i in (bestidx-1, bestidx, bestidx+1):
    if i < 0 or i > len(Ts):
      continue
    for fname in modelfiles[Ts[i]]:
      filenames.append(fname)
  #temperature = modelfiles.keys()[numpy.argmin((modelfiles.keys() - p_temp)**2)]
  #print p_temp, temperature

  #filenames = modelfiles[temperature]
  
  """
  #Find the best logg
  best_logg = 9e9
  indices = []
  for fname in modelfiles[temperature]:
    logg_val = float(fname.split("lte")[-1].split("-")[1][:3])
    if numpy.abs(logg_val - logg) < numpy.abs(best_logg - logg):
      best_logg = logg_val
  for i, fname in enumerate(modelfiles[temperature]):
    logg_val = float(fname.split("lte")[-1].split("-")[1][:3])
    if logg_val == best_logg:
      indices.append(i)
  filenames = []
  for i in indices:
    filenames.append(modelfiles[temperature][i])
  """

  chisquareds = []
  models = []
  for modelfile in filenames:
    #Read in model
    print "Model file: %s" %modelfile
    x,y = numpy.loadtxt(modelfile, usecols=(0,1), unpack=True)
    x = x*Units.nm/Units.angstrom
    y = 10**y
    model = DataStructures.xypoint(x=x, y=y)

    orders = GetModel(list(orders_original), model, vel=vel*Units.cm/Units.km, vsini=vsini*Units.cm/Units.km)
    chisq = 0.0
    i = 0
    for model, data in zip(orders, orders_original):
      if i != 15 and i != 16 and i != 44:
        data.cont = FittingUtilities.Continuum(data.x, data.y, lowreject=2, highreject=2)
        chisq += numpy.sum((data.y - model.y*data.cont/model.cont)**2 / data.err**2) / float(data.size())
        i += 1
    chisquareds.append(chisq)
    models.append(orders)

  best_index = numpy.argmin(chisquareds)
  orders = models[best_index]
  for x2, fname in zip(chisquareds, filenames):
    print fname.split("/")[-1], "\t", x2

  
  if "-" in datafile:
    i = 0
    itemp = 0
    while i >= 0:
      i = datafile.find("-", i+1, len(datafile))
      if i > 0:
        itemp = i
    idx = int(datafile[itemp:].split(".fits")[0])
    outfilename = "%s-%i.fits" %(datafile[:itemp], idx+1)
  else:
    outfilename = "%s-0.fits" %(datafile.split(".fits")[0])
  print "Outputting to %s" %outfilename

  orders_original = orders_original[::-1]
  orders = orders[::-1]

  column_list = []
  Modeler = MakeModel.Modeler()
  for i, original in enumerate(orders_original):
    """
    #Make a telluric model for this order
    lowfreq, highfreq = 1e7/original.x[-1], 1e7/original.x[0]
    telluric = Modeler.MakeModel(humidity=50.0, lowfreq=lowfreq, highfreq=highfreq)
    telluric = FittingUtilities.ReduceResolution(telluric, 20000)
    telluric = FittingUtilities.RebinData(telluric, original.x)
    """

    original.y *= tellurics[i].y
    
    #original = orders_original[i+2]
    original.cont = FittingUtilities.Continuum(original.x, original.y, lowreject=2, highreject=2)
    model = orders[i]
    if 819 > original.x[0] and 819 < original.x[-1]:
      pylab.figure(1)
      pylab.plot(original.x, original.y/original.cont, 'k-')
      pylab.plot(model.x, model.y/model.cont, 'r-')
      #pylab.plot(telluric.x, telluric.y, 'g-')
      pylab.figure(2)
      pylab.plot(original.x, original.y/(original.cont * model.y/model.cont), 'k-')
      #pylab.plot(telluric.x, telluric.y, 'g-')
    original.y /= model.y/model.cont
    
    columns = {"wavelength": original.x,
               "flux": original.y,
               "continuum": original.cont,
               "error": original.err}
    column_list.append(columns)
  
  
  # Output
  HelperFunctions.OutputFitsFileExtensions(column_list, datafile, outfilename, mode="new")
    
  #pylab.figure(1)
  #pylab.savefig("Comparison.pdf")
  #pylab.figure(2)
  #pylab.savefig("Corrected.pdf")
  pylab.show()



if __name__ == "__main__":
  import os
  import sys
  if len(sys.argv) > 1:
    p_spt="A0"
    vsini = 15
    vel = 0.0
    metal = 0.0
    fileList = []
    for arg in sys.argv[1:]:
      if "primary" in arg:
	p_spt = arg.split("=")[-1]
      elif  "vsini" in arg:
	vsini = float(arg.split("=")[-1])
      elif "rv" in arg:
	vel = float(arg.split("=")[-1])
      elif "metal" in arg:
	metal = float(arg.split("=")[-1])
      else:
	fileList.append(arg)
    
    for fname in fileList:
      main1(fname, p_spt=p_spt, vsini=vsini, vel=vel, metal=metal)

  else:
    infile = open("par.list")
    lines = infile.readlines()
    infile.close()
    for line in lines:
      if line.startswith("#"):
	continue

      fields = line.split()
      fname = fields[0]
      p_spt = fields[1]
      vel = float(fields[2])
      vsini = float(fields[3])
      main1(fname, p_spt=p_spt, vsini=vsini, vel=vel, metal=0.0)

