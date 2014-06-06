import matplotlib.pyplot as plt
import numpy
from scipy.optimize import leastsq
import SpectralTypeRelations
import HelperFunctions
import pyspeckit

plt.ion()

dy = 0.25

"""
  Function to get the equivalent width of the Na doublet lines (returns a list of size 2)
"""
def EW2(spectrum):
  fitfcn = lambda x, mu1, sig1, amp1, mu2, sig2, amp2: 1.0 - HelperFunctions.Gauss(x, mu1, sig1, amp1) - HelperFunctions.Gauss(x, mu2, sig2, amp2)
  errfcn = lambda pars, x, y, e: (y - fitfcn(x, *pars))**2 / e**2
  pars = [818.2,    #mu1
          0.05,     #sig1
          0.5,      #amp1
          819.4,    #mu2
          0.05,     #sig2
          0.5]      #amp2
  pars, success = leastsq(errfcn, pars, args=(spectrum.x, spectrum.y/spectrum.cont, spectrum.err/spectrum.cont))
  print pars
  fitted = fitfcn(spectrum.x, *pars)
  plt.plot(spectrum.x, spectrum.y/spectrum.cont, 'k-')
  plt.plot(spectrum.x, fitted, 'r--')
  plt.show()



def EW(spectrum, figure, title):
  x = pyspeckit.units.SpectroscopicAxis(spectrum.x, units='nm')
  spec = pyspeckit.Spectrum(data=spectrum.y/spectrum.cont, xarr=x, error=spectrum.err/spectrum.cont, units='nm')
  pl = spec.plotter(autorefresh=True, figure=figure)
  plt.title(title)

  # Fit the baseline in this region (continuum)
  spec.baseline(subtract=False, order=1, interactive=True)
  plt.show()
  done = raw_input("hit enter ")

  # Fit voigt profiles to the lines
  fitguesses=[-0.5, 818.3, 0.08, 0.0,
              -0.5, 819.4, 0.08, 0.0]
  spec.specfit(fittype='voigt', multifit=True, guesses=fitguesses)
  plt.draw()

  # Grab the fitted parameters and their errors
  fitpars = spec.specfit.modelpars
  fiterrs = spec.specfit.modelerrs

  # Sometimes, the lines 'switch' places, so check that
  if fitpars[1] > fitpars[5]:
    fitpars = fitpars[4:] + fitpars[:4]
    fiterrs = fiterrs[4:] + fiterrs[:4]

  # Save the parameters as a dictionary for each line
  line1 = {"Amplitude": [fitpars[0], fiterrs[0]],
           "Wavelength": [fitpars[1], fiterrs[1]],
           "Gaussian Width": [fitpars[2], fiterrs[2]],
           "Lorentzian Width": [fitpars[3], fiterrs[3]]}
  line2 = {"Amplitude": [fitpars[4], fiterrs[4]],
           "Wavelength": [fitpars[5], fiterrs[5]],
           "Gaussian Width": [fitpars[6], fiterrs[6]],
           "Lorentzian Width": [fitpars[7], fiterrs[7]]}

  #Determine the start and end of each line
  start1 = line1["Wavelength"][0] - 7*numpy.sqrt(line1["Gaussian Width"][0]**2 + line1["Lorentzian Width"][0]**2)
  end1 = line1["Wavelength"][0] + 7*numpy.sqrt(line1["Gaussian Width"][0]**2 + line1["Lorentzian Width"][0]**2)
  start2 = line2["Wavelength"][0] - 7*numpy.sqrt(line2["Gaussian Width"][0]**2 + line2["Lorentzian Width"][0]**2)
  end2 = line2["Wavelength"][0] + 7*numpy.sqrt(line2["Gaussian Width"][0]**2 + line2["Lorentzian Width"][0]**2)

  print start1, end1
  print start2, end2

  #Find the equivalent width for both lines
  ews = []
  for line in [[start1, end1], [start2, end2]]:
    xmin=numpy.searchsorted(spectrum.x, line[0])
    xmax=numpy.searchsorted(spectrum.x, line[1])
    ew = spec.specfit.EQW(xmin=xmin, xmax=xmax, fitted=False, continuum=1.0, plot=True)
    ews.append(ew)
  done = raw_input("hit enter ")
  print "\n\n"

  # Save output as a dictionary for each line
  line1["EQW"] = ews[0]
  line2["EQW"] = ews[1]
  return line1, line2
  

if __name__ == "__main__":
  #log = open("logfile.dat", "a")
  #log.write("#" + " "*20 + "\t\t\t\t   \t\t\t\t\tLine 1\t\t\t\t\t\t\t\t\t\t\t|\t\t\t\tLine 2\n")
  #log.write("#%s\t\t\t\tSpT\t\t\tAmplitude\t\t\tGaussWidth\t\t\tLorentzWidth\t\tEQW_1\t\t|\t\tAmplitude\t\t\tGaussWidth\t\t\tLorentzWidth\t\tEQW_2\n" %("filename".ljust(20)))
  #log.write("="*87 + "\n")

  # Set up some things
  MS = SpectralTypeRelations.MainSequence()
  fig = plt.figure()
  
  
  infile = open("stars.list")
  lines = infile.readlines()
  infile.close()

  spts = []
  data = []
  fnames = []
  for line in lines:
    if line.startswith("#"):
      continue

    fields = line.split()
    fname = fields[0]
    print "Reading %s with spt %s" %(fname, fields[1])
    spts.append(fields[1])
    rv = float(fields[2])
    orders = HelperFunctions.ReadFits(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
    for order in orders:
      if order.x[0] < 819 and order.x[-1] > 819:
        break
    order.x *= (1-rv/3e5)
    data.append(order)
    fnames.append(fname)

  #sorter = numpy.argsort(spts, key = lambda s: MS.SpT_To_Number(s))
  sorter = [i[0] for i in sorted(enumerate(spts), key=lambda s: MS.SpT_To_Number(s[1]))]
  sorted_spts = []
  sorted_data = []
  for i, idx in enumerate(sorter):
    spt = spts[idx]
    order = data[idx]
    sorted_spts.append(spt)
    sorted_data.append(order)

    #Plot the lines
    left = numpy.searchsorted(order.x, 817.5)
    right = numpy.searchsorted(order.x, 820)
    order = order[left:right]
    idx = numpy.argmin(order.y/order.cont)
    shift = order.x[idx] - 819.5
    order.x -= shift
    #fig.suptitle(fnames[idx])
    #fig.clf()

    
    #line1, line2 = EW(order, fig, fnames[idx])
    #log.write("%s\t\t%s\t\t%.3e +/- %.3e\t%.3e +/- %.3e\t\t%.3e +/- %.3e\t\t%.3e\t|\t%.3e +/- %.3e\t%.3e +/- %.3e\t\t%.3e +/- %.3e\t\t%.3e\n" \
    #          %(fnames[idx].ljust(20), 
    #            spt, 
    #            line1["Amplitude"][0], line1["Amplitude"][1],
    #            line1["Gaussian Width"][0], line1["Gaussian Width"][1],
    #            line1["Lorentzian Width"][0], line1["Lorentzian Width"][1],
    #            line1["EQW"],
    #            line2["Amplitude"][0], line2["Amplitude"][1],
    #            line2["Gaussian Width"][0], line2["Gaussian Width"][1],
    #            line2["Lorentzian Width"][0], line2["Lorentzian Width"][1],
    #            line2["EQW"],))
    
    print spt
    plt.plot(order.x, order.y/order.cont + i*dy, 'k-', lw=2)
    plt.text(817.23, 1.0+i*dy, spt, fontsize=14, color='black')
  plt.xlabel("Wavelength (nm)", fontsize=16)
  plt.ylabel("Normalized Flux + Constant", fontsize=16)
  plt.show()
  raw_input("Press enter ")
    
