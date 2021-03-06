import numpy
import sys
import os
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import TelluricFitter
import DataStructures
import Units
from astropy import units, constants
import HelperFunctions
import FittingUtilities

badregions = [[587, 590],
              [720, 720.3],
              [720.5, 721],
              [721.5, 721.6],
              [722.2, 722.4],
              [724.3, 724.5],
              [725, 725.4],
              [726, 726.22],
              [732.4, 732.8],
              [734.2, 737],
              [758,756.8],
              [806,810.7],
              [811.5,811.8],
              [817.56, 820.6]]
              
              
              
namedict = {"pressure": ["PRESFIT", "PRESVAL", "Pressure"],
                  "temperature": ["TEMPFIT", "TEMPVAL", "Temperature"],
                  "angle": ["ZD_FIT", "ZD_VAL", "Zenith Distance"],
                  "resolution": ["RESFIT", "RESVAL", "Detector Resolution"],
                  "h2o": ["H2OFIT", "H2OVAL", "H2O abundance"],
                  "co2": ["CO2FIT", "CO2VAL", "CO2 abundance"],
                  "o3": ["O3FIT", "O3VAL", "O3 abundance"],
                  "n2o": ["N2OFIT", "N2OVAL", "N2O abundance"],
                  "co": ["COFIT", "COVAL", "CO abundance"],
                  "ch4": ["CH4FIT", "CH4VAL", "CH4 abundance"],
                  "o2": ["O2FIT", "O2VAL", "O2 abundance"],
                  "no": ["NOFIT", "NOVAL", "NO abundance"],
                  "so2": ["SO2FIT", "SO2VAL", "SO2 abundance"],
                  "no2": ["NO2FIT", "NO2VAL", "NO2 abundance"],
                  "nh3": ["NH3FIT", "NH3VAL", "NH3 abundance"],
                  "hno3": ["HNO3FIT", "HNO3VAL", "HNO3 abundance"]}              
              
              

def FindAtmosphereFile(indexfile, date, time):
  index = numpy.recfromtxt(indexfile)
            
  year = int(date.split("-")[0])
  month = int(date.split("-")[1])
  day = int(date.split("-")[2])
  hour = int(time.split(":")[0])
  minute = int(time.split(":")[1])
  second = int(time.split(":")[2])
  obs_time = second + 60*minute + 3600*hour + 3600*24*day + 3600*24*30*month + 3600*24*365*year
  
  diff = 9e30
  filename = index[0]['f0']
  for i, entry in enumerate(index):
    hour = (entry['f4'] + entry['f5'])/2.0
    entry_time = 3600*hour + 3600*24*entry['f3'] + 3600*24*30*entry['f2'] + 3600*24*365*entry['f1']
    if abs(obs_time - entry_time) < diff:
      diff = obs_time - entry_time
      filename = entry['f0']
  return filename
  
  
              

if __name__ == "__main__":
  #Initialize fitter
  fitter = TelluricFitter.TelluricFitter(debug=False, debug_level=5)
  logfile = open("fitlog.txt", "w")
 
  fileList = []
  start = 0
  makenew = True
  exists = True
  edit_atmosphere = False
  for arg in sys.argv[1:]:
    if "-start" in arg:
      makenew = False
      start = int(arg.split("=")[-1])
    elif "-atmos" in arg:
      edit_atmosphere = True
      #atmosphere_fname = arg.split("=")[-1]
    else:
      fileList.append(arg)


  #START LOOPING OVER INPUT FILES
  for fname in fileList:
    logfile.write("Fitting file %s\n" %(fname))
    name = fname.split(".fits")[0]
    outfilename = "Corrected_%s.fits" %name
    if outfilename not in os.listdir("./"):
      exists = False

    orders = HelperFunctions.ReadFits(fname, errors="error", extensions=True, x="wavelength", y="flux")
    header = pyfits.getheader(fname)
    
    observatory = {"latitude": header['SITELAT'],
                   "altitude": header['SITEALT']/1000.0}
    fitter.SetObservatory(observatory)
    temperature = header['TEMPOUTS']+273.15
    pressure = 764.5
    humidity = 40.0
    resolution = 40000.0
    angle = numpy.arccos(1.0/header["AIRMASS"])*180.0/numpy.pi

  

    if edit_atmosphere:
      #Find the appropriate atmosphere profile from the date/time
      date = header["UT-DATE"]
      time = header["UT-START"] 
      atmosphere_fname = FindAtmosphereFile("Atmosphere_profile_directory.txt", date, time)
      
    
      #Read in NAM atmosphere profile information
      Pres,height,Temp,dew = numpy.loadtxt(atmosphere_fname, usecols=(0,1,2,3), unpack=True)
      sorter = numpy.argsort(height)
      height = height[sorter]
      Pres = Pres[sorter]
      Temp = Temp[sorter]
      dew = dew[sorter]
      
      #Convert dew point temperature to ppmv
      Pw = 6.116441 * 10**(7.591386*Temp/(Temp + 240.7263))
      h2o = Pw / (Pres-Pw) * 1e6
      humidity = 100.0* 10**(7.591386*(dew/(dew+240.7263) - Temp/(Temp+240.7263)))
      
      height /= 1000.0
      Temp += 273.15
      fitter.EditAtmosphereProfile("Temperature", height, Temp)
      fitter.EditAtmosphereProfile("Pressure", height, Pres)
      fitter.EditAtmosphereProfile("H2O", height, h2o)
      
      #Re-set the pressure and temperature guesses to the closest values
      idx = numpy.argmin(abs(observatory['altitude'] - height))
      pressure = Pres[idx]
      humidity = humidity[idx]
      
      

      
    #Adjust fitter values
    fitter.FitVariable({"h2o": humidity, 
                        "o2": 2.12e5,
                        "temperature": temperature})
    fitter.AdjustValue({"angle": angle,
                        "pressure": pressure,
                        "resolution": resolution})
    fitter.SetBounds({"h2o": [1.0, 99.9],
                      "o2": [5e4, 1e6],
                      "temperature": [temperature-50, temperature+50],
                      "resolution": [40000, 80000]})
    fitter.IgnoreRegions(badregions)
    models = []
    
    # Fit orders 22 and 23 to get the appropriate parameters.
    h2o = []
    o2 = []
    temperature = []
    resolution = []
    waveshifts = []
    wave0 = []
    chisquared = []
    start = 28
    for i, order in enumerate(orders[start:start+1]):
      print "\n***************************\nFitting order %i: " %(i+start)
      
      fitter.AdjustValue({"wavestart": order.x[0] - 10.0,
                          "waveend": order.x[-1] + 10.0})
      order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=2, highreject=10)
      fitter.ImportData(order)
      logfile.write("Fitting order %i with guassian line profiles\n" %(i+start)) 
      model, R = fitter.Fit(resolution_fit_mode="gauss", fit_primary=False, adjust_wave="model", continuum_fit_order=3, return_resolution=True)
     
      fitmodel = model.copy()
      
      # Save fitted parameters
      resolution.append(R)
      waveshifts.append(fitter.shift)
      wave0.append(fitter.data.x.mean())
      idx = numpy.where(numpy.array(fitter.parnames) == "h2o")[0]
      h2o.append(fitter.const_pars[idx])
      idx = numpy.where(numpy.array(fitter.parnames) == "o2")[0]
      o2.append(fitter.const_pars[idx])
      idx = numpy.where(numpy.array(fitter.parnames) == "temperature")[0]
      temperature.append(fitter.const_pars[idx])
      chisquared.append(1.0/fitter.chisq_vals[-1])
      
      
    # Take the average parameter values
    logfile.write("O2, H2O, T, R: \n")
    for x1, x2, x3, x4 in zip(o2, h2o, temperature, resolution):
      logfile.write("%g\t%g\t%g\t%g\n" %(x1, x2, x3, x4)  )
    print o2
    print h2o
    print temperature
    print resolution, '\n\n'
    chi2 = numpy.array(chisquared)
    o2 = numpy.array(o2)
    h2o = numpy.array(h2o)
    temperature = numpy.array(temperature)
    resolution = numpy.array(resolution)
    o2 = numpy.sum(o2*chi2)/numpy.sum(chi2)
    h2o = numpy.sum(h2o*chi2)/numpy.sum(chi2)
    temperature = numpy.sum(temperature*chi2)/numpy.sum(chi2)
    resolution = numpy.sum(resolution*chi2)/numpy.sum(chi2)
    
    # The waveshift difference is actually an error in the vacuum-->air correction
    # It looks like a velocity shift
    print waveshifts
    print wave0
    velshifts = 3e5*numpy.array(waveshifts)/numpy.array(wave0)
    print velshifts
    vel = numpy.sum(velshifts*chi2)/numpy.sum(chi2)
    #logfile.write("vshift = %g and %g" %(velshifts[0], velshifts[1]))
    logfile.write("vshift = %g\n" %vel)
    print "velocity shift: ", vel
    
    #Now, apply these values to all orders
    fitter.shift = 0.0
    for i, order in enumerate(orders):
      if i == start:
        model = fitmodel
      else:
        fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                            "waveend": order.x[-1] + 20.0,
                            "o2": o2,
                            "h2o": h2o,
                            "temperature": temperature,
                            "resolution": resolution})
        fitpars = [fitter.const_pars[j] for j in range(len(fitter.parnames)) if fitter.fitting[j] ]
      
        order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=2, highreject=10)
      
        fitter.ImportData(order)
        fitter.resolution_fit_mode = "gauss"
        model = fitter.GenerateModel(fitpars, nofit=True)
        model.x /= (1.0 + vel/3e5)
        model = FittingUtilities.RebinData(model, order.x)
      
      data = order.copy()
      primary = DataStructures.xypoint(x=order.x, y=numpy.ones(order.x.size))
      
      # Output
      #Set up data structures for OutputFitsFile
      columns = {"wavelength": data.x,
                 "flux": data.y,
                 "continuum": data.cont,
                 "error": data.err,
                 "model": model.y,
                 "primary": primary.y}
      
      header_info = []
      numpars = len(fitter.const_pars)
      for j in range(numpars):
        try:
          parname = fitter.parnames[j]
          parval = fitter.const_pars[j]
          fitting = fitter.fitting[j]
          header_info.append([namedict[parname][0], fitting, namedict[parname][2] ])
          header_info.append([namedict[parname][1], parval, namedict[parname][2] ])
        except KeyError:
          print "Not saving the following info: %s" %(fitter.parnames[j])
      
      
      if (i == 0 and makenew) or not exists:
        HelperFunctions.OutputFitsFileExtensions(columns, fname, outfilename, headers_info=[header_info,], mode="new")
        exists = True
      else:
        HelperFunctions.OutputFitsFileExtensions(columns, outfilename, outfilename, headers_info=[header_info,], mode="append")
     
      
      
    
    """
    
    Old method!
    
    #Make a test model, to determine whether/how to fit each value
    fitter.AdjustValue({"wavestart": orders[start].x[0]-20,
                        "waveend": orders[-1].x[-1]+20})
    fitpars = [fitter.const_pars[j] for j in range(len(fitter.parnames)) if fitter.fitting[j] ]
    test_model = fitter.GenerateModel(fitpars, nofit=True)
    

    #START LOOPING OVER ORDERS
    print "There are %i orders" %len(orders)
    for i, order in enumerate(orders[start:]):        
    
      print "\n***************************\nFitting order %i: " %(i+start)
      fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                          "waveend": order.x[-1] + 20.0})
      
      order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=2, highreject=10)
      primary = DataStructures.xypoint(x=order.x, y=numpy.ones(order.x.size))
      
      fitter.ImportData(order)

      #Determine how to fit the data from the initial model guess
      left = numpy.searchsorted(test_model.x, order.x[0])
      right = numpy.searchsorted(test_model.x, order.x[-1])
      
      model = DataStructures.xypoint(x=test_model.x[left:right], y=test_model.y[left:right])
      model_amplitude = 1.0 - min(model.y)
      print "Model amplitude: %g" %model_amplitude
      if i+start == 11 or i+start == 12:
        logfile.write("Fitting order %i with guassian line profiles\n" %(i+start)) 
        print "Fitting line profiles with gaussian profile"
        try:
          model = fitter.Fit(resolution_fit_mode="gauss", fit_primary=False, adjust_wave="model", continuum_fit_order=3)
        except ValueError:
          model = DataStructures.xypoint(x=order.x.copy(), y=numpy.ones(order.x.size))
        
        models.append(model)
        data = fitter.data
        
      else:
        logfile.write("Skipping order %i\n" %(i+start))
        print "Skipping order %i" %(i+start)
        data = order.copy()
        fitter.resolution_fit_mode = "gauss"
        model = fitter.GenerateModel(fitpars)
        #model = DataStructures.xypoint(x=order.x.copy(), y=numpy.ones(order.x.size))
        primary = model.copy()
      

      logfile.write("Array sizes: wave, flux, cont, error, model, primary\n")
      logfile.write("%i\n%i\n%i\n%i\n%i\n%i\n\n\n" %(data.x.size, data.y.size, data.cont.size, data.err.size, model.y.size, primary.y.size))
      #Set up data structures for OutputFitsFile
      columns = {"wavelength": data.x,
                 "flux": data.y,
                 "continuum": data.cont,
                 "error": data.err,
                 "model": model.y,
                 "primary": primary.y}
      namedict = {"pressure": ["PRESFIT", "PRESVAL", "Pressure"],
                  "temperature": ["TEMPFIT", "TEMPVAL", "Temperature"],
                  "angle": ["ZD_FIT", "ZD_VAL", "Zenith Distance"],
                  "resolution": ["RESFIT", "RESVAL", "Detector Resolution"],
                  "h2o": ["H2OFIT", "H2OVAL", "H2O abundance"],
                  "co2": ["CO2FIT", "CO2VAL", "CO2 abundance"],
                  "o3": ["O3FIT", "O3VAL", "O3 abundance"],
                  "n2o": ["N2OFIT", "N2OVAL", "N2O abundance"],
                  "co": ["COFIT", "COVAL", "CO abundance"],
                  "ch4": ["CH4FIT", "CH4VAL", "CH4 abundance"],
                  "o2": ["O2FIT", "O2VAL", "O2 abundance"],
                  "no": ["NOFIT", "NOVAL", "NO abundance"],
                  "so2": ["SO2FIT", "SO2VAL", "SO2 abundance"],
                  "no2": ["NO2FIT", "NO2VAL", "NO2 abundance"],
                  "nh3": ["NH3FIT", "NH3VAL", "NH3 abundance"],
                  "hno3": ["HNO3FIT", "HNO3VAL", "HNO3 abundance"]}
      header_info = []
      numpars = len(fitter.const_pars)
      for j in range(numpars):
        try:
          parname = fitter.parnames[j]
          parval = fitter.const_pars[j]
          fitting = fitter.fitting[j]
          header_info.append([namedict[parname][0], fitting, namedict[parname][2] ])
          header_info.append([namedict[parname][1], parval, namedict[parname][2] ])
        except KeyError:
          print "Not saving the following info: %s" %(fitter.parnames[j])
      
      
      if (i == 0 and makenew) or not exists:
        HelperFunctions.OutputFitsFileExtensions(columns, fname, outfilename, headers_info=[header_info,], mode="new")
        exists = True
      else:
        HelperFunctions.OutputFitsFileExtensions(columns, outfilename, outfilename, headers_info=[header_info,], mode="append")
    """
  logfile.close()
