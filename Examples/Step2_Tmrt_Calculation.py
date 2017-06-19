# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 16:23:35 2017
@author: Tiffany Sin 2017

TMRT for idealized Matrices. This example is dependent on Step1_Model_SetUp. 

"""
#import pyliburo
#import ExtraFunctions
import datetime
import time
from datetime import datetime
from datetime import datetime, date, time

current_path = os.path.dirname("__file__")
parent_path = os.path.abspath(os.path.join(current_path, os.pardir))

install_and_import('easygui')

import easygui
print("Navigate to thermal comfort file ")
thermalcomfortpath=easygui.fileopenbox()
import imp #importing is causing heaps of problems so we use imp to help us keep our directories straight. 
#thermalcomfort = imp.load_source('thermalcomfort',parent_path+'\\thermalcomfort.py')
#reload(thermalcomfort) #still causing us problems so we reload thermalcomfort just in case

import time
#%%

#%% Ped Locations 
simdate = 'Today'
def install_and_import(package):
    import importlib
    try:
        importlib.import_module(package)
        print "Checking for package.."
    except ImportError:
        import pip
        pip.main(['install', package])
        print "Package not available. Importing package.."
    finally:
        globals()[package] = importlib.import_module(package)
        print "Package installed"


#modelinputs_pathfile = raw_input('Enter the full path of the model input file ')
#model_inputs=pd.read_csv(modelinputs_pathfile)  

latitude = model_inputs.latitude[0]
longitude = model_inputs.longitude[0]
casetime = model_inputs.time[0]

install_and_import('tzlocal')
from tzlocal import get_localzone
tz = get_localzone() # Prints Asia/Singapore
d = datetime.now(tz) # Prints date and time
utc_offset = d.utcoffset().total_seconds()
timezone=-(utc_offset/3600.0)
print timezone

# When calculating spatial variation of thermal comfort (i.e. not temporal variation), solar parameters only need to be calculated once
time1 = time.clock()
solarparam = thermalcomfort.calc_solarparam(casetime,latitude,longitude) 
time2 = time.clock()
print 'solar_parameters() CALCULATION TIME: ',(time2-time1)/60.0, 'minutes'

cases = [myexperiment]
for config in cases:
    #Initializing... 
    pedkeys  = config['pedkeys'] #pedkeys are the coordinates at which Tmrt will be calculated.
    compound = config['model']
    pdTs = config["Tsurf"]
    pdReflect = config["Refl"]
    #The air temperature provided is a bulk value. Repeat the value across all coordinates to create a pdcoord for Tair. 
    pdTa = thermalcomfort.pdcoords_from_pedkeys(pedkeys, np.array([config["Tair"]]*len(pedkeys))) 
    
    #initialize a pdcoord for Tmrt that is filled with zeros
    config['TMRT'] =thermalcomfort.pdcoords_from_pedkeys(pedkeys) 
    for index, row in config['TMRT'].data.iterrows(): #For each pedestrian coordinate...
       
        time1 = time.clock()
        pedkey = (row.x,row.y,row.z) #retrieve the pedestrian's coordinate
        results = thermalcomfort.all_mrt(pedkey,compound,pdTa,pdReflect,pdTs,solarparam,model_inputs,ped_properties,gridsize=3) #this calculates all steps necessary for MRT calculation.
        #see code in part 2 of thermalcomfort.py to see step-by-step explanation of the calculation.        
        row.v = results.TMRT[0]
        time2 = time.clock()
        tottime = (time2-time1)/60.0
        print  index, results.TMRT[0], tottime
    
    #Save results to a csv file    
    config['TMRT'].data.to_csv(config['name']+ '_'+simdate +'_TMRT.csv')



