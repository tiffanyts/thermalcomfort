# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 09:43:35 2017
@author: Tiffany Sin 2017

After calculating TMRT and importing initial information, SET can be calculated at every pedestrian location. This example is dependent on previous steps. 
"""
import time
import numpy

thermalcomfortpath=raw_input("Enter the path of thermal comfort file")
extrafunctionspath=raw_input("Enter the path of extra functions file")
import imp #importing is causing heaps of problems so we use imp to help us keep our directories straight. 
thermalcomfort = imp.load_source('thermalcomfort',thermalcomfortpath)
ExtraFunctions = imp.load_source('ExtraFunctions',extrafunctionspath)
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
#%% 
#import Importing
simdate = 'Today'


#%% Microclimate
cases = [myexperiment]
for config in cases:
    Tmrt = config['TMRT']
    pedkeys  = config['pedkeys'] #imported previously. 
    Ta = config['Tair'] #constant
    V = config['wind'] 
    RH = 50
    config['SET']  = thermalcomfort.pdcoords_from_pedkeys(pedkeys) #initialize a pdcoord for SET that is filled with zeros

    for index,row in config['SET'].data.iterrows(): #calculate SET line by line along the pdcoord (i.e. for each coordinate on the grid)
        time1 = time.clock()
        pedkey = (row.x,row.y,row.z)
        microclimate = pd.DataFrame({
        'T_air':[Ta],
        'wind_speed':[np.mean(abs(V.val_at_coord(pedkey, radius = 0.2).v))],
       'mean_radiant_temperature': [np.mean(Tmrt.val_at_coord(pedkey, radius = 1).v)-273.15],
        'RH':[50],  
        })
    
        row.v = thermalcomfort.calc_SET(microclimate,ped_properties)    
        time2 = time.clock()
        tottime = (time2-time1)/60.0
        print  ' row ' , index,  ' | SET is ', row.v ,'. TIME TAKEN', tottime

    config['SET'].data.to_csv(config['name'] + '_'+simdate +'_SET.csv')



