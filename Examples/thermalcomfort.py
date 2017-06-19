# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 09:28:56 2016

@author: Tiffany Sin 2017

Functions for calculating Mean Radiant Temperature and Standard Effective Temperature. 
Please see examples for the implementation of OTC3D. 

This code is in three parts. After importing the necessary modules, the first section determines the pdcoord helper class, which is used to pass spatial data around.
The second section provides the mean radiant temperature calculation. It is comprised of several smaller functions that are combined in sequence in the function all_mrt()
The final section consists of the SET calculation, which relies on inputs as defined in the example.  

"""
#
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
from matplotlib.path import Path
import matplotlib.patches as patches

import numpy as np
import pandas as pd
#import pyliburo
#import pvlib

#import pvlib

import datetime
import time

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

install_and_import('pyliburo')
install_and_import('pvlib')
install_and_import('scipy')
install_and_import('pandas')
install_and_import('numpy')

from OCC.Display import OCCViewer
#from ExtraFunctions import *


#%% Part 1)  Handling input and output data
# To deal with spatial variability, we introduce the helper class 'pdcoord', which is a Pandas DataFrame of four columns (x, y, z, and value) columns.
# This allows us to import and pass around detailed spatial data, like surface temperatures, air temperature, and wind speeds, and also to calculate SET and Tmrt in detail.  
# We also include several functions to help visualize and manipulate our data.
# A crash course on Pandas DataFrames can be found here: https://www.datacamp.com/community/tutorials/pandas-tutorial-dataframe-python#gs._wHPqvc

def read_pdcoord(csv_inputfile,separator = ','):
    """ Function to import data into dataframe structure. 
    To create an empty dataframe, input the string 'Empty' instead. """
    if type(csv_inputfile)== pd.core.frame.DataFrame: #If inputdata is already a dataframe, rename columns to ('x','y','z','v')
        csv_input = csv_inputfile 
        csv_input.columns =['x','y','z','v']
    elif csv_inputfile == 'Empty':  #Create an empty dataframe. 
        csv_input = pd.DataFrame(columns = ['x','y','z','v'])
    else:
        try:     # if it's a csv file, read csv directly into a dataframe
            csv_input= pd.read_csv(csv_inputfile,sep=separator,skiprows=[0],header=None,names = ['x','y','z','v']) # r"\s+" is whitespace
        except IOError:   # If it's a numpy array, convert to pdframe
            try:
                csv_input = pd.DataFrame(np.array(csv_inputfile), columns = ['x','y','z','v'])
            except IOError:
                print "Incorrect Input: read_pdcoord requires the input of a 4-column csv file, numpy array, or pandas dataframe." 
                return None
    csv_input.sortlevel(axis=0,inplace=True,sort_remaining=True)
    return csv_input

class pdcoord(object):
    """ Helper class for all x,y,z,v input files used in thermal comfort analysis
    To initialize a new pdcoord with input data, use pdcoord(csv_input)
    To create a new pdcoord based on only coordinates, use pdcoords_from_pedkeys() (see below)
    """
    
    def __init__(self, csv_input, sep =','):
        self.data = read_pdcoord(csv_input,separator = sep)
    
    def recenter(self,origin=(0,0)):
        """ Shifts coordinate data such that the bottom left corner is at the origin """
        print 'shifted by x:' ,min(self.data.x)+ origin[0], ' and y:',min(self.data.y) + origin[1]
        self.data['x'] = self.data['x'] - min(self.data.x) + origin[0]
        self.data['y'] = self.data['y'] - min(self.data.y) + origin[1]
        return self
    
    def repeat_outset(self,unit=1):
        """ extends in 8 directions for repeating data. Used to repeat surface data to account for more buildings in the area of study"""
        shiftLR = max(self.data.x) - min(self.data.x) + unit
        shiftUD = max(self.data.y) - min(self.data.y) + unit
        
        #clockwise
        Q1 = pd.DataFrame(np.array([self.data.x-shiftLR, 
                           self.data.y+shiftUD,
                           self.data.z, 
                           self.data.v]).T, 
                           columns = ['x','y','z','v'])
        Q2 = pd.DataFrame(np.array([self.data.x, 
                           self.data.y+shiftUD,
                           self.data.z, 
                           self.data.v]).T, 
                           columns = ['x','y','z','v'])
        Q3 = pd.DataFrame(np.array([self.data.x+shiftLR, 
                           self.data.y+shiftUD,
                           self.data.z, 
                           self.data.v]).T, 
                           columns = ['x','y','z','v'])
        Q4 = pd.DataFrame(np.array([self.data.x+shiftLR, 
                           self.data.y,
                           self.data.z, 
                           self.data.v]).T, 
                           columns = ['x','y','z','v'])
        Q5 = pd.DataFrame(np.array([self.data.x+shiftLR, 
                           self.data.y-shiftUD,
                           self.data.z, 
                           self.data.v]).T, 
                           columns = ['x','y','z','v'])
        Q6 = pd.DataFrame(np.array([self.data.x, 
                           self.data.y-shiftUD,
                           self.data.z, 
                           self.data.v]).T, 
                           columns = ['x','y','z','v'])
        Q7 = pd.DataFrame(np.array([self.data.x-shiftLR, 
                           self.data.y-shiftUD,
                           self.data.z, 
                           self.data.v]).T, 
                           columns = ['x','y','z','v'])
        Q8 = pd.DataFrame(np.array([self.data.x-shiftLR, 
                           self.data.y,
                           self.data.z, 
                           self.data.v]).T, 
                           columns = ['x','y','z','v'])

        for Q in [Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8]:
            self.data =self.data.append(Q,ignore_index=True)
        return self

    def val_at_coord(self,listcoord, radius = 1):
        """ Returns the values of a pdcoord at a target coordinate, within a distance of the radius. 
        Input coordinates as a list ([X,Y,Z], [X], or [X,Y]). If only X or only X,Y are given, returns selected dataframe. 
        Range of selection can be widened or narrowed with the radius.
        For a single value result, take the mean of the resulting pdcoord. e.g. surfpdcoord.val_at_coord(listcoord,radius = 5).v.mean()"""
        minx, miny, minz = map(lambda a: a-radius, listcoord) #make bins about the target listcoord.
        maxx, maxy, maxz = map(lambda a: a+radius, listcoord) 
    
        if len(listcoord) == 1: return self.data[(self.data.x <=maxx) & (self.data.x >= minx)] #select values that lie within the location bins.
        elif len(listcoord) == 2: return self.data[(self.data.x <=maxx) & (self.data.x >= minx) & (self.data.y <=maxy) & (self.data.y >= miny) ]
        elif len(listcoord) == 3: 
            try: return self.data[(self.data.x <=maxx) & (self.data.x >= minx) & (self.data.y <=maxy) & (self.data.y >= miny) &(self.data.z <=maxz) & (self.data.z >= minz)]
            except ValueError: print "No data found at that coordinate" 

    def scatter3d(self,title='',size=40,model=[]):
        """Returns a scatterplot of the pdcoord. Useful for visualizing 3D surface data """
        font = {'weight' : 'medium',
                'size'   : 22}
        #plt.rc('font', **font)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d'); ax.pbaspect = [1, 1, 1] #always need pbaspect
        ax.set_title(title)
        p = ax.scatter(list(self.data.x), list(self.data.y),list(self.data.z),c = list(self.data.v),edgecolors='none', s=size, marker = ",", cmap ='jet')
        ax.view_init(elev=90, azim=-89)
        ax.set_xlabel('X axis'); ax.set_ylabel('Y axis'); ax.set_zlabel('Z axis')
        fig.colorbar(p)
        
        plt.draw()
        try:
            vertices = [(vertex.X(), vertex.Y(),vertex.Z()) for vertex in pyliburo.py3dmodel.fetch.vertex_list_2_point_list(pyliburo.py3dmodel.fetch.topos_frm_compound(model)["vertex"])]
            V1,V2,V3 = zip(*vertices)            
            p = ax.plot_wireframe(V1,V2,V3 )

        except TypeError:
            pass        
        
        return fig
        
    def contour(self,title='',cbartitle = '',model=[], zmax = None, zmin = None, filename = None, resolution = 1, unit_str = '', bar = True):
        """ Returns a figure with contourplot of 2D spatial data. Insert filename to save the figure as an image. Increase resolution to increase detail of interpolated data (<1 to decrease)"""

        font = {'weight' : 'medium',
                'size'   : 22}
        
        xi = np.linspace(min(self.data.x), max(self.data.x),len(set(self.data.x))*resolution)
        yi = np.linspace(min(self.data.y), max(self.data.y),len(set(self.data.y))*resolution)
        
    
        zi = ml.griddata(self.data.x, self.data.y, self.data.v.interpolate(), xi, yi,interp='linear')
        
        fig = plt.figure()
        plt.rc('font', **font)
        plt.title(title)
        plt.contour(xi, yi, zi, 15, linewidths = 0, cmap=plt.cm.bone)
        plt.pcolormesh(xi, yi, zi, cmap = plt.get_cmap('rainbow'),vmax = zmax, vmin = zmin)
        if bar: cbar = plt.colorbar(); cbar.ax.set_ylabel(cbartitle)

        plt.absolute_import
        try:
            vertices = [(vertex.X(), vertex.Y()) for vertex in pyliburo.py3dmodel.fetch.vertex_list_2_point_list(pyliburo.py3dmodel.fetch.topos_frm_compound(model)["vertex"])]
            shape = patches.PathPatch(Path(vertices), facecolor='white', lw=0)
            plt.gca().add_patch(shape)
        except TypeError:
            pass
        plt.show()
        
        try:
            fig.savefig(filename)
        except TypeError:
            return fig
    
#    def plot_along_line(self,X,Y, tick_list):
#        V = self.data.v
#        plt.plot(heights, SVFs_can, label='Canyon')
    
def pdcoords_from_pedkeys(pedkeys_np, values = np.zeros(0)):
    """ Creates a pdcoord using coordinates from a numpy array (np.array([(X1,Y1,Z1),(X2,Y2,Z2),...])) and values from a 1D numpy array of the same size. Can be thought of as a 'zip' function to attach a 1D series of data to coordinates  """
    if not values.size:
        fill = [np.nan]*len(pedkeys_np)
    else: fill = values
    empty = pdcoord(zip(pedkeys_np.transpose()[0], pedkeys_np.transpose()[1], pedkeys_np.transpose()[2],fill))
    return empty

#%% Part 2) Radiation Model Functions

#1) Calculate solar parameters

#For each pedestrian location:
#2) Calculate shadow
#3) Use fourpiradiation() to calculate sky view factor, list of visible intecept,
#    and number of remaining directions that hit the ground (groundN)
#4) Use call_values to call values of Ts, reflect, wall_albedo, wall_emissivity 
#    at intecepted locations on building surfaces
#5) Use calc_radiation_from_values() to return longwave & shortwave radiation from surfaces
#6) Use wall_emissivity*sigma*Ts**4*groundN/Ndir to calculate remaining ground longwave radiation
#7) Use E_dif*solarvf + E_sol*solarvf*shadowint for direct and diffuse solar radiation
#8) Tmrt = ((Eshort*(1-ped_albedo)+Elong)/sigma)**(1/4.)

Ndir = 1000
unitball = pyliburo.skyviewfactor.tgDirs(Ndir)
sigma =5.67*10**(-8)  

#%%

def calc_solarparam(time_str,latitude,longitude, UTC_diff=0, groundalbedo=0.18):
    """ This function uses PVLib to calculate solar parameters. 
    Returns a DataFrame of solar vector, solar view factor, 
    direct solar radiation intensity, and diffuse solar radiation 
    intensities from the sky and the ground """

    time_shift = datetime.timedelta(hours=UTC_diff) #SGT is UTC+8    
    thistime =pd.DatetimeIndex(start=time_str, end=time_str, freq='1min') # pd.to_datetime(pd.Timestamp(casetime)) # pd.DatetimeIndex([pd.Timestamp(np.datetime64(datetime.datetime(y,mo,d,h,mi) + time_shift), tz='UTC')])  
    #thistime = pd.DatetimeIndex([pd.to_datetime(time_str,format= '%B %d %Y %H:%M')])
    thisloc = pvlib.location.Location(51.4826,  0.0077, tz='UTC', altitude=0, name=None)
    solpos = thisloc.get_solarposition(thistime)    
    
    sunpz = np.sin(np.radians(solpos.elevation[0])); hyp = np.cos(np.radians(solpos.elevation[0]))
    sunpy = hyp*np.cos(np.radians(solpos.azimuth[0]))
    sunpx = hyp*np.sin(np.radians(solpos.azimuth[0]))
    
    solar_pmt =  thisloc.get_clearsky(thistime,model='ineichen') #
    E_sol= solar_pmt.dni[0] #direct normal solar irradiation  [W/m^2]
    Ground_Diffuse = pvlib.irradiance.grounddiffuse(90,solar_pmt.ghi,albedo =groundalbedo)[0]  #Ground Reflected Solar Irradiation - vertical asphalt surface [W/m^2]
    Sky_Diffuse =  pvlib.irradiance.isotropic(90, solar_pmt.dhi)[0] #Diffuse Solar Irradiation - vertical surface[W/m^2].        
    
    #Formula 9 in Huang et. al. for a standing person, largely independent of gender, body shape and size. For a sitting person, approximately 0.25
    solarvf=abs(0.0355*np.sin(solpos.elevation[0])+2.33*np.cos(solpos.elevation[0])*(0.0213*np.cos(solpos.azimuth[0])**2+0.00919*np.sin(solpos.azimuth[0])**2)**(0.5)); 
    results = pd.DataFrame({
    'solarvector':[(sunpx,sunpy,sunpz)],
    'solarviewfactor':[solarvf],
    'direct_sol':[E_sol],
    'diffuse_frm_sky':[Sky_Diffuse],
    'diffuse_frm_ground':[Ground_Diffuse]
    })    
    return results
#%%

def check_shadow(key, model, solarvector):
    occ_interpt, occ_interface = pyliburo.py3dmodel.calculate.intersect_shape_with_ptdir(model,key,solarvector)
    if occ_interpt != None: return 0
    else: return 1 

def get_shadow(pedestrian_keys, model,solar_vector):
    """ Returns a dataframe of shadowed (0) and sunlit (1) locations. Ignores points that are on the wall (treats them as not shadowed)  """
    shadow = pdcoord(zip(pedestrian_keys.transpose()[0], pedestrian_keys.transpose()[1], pedestrian_keys.transpose()[2],np.zeros(len(pedestrian_keys)) ))
    shadow.data['v'] = shadow.data.apply(lambda row: check_shadow((row['x'], row['y'], row['z']),model, solar_vector), axis=1)
    return shadow

def fourpiradiation(key, model):
    """ returns SVF, number of ground points (N), and list of intercepts. 
    For uniform ground temperature, do not include ground surface in model. Longwave irradiance from ground can be calculated as emissivity*sigma*groundtemp**4*N/Ndir 
    If ground temperature is not uniform, include the ground in the model, and radiation will be calculated with the other surfaces. """ 
    sky=0.; ground = 0.; intercepts=[]
    for direction in unitball.getDirUpperHemisphere():
        (X,Y,Z) = (direction.x,direction.y,direction.z)
        occ_interpt, occ_interface = pyliburo.py3dmodel.calculate.intersect_shape_with_ptdir(model,key,(X,Y,Z))
        if occ_interpt != None: intercepts.append([occ_interpt.X(), occ_interpt.Y(), occ_interpt.Z()])
        else: sky +=1.0
    for direction in unitball.getDirLowerHemisphere():
        (X,Y,Z) = (direction.x,direction.y,direction.z)
        occ_interpt, occ_interface = pyliburo.py3dmodel.calculate.intersect_shape_with_ptdir(model,key,(X,Y,Z))
        if occ_interpt != None: 
            intercepts.append([occ_interpt.X(), occ_interpt.Y(), occ_interpt.Z()])
            #ground += int(int(occ_interpt.Z())==0) if ground is included in model 
        else: ground +=1.0
    SVF = (sky)/(len(unitball.getDirUpperHemisphere()));
    GVF = (ground)/(len(unitball.getDirLowerHemisphere()));
    return SVF, GVF, np.array(intercepts)

def call_values(intercepts, surfpdcoord, gridsize):
    """ Given a list of intercepts, a pdcoord of surface values, and the grid size, a list of values is returned """
    visibletemps = [surfpdcoord.val_at_coord(target,gridsize).v.mean() for target in intercepts]
    return np.array(visibletemps)
    
def calc_radiation_from_values(SurfTemp, SurfReflect, SurfEmissivity):
    """ List of values for visible surface parameters. returns long and shortwave radiative components. Assumes that lists are in order and of the same length"""
    sigma =5.67*10**(-8)    
    longwave =  sum([emissivity*sigma*temp**4/Ndir for temp, emissivity in zip(SurfTemp,SurfEmissivity)])
    shortwave =  sum([reflect/Ndir for reflect in SurfReflect])
    return longwave, shortwave

def calc_Esky_emis(Ta,RH):
    """ returns scalar of longwave radiation from the sky, that needs to be factored by SVF  """
    sigma =5.67*10**(-8)
    TaK =  Ta+273.15
    vp = RH*6.1121*np.exp((18.678-(Ta)/234.4)*(Ta)/(Ta+257.14))/1000
    skyemis = 1.24*(vp/Ta)**(1/7.)
    Esky = sigma*skyemis*(TaK**4)*(0.82-0.25*10**(-0.0945*vp))
    return Esky

def meanradtemp(Esky,Esurf, Eground,Ereflect, solarparam, SVF, GVF,  pedestrian_albedo, shadow=False):
    """ calculates Stefan-Boltzmann equation for mean radiant temperature, from different sources of radiation in the urban environment """
    sigma =5.67*10**(-8)
    Eshort =  solarparam.diffuse_frm_sky[0]*SVF/2 + solarparam.diffuse_frm_ground[0]*GVF/2 + solarparam.direct_sol[0]*solarparam.solarviewfactor[0]*shadow+ Ereflect 
    Elong = Esky*SVF/2+Esurf+Eground
    t_mrt= ((Eshort*(1-pedestrian_albedo)+Elong)/sigma)**(1/4.) - 273.15  

    results = pd.DataFrame({
    'TMRT':[t_mrt],
    'Elong':[Elong],
    'Eshort':[Eshort]
    })      
    
    return results

def all_mrt(key,compound,pdAirTemp,pdReflect,pdSurfTemp,solarparam,model_inputs,ped_properties,gridsize=1):
    """ Accepts dataframe of solar parameters, model inputs. This function calls all the previous functions in order to calculate each component needed in the Stefan-Boltzmann equation for mean radiant temperature."""
    sigma =5.67*10**(-8)

    RH = model_inputs.RH[0]
    
    try: Esky =np.mean(calc_Esky_emis(pdAirTemp.val_at_coord(key).v, RH)) #Calculation of sky irradiance if air temperature is a pdcoord
    except AttributeError: Esky = calc_Esky_emis(pdAirTemp, RH) #calculation of Esky if air temperature is a bulk value
  
    
    svf, gvf, intercepts = fourpiradiation(key, compound) #Calculate Sky view factor, ground view factor, and locations ('intercepts') on the wall at which WVF and wall temperatures need to be retrieved. 

    shadowint = check_shadow(key, compound,solarparam.solarvector[0]) #Check if the pedestrian is in a shaded area

    try: SurfTemp =call_values(intercepts, pdSurfTemp, gridsize) # Retrieve surface temperatures at the wall intercepts if surface temperature is a pdcoord
    except AttributeError: SurfTemp = [pdSurfTemp]*len(intercepts) # If not, and surface temperature is constant, get a list of the bulk surface temperature according to the number of intercepts (this is like a wall view factor)
    
    try: SurfReflect =call_values(intercepts, pdReflect, gridsize) #Same as surface temperature. 
    except AttributeError: SurfReflect = [pdReflect]*len(intercepts)

    if np.isnan(SurfTemp).any(): #If any of the intercept locations do not have a corresponding surface temperature in the input data, warn the user; treat it as sky.
        print 'Warning: ' , sum(np.isnan(SurfTemp)), ' intercepts do not have values. Treated as Sky'
        SurfTemp = SurfTemp[~np.isnan(SurfTemp)]
        SurfReflect = SurfReflect[~np.isnan(SurfReflect)]
        svf+=sum(np.isnan(SurfTemp))/Ndir
        
    SurfAlbedo, SurfEmissivity =  [[x]*len(SurfTemp) for x in [model_inputs.wall_albedo[0], model_inputs.wall_emissivity[0]]] # Repeat surface albedo and emissivity for the same number (N_intercepts) of intercepts. These could be coded differently to be treated as detailed pdcoords. 
    Elwall, Eswall = calc_radiation_from_values(SurfTemp, SurfReflect, SurfEmissivity) #Calculate radiation based on the four parameters. Surf Albedo not included since SurfReflect already accounts for albedo. 
    Eground = model_inputs.ground_emissivity[0]*sigma*gvf/2*model_inputs.groundtemp[0]**4 # Calculate radiation from the ground based on the ground temperature.
    
    mrtresults =  meanradtemp(Esky,Elwall, Eground,Eswall, solarparam, svf,gvf, ped_properties.body_albedo[0], shadow=shadowint) #calculate Tmrt according to stefan-boltzmann. 
                         
    results = pd.DataFrame({
    'TMRT':[mrtresults.TMRT[0]],
    'Elong':[mrtresults.Elong[0]],
    'Eshort':[mrtresults.Eshort[0]],
    'SVF':[svf],
#    'Elwall':[Elwall],
#    'Eswall':[Eswall],
#    'Eground':[Eground],
    'Esky':[Esky*svf/2],
    'shadow':[bool(shadowint)]
    })    
    return results

#%% SET Calculations 

def calc_SET(microclimate,ped_properties):
    """
    Parameters
    ---------
    ped_properties:  DataFrame with columns 
    ped_constants: Properties of a typical standing person. Dataframe with columns       'eff_radiation_surface_area_ratio'
    microclimate: DataFrame with columns        'air_temperature','wind_speed','mean_radiant_temperature','RH'
    
    See Gagge 1986 and the thesis that accompanies this GitHub (Sin 2017) for details on each variable. 
    """
    k=0.155; #unit conversion factor for 1clo to m^2*K/W
    pt=101.325;   #local atmosphere presssure in kPa
    sigma =5.67*10**(-8)
    
    wind_speed = abs(microclimate['wind_speed'][0])    
    dubois_area = 0.202*ped_properties['mass']**0.425*ped_properties['height']**0.725
    
    H = ped_properties['met']*58.2 - ped_properties['work']   #1 met = 58.2
 
    body_mu = ped_properties['work']/ped_properties['met']/58.2 #1 met = 58.2
    heat_produced = ped_properties['met']*(1-body_mu)*58.2 #1 met = 58.2
    Tsk = ped_properties['T_skin'] = 35.7 - 0.032*heat_produced/dubois_area #Auliciems and Szokolay pg 19
    
    Ta = microclimate['T_air'][0]
#    interm=(np.log(2)+((17.27*(Ta)/(237.2+(Ta)))))/17.27
#    dpt=(237.3*interm)/(1-interm) #dew point and dry bulb are different
    
    Tdp = (243.04*(np.log(microclimate['RH']/100.) + (17.625 *Ta)/(243.04 +Ta)))/(17.625-np.log(microclimate['RH']/100.)-17.625 *Ta/(243.04+Ta))
        
    vp = microclimate['vapor_pressure'] = microclimate['RH']/100*0.133322*np.exp(20.386-5132/(Ta+273.15)) # water vapor pressure at air temp; units in kPa antoine equation
    pssk=  np.exp(20.386-5132/(ped_properties['T_skin']+273.15))*.133322368 #water vapor pressure at skin; units in kPa
    pdpt= np.exp(20.386-5132/(Tdp+273.15))*.133322368

    #Hsk =  ped_properties['met']*58.2 - ped_properties['work'] - 0.0173*ped_properties['met']*58.2*(5.87-vp)- 0.0014*ped_properties['met']*58.2*(34-Ta)
    
    
    Icl = k*ped_properties['iclo']
    ped_properties['T_clothing']  = ped_properties['T_skin'] \
        - 0.0275*(H)\
        - Icl*((H)-3.05*(5.73-0.007*(H)-vp) \
        - 0.42*((H)-58.15) \
        - 0.0173*ped_properties['met']*58.2*(5.87-vp) \
        - 0.0014*ped_properties['met']*58.2*(34-Ta)) #other equation in Ye et al 2003, Doherty 1998
       
    
    # heat transfer coefficients and operative temperature, pressure
    if 5.66 *(ped_properties['met'][0] - 0.85)**0.39 >= 8.9*wind_speed**0.5:
        hsc = 5.66 *(ped_properties['met'] - 0.85)**0.39
    else:
        hsc = 8.6*wind_speed**0.53 # Gagge 1986, ASHRAE convective heat transfer coeff
    
    hsc = hsc*((vp+pt)/ 101.33)**0.55
    
    lr = ped_properties['Lewis_ratio'] =  15.15 *(ped_properties['T_clothing'] + 273.2)/273.2 #[K/kPa] De Dear 1996
    he=lr*hsc;  #evaporate heat transfer coefficient
    hesp=he*(101.33/(vp+pt))**0.45; ##standard evaporate heat transfer coefficient, from Gagge,1986 and ASHRAE   
    
    hr=4*ped_properties['body_emis']*sigma*ped_properties['eff_radiation_SA_ratio']*(273.15+(ped_properties['T_clothing']+microclimate['mean_radiant_temperature'])/2.)**3;   #radiative heat transfer coefficient, ASHRAE handbook
    hz=hr+hsc;
    
    # CLOTHING PARAMETERS
    Ia=1./(hz*ped_properties['fcl']);   #intrinsic insulation of the air layer, Gagge 1986
    hp=1./(Ia+Icl);   #Sensible Heat Transfer Coefficient, Gagge 1986
    hsp=hp+hr; #overall sensible heat transfer caefficient

    Rea=1/(lr*ped_properties['fcl']*hsc); # Gagge 1986
    Recl=Icl/(lr*ped_properties['icl'])  
    hep=1/(Rea+Recl);    #insensible heat transfer coefficient,Gagge 1986
        
    #Operative temperature and pressure
    v0 =0.08 #reference wind speed
    #to= Ta+ (1. - 0.73*wind_speed**0.2)*(microclimate['mean_radiant_temperature'][0]- Ta)   # ISO 7726 , qtd in Kazkaz 2012
    if wind_speed < 0.08:
        to =  (hr*(microclimate['mean_radiant_temperature'])+hsc*Ta)/(hr+hsc) # ASHRAE, no wind correction
    else: to=(hr*(microclimate['mean_radiant_temperature'])+hsc*(Ta*(wind_speed/v0)**0.5 - ped_properties['T_skin']*(wind_speed/v0-1)**0.5))/(hr+hsc) #operative temperature, Auliciems and Szokolay
    ttso=microclimate['Tso']=(hp/hsp)*to+(1-hp/hsp)*(ped_properties['T_skin'])  #standard operative temperature[C]
    ppso=(hep/hesp)*vp+(1.-hsp/hesp)*pssk; #standard operative pressure
    
    # skin wetness
    #w = ped_properties['skin_wetness'] = (Hsk - hp*(Tsk - to))/(hep*(pssk - vp))
    func = lambda st : (ttso - st + (0.088*(hesp+hsc)/hsp)*(ppso-microclimate['RH']/100*.133322368*np.exp(20.386-5132/(st+273.15))))
    try:
        s_set = microclimate['SET']= fsolve(func,0)[0]
    except NameError: s_set = np.nan
    return s_set



