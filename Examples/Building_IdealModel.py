# -*- coding: utf-8 -*-
"""
Created on Mon May 08 00:34:26 2017
@author: Tiffany Sin 2017

Example construction of a neighborhood of square buildings on Pyliburo ('The model'). The location of the model should be compatible with the imported data (the easiest way to do so is to recenter data coordinates to the origin, and build the model about the origin). 

NOTE: Building footprint must be square; height and streetwidth can be variable

"""
import os
import sys
import pyliburo
import ExtraFunctions

origin = ExtraFunctions.quickline((0,0,0),3)
#==============METHOD #1: From a  CSV File =====================================

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

install_and_import('easygui')

import easygui

print("Open the topology input file from the dialog box ")

testfile = easygui.fileopenbox()
#current_path = os.path.dirname("__file__")
#parent_path = os.path.abspath(os.path.join(current_path, os.pardir))
#testfile = os.path.join(parent_path,'Examples','Input_Data','Topology.csv') # This is a CSV of 0 and height values that determine the location and size of each cube. The delimiter is a whitespace.
moddict1 = ExtraFunctions.makemodel_frmcsv(testfile,0.5,delimiter_str=' ') #the function saves a dictionary of the model and some other attributes. 
#Visualization
pyliburo.py3dmodel.construct.visualise([[origin]]+[[moddict1['model']]],['BLACK','BLUE'])
#%%
#==============METHOD #2: From dimensions =====================================
[street,width,height] = [5.,3.,4.]
moddict2 = ExtraFunctions.makemodelmatrix((3,5),street,width,height)
pyliburo.py3dmodel.construct.visualise([[origin]]+[[moddict2['model']]],['BLACK','RED'])

#%%
#================ SVf calculations=====================================
pedkey = (15.,8.,0.) #location of the pedestrian
pedline = ExtraFunctions.quickline(pedkey,1.5)
print 'SVF at ', pedkey, ' is ', pyliburo.skyviewfactor.calc_SVF(pedkey,moddict2['model'])
pyliburo.py3dmodel.construct.visualise([[pedline]]+[[moddict2['model']]],['BLACK','RED'])
