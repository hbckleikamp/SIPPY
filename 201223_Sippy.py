# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 11:27:04 2023

@author: ZR48SA
"""




#%% Workflow


### Modules ###
from pathlib import Path
import os
from inspect import getsourcefile
import matplotlib.pyplot as plt

# change directory to script directory (should work on windows and mac)
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
sippy_dir=os.getcwd()
basedir=os.getcwd()

#Import submodules
os.chdir("fun_version") #ToDo, rename to functions, and put everything together inside a simple script instead of a folder
from IDFileParser_fun import *
from EleCounter_fun import *
from IsotopePredictor_fun import *
from IsotopeMassSeries_fun import *
from MassGrabber_fun import *


# try:
#     from IPython import get_ipython
#     get_ipython().magic('clear')
# except:
#     pass
os.chdir(sippy_dir)

#%%




#mzMLFiles
#IDFiles


### Parameters ###
isotope="C13"
IDFile="E:\Sippy/Datasets/PXD024174 (Ecoli Spiking)/MSFragger/simple_mods_hit1/Run1_MockU2_EcoliR1_5_2000ng.pepXML"
mzMLFile="E:\Sippy/Datasets/PXD024174 (Ecoli Spiking)/mzML/Run1_MockU2_EcoliR1_5_2000ng.mzML"

IDFiles=["E:\Sippy/Datasets/PXD024174 (Ecoli Spiking)/MSFragger/simple_mods_hit1/Run1_MockU2_EcoliR3_10_2000ng.pepXML",
"E:\Sippy/Datasets/PXD024174 (Ecoli Spiking)/MSFragger/simple_mods_hit1/Run1_MockU2_EcoliR1_1_2000ng.pepXML",
"E:\Sippy/Datasets/PXD024174 (Ecoli Spiking)/MSFragger/simple_mods_hit1/Run1_MockU2_EcoliR1_5_2000ng.pepXML",
"E:\Sippy/Datasets/PXD024174 (Ecoli Spiking)/MSFragger/simple_mods_hit1/Run1_MockU2_EcoliR1_10_2000ng.pepXML",
"E:\Sippy/Datasets/PXD024174 (Ecoli Spiking)/MSFragger/simple_mods_hit1/Run1_MockU2_EcoliR2_1_2000ng.pepXML",
"E:\Sippy/Datasets/PXD024174 (Ecoli Spiking)/MSFragger/simple_mods_hit1/Run1_MockU2_EcoliR2_5_2000ng.pepXML",
"E:\Sippy/Datasets/PXD024174 (Ecoli Spiking)/MSFragger/simple_mods_hit1/Run1_MockU2_EcoliR2_10_2000ng.pepXML",
"E:\Sippy/Datasets/PXD024174 (Ecoli Spiking)/MSFragger/simple_mods_hit1/Run1_MockU2_EcoliR3_1_2000ng.pepXML",
"E:\Sippy/Datasets/PXD024174 (Ecoli Spiking)/MSFragger/simple_mods_hit1/Run1_MockU2_EcoliR3_5_2000ng.pepXML"]

mzMLFiles=["E:\Sippy/Datasets/PXD024174 (Ecoli Spiking)/mzML/Run1_MockU2_EcoliR3_10_2000ng.mzML",
"E:\Sippy/Datasets/PXD024174 (Ecoli Spiking)/mzML/Run1_MockU2_EcoliR1_1_2000ng.mzML",
"E:\Sippy/Datasets/PXD024174 (Ecoli Spiking)/mzML/Run1_MockU2_EcoliR1_5_2000ng.mzML",
"E:\Sippy/Datasets/PXD024174 (Ecoli Spiking)/mzML/Run1_MockU2_EcoliR1_10_2000ng.mzML",
"E:\Sippy/Datasets/PXD024174 (Ecoli Spiking)/mzML/Run1_MockU2_EcoliR2_1_2000ng.mzML",
"E:\Sippy/Datasets/PXD024174 (Ecoli Spiking)/mzML/Run1_MockU2_EcoliR2_5_2000ng.mzML",
"E:\Sippy/Datasets/PXD024174 (Ecoli Spiking)/mzML/Run1_MockU2_EcoliR2_10_2000ng.mzML",
"E:\Sippy/Datasets/PXD024174 (Ecoli Spiking)/mzML/Run1_MockU2_EcoliR3_1_2000ng.mzML",
"E:\Sippy/Datasets/PXD024174 (Ecoli Spiking)/mzML/Run1_MockU2_EcoliR3_5_2000ng.mzML"]

IDFiles=[
    "E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R1a_0_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R3a_10_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R3a_5_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R3a_1_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R3a_0-25_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R3a_0-025_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R3a_0-1_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R3a_0-01_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R3a_0_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R2a_10_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R2a_5_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R2a_1_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R2a_0-25_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R2a_0-025_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R2a_0-1_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R2a_0-01_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R2a_0_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R1a_10_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R1a_5_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R1a_1_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R1a_0-25_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R1a_0-025_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R1a_0-1_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/centroid/Run4_Ecoli268_R1a_0-01_13C_400ng.pepXML",


"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R3a_10_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R1a_0_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R1a_0-01_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R1a_0-1_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R1a_0-025_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R1a_0-25_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R1a_1_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R1a_5_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R1a_10_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R2a_0_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R2a_0-01_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R2a_0-1_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R2a_0-025_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R2a_0-25_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R2a_1_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R2a_5_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R2a_10_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R3a_0_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R3a_0-01_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R3a_0-1_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R3a_0-025_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R3a_0-25_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R3a_1_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R3a_5_13C_400ng.pepXML"]

mzMLFiles=[
           
           
  "E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R1a_0_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R1a_0-01_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R1a_0-1_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R1a_0-025_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R1a_0-25_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R1a_1_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R1a_5_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R1a_10_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R2a_0_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R2a_0-01_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R2a_0-1_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R2a_0-025_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R2a_0-25_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R2a_1_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R2a_5_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R2a_10_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R3a_0_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R3a_0-01_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R3a_0-1_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R3a_0-025_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R3a_0-25_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R3a_1_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R3a_5_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/centroid/Run4_Ecoli268_R3a_10_13C_400ng.mzML",


"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R1a_1_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R1a_5_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R1a_10_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R2a_0_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R2a_0-01_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R2a_0-1_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R2a_0-025_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R2a_0-25_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R2a_1_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R2a_5_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R2a_10_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R3a_0_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R3a_0-01_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R3a_0-1_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R3a_0-025_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R3a_0-25_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R3a_1_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R3a_5_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R3a_10_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R1a_0_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R1a_0-01_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R1a_0-1_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R1a_0-025_13C_400ng.mzML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R1a_0-25_13C_400ng.mzML"]





IDFiles=[
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R2a_0_13C_400ng.pepXML",    
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R2a_0-01_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R2a_0-025_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R2a_0-1_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R2a_0-25_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R2a_1_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R2a_5_13C_400ng.pepXML",
"E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R2a_10_13C_400ng.pepXML"]

mzMLFiles=[
'E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R2a_0_13C_400ng.mzML',
'E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R2a_0-01_13C_400ng.mzML',
'E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R2a_0-025_13C_400ng.mzML',
'E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R2a_0-1_13C_400ng.mzML',
'E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R2a_0-25_13C_400ng.mzML',
'E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R2a_1_13C_400ng.mzML',
'E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R2a_5_13C_400ng.mzML',
'E:\Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R2a_10_13C_400ng.mzML']

# IDFiles=["/Volumes/One Touch/Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R1a_10_13C_400ng.pepXML"]
# mzMLFiles=["/Volumes/One Touch/Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R1a_10_13C_400ng.mzML"]

# IDFiles=["/Volumes/One Touch/Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/MSFragger/profile/Run4_Ecoli268_R2a_0-1_13C_400ng.pepXML"]
# mzMLFiles=["/Volumes/One Touch/Sippy/Sippyv3/Datasets/PXD023693 Ecoli C1-6/mzML/profile/Run4_Ecoli268_R2a_0-1_13C_400ng.mzML"]

# IDFiles.sort()
# mzMLFiles.sort()
#%%
for ix_f,IDFile in enumerate(IDFiles):
    
    print(IDFile)
    mzMLFile=mzMLFiles[ix_f]

    
    mdf=IDFileParser(IDFile)
    mdf=EleCounter(mdf)
    mdf=IsotopePredictor(mdf)
  
    #Isotope parse
    
    def ParseIsotope(isotope):
        #load in isotopes to get "Elements", and also parse isotope information
        na_table_path=    str(Path(Path(basedir).parents[0],"utils","natural_isotope_abundances.tsv"))
        isotope_table=pd.read_csv(na_table_path,sep="\t",index_col=[0])
        isotope_table["mass_number"]=isotope_table["mass_number"].astype(int)
    
        ### parse isotopic data
        element=([i for i in isotope if i.isupper()]+[i for i in isotope if i.islower()])[0]
        atom_number=int("".join([i for i in isotope if i.isdigit()]))
        isotopes=isotope_table[isotope_table["symbol"]==element]
        natural_isotope=isotopes[isotopes["Standard Isotope"]]
        labelled_isotope=isotopes[isotopes["mass_number"]==atom_number]
        labelled_isotope_natural_abudandance=labelled_isotope["Isotopic  Composition"].values[0]
        isotope_mass=(labelled_isotope['Relative Atomic Mass'].values-natural_isotope['Relative Atomic Mass'].values)[0]
        neutron_count=math.ceil(isotope_mass/neutron_mass)
        
        
        #isotope_dict=dict(((k, eval(k)) for k in ('element', 'labelled_isotope_natural_abudandance','isotope_mass','neutron_count')))

        isotope_dict=pd.Series([element,   labelled_isotope_natural_abudandance,  isotope_mass,  neutron_count],
                         index=['element', 'labelled_isotope_natural_abudandance','isotope_mass','neutron_count'])

        return isotope_dict
    
    def IsotopeMassSeries(mdf, #dataframe, should contain the columns scan, charge and mass
                          isotope_dict,
                          isotope_series=10): #isotopes that will be included in search, (either supply an integer x, which defaults to range(x) or a list)
    
        # parse isotopic data
        neutron_count,isotope_mass=isotope_dict.get("neutron_count"),isotope_dict.get("isotope_mass")
        
        # parse input dataframe
        if "monoisotopic_mass" in mdf.columns:
            mass_col="monoisotopic_mass"
        else:
            mass_col=[i for i in mdf.columns if "mass" in i][0]
        mdf[["charge",mass_col]]=mdf[["charge",mass_col]].astype(float)
        
        #compute isotope mass series
        if type(isotope_series)==int or type(isotope_series)==float:
            isotope_series=range(int(isotope_series))
        edf=pd.DataFrame()
        for i in isotope_series:
            edf["mass_"+str(i*neutron_count)]=(mdf[mass_col]+proton_mass*mdf["charge"]+isotope_mass*i)/mdf["charge"]
        edf["scan"]=mdf["scan"]
        
        return edf
        
    
    
    isotope_dict=ParseIsotope(isotope)
    masses=IsotopeMassSeries(mdf,isotope_dict,isotope_series=10)
    
    
    
    raw_intensities=MassGrabber(masses,mzMLFile,group="max",
                                neighbouring_scans=6,
                                ppm=10)
    
    
    
    from numpy.linalg import norm
    def IsotopeFilter(its,
                           minimum_cosine_similarity=0.0, # filter on minimum cosine similarity of isotope to monoisotopic peak (between 0-1)
                           default="monoisotopic"       , # which peak to default to for cosine similarity filtering (monoisotopic, most_intense, first_most_nonzero)
                           minimum_isotopes=3           , # filter on minimum isotopes that should be left after cosine similarity filtering
                           minimum_intensity=0          , # filter on minimum total intensity of isotopic peaks
                           top_intensity_fraction=0.7   , # filter on minimum intensity fraction of most intense neighbouring scan
                           remove_gapped=False          , # remove gapped isotope patterns that miss one or more isotopes in the middle of the detected envelope 
                           combine=False                , # sum intensities of neighbours                 
                           normalize=True):               # normalize intensities
    
        icols=np.array([i for i in its.columns if i.startswith("intensity_")])
    
    
        #cosine similarity filtering
        if minimum_cosine_similarity>0:
            print("cosine similarity filtering")
            ux=its["index"].nunique()
            groups=its.groupby("index")
            gs=[]
            for n,g in groups:
         
                if n%10000==0:
        
                    print("fraction completed: "+str(round(n/ux,3)))
    
                #which peak to default to for cosine similarity filtering?
                if default=="biggest_drop":                      #biggest drop biggest drop in intensity, should work better with co-eluting peptides
                    A=g[icols[g[icols].sum().diff().argmin()-1]]
                elif default=="first_most_nonzero":              #if you pick first most nonzero you might select more background noise?
                    A=g[(g[icols]>0).sum().idxmax()] 
                elif default=="monoisotopic":                    #if you pick monoisotopic, you will detect less at higher labelling and longer peptides= 
                    A=g[icols[0]] 
                elif default=="most_intense":
                    A=g[g[icols].sum().idxmax()]                 #if you pick most intense you might select more contamination from co-eluting peptides?
                    
                    
                nA=norm(A)
                if nA:
                    cosines=[]
                    for i in icols:
                        B=g[i]
                        nB=norm(B)
                        if nB:
                            cosines.append(np.dot(A,B)/(nA*nB))
                        else:
                            cosines.append(0)
                            
                    g[icols[(np.argwhere(np.array(cosines)<minimum_cosine_similarity)).flatten().tolist()]]=0
                    gs.append(g.values)
    
            gs=np.vstack(gs)
            its=pd.DataFrame(gs,columns=its.columns)
            print("done")
    
        its=its[~((its[icols]>0).sum(axis=1)<minimum_isotopes)]                                                                                # filter on minimum number of isotopic peaks
        its=its[~(its[icols].sum(axis=1)<minimum_intensity)]                                                                                   # filter on minimum total intensity
        its=its[(its[icols].sum(axis=1) /  its[icols].sum(axis=1).groupby(its["index"],sort=False).transform(np.max))>=top_intensity_fraction] # pick most intensity neighbouring scans
        if remove_gapped: its=its[its[icols].apply(np.flatnonzero,axis=1).apply(np.diff).apply(max)==1]                                        # remove gapped isotopic envelopes
        if combine: its=pd.concat([its.groupby(["index","scan"])["neighbour_scan","neighbour_scan_RT"].agg(set).applymap(list), its.groupby(["index","scan"])[icols].sum()],axis=1)
        if normalize: its[icols]=its[icols].divide(its[icols].sum(axis=1),axis=0)
    
        return its
 

    intensities=IsotopeFilter(raw_intensities,minimum_cosine_similarity=0,remove_gapped=False) #cosine similarity filtering
    #intensities=IsotopeFilter(raw_intensities,minimum_cosine_similarity=0.6,remove_gapped=True)
    #Check if remove gapped improves situation
    
    
    ##%% Binomial fitting
    
    
    #load in precomputed binomial
    precompute_binom_path=    str(Path(Path(basedir).parents[0],"utils","precomputed_binomial_space.npy"))
    binomial_space=pd.DataFrame(np.load(precompute_binom_path))
    binomial_space.columns=["element_count","labelling_rate"]+[i for i in range(binomial_space.shape[1]-2)] #["q_"+str(i) for i in range(binomial_space.shape[1]-2)]
    binomial_space=binomial_space.set_index("element_count")
    
    
    def ComputeNewBinomialSpace(binomial_space,to_fit,isotope_dict):
        
        
        element=isotope_dict.get('element')
    
        #first add missing element counts
        missing_element_counts=np.array(list(set(to_fit[element])-set(binomial_space.index)))
        if len(missing_element_counts):
            
            counts=range(missing_element_counts.min(),missing_element_counts.max()+1)
            occurrences=[int(i.replace("q_","")) for i in binomial_space.columns if i.startswith("q_")]
            samplepoints=binomial_space.groupby(binomial_space.index).size().max()
            chances=np.linspace(binomial_space.labelling_rate.min(),binomial_space.labelling_rate.max(),samplepoints) 
            
            a=[]
            for n in counts:
                for c in chances:
                    a.append([n,c]+[binom.pmf(r, n, c) for r in occurrences]) 
                    
            a=np.vstack(a)
            a[:,2:]=a[:,2:]/np.sum(a[:,2:],axis=1).reshape(-1,1)
            a_df=pd.DataFrame(a,columns=["element_count","labelling_rate"]+["q_"+str(i) for i in range(a.shape[1]-2)]).set_index("element_count")
            binomial_space=pd.concat([binomial_space,a_df],axis=0)
            
        #then add missing isotopic columns
        missing_isotopes=list(set(binomial_space.columns)-set(to_fit.columns))
        
        if len(missing_isotopes):
            
            counts=range(binomial_space.index.min().min(),binomial_space.index.min()+1)
            occurrences=[int(i.replace("q_","")) for i in missing_isotopes if i.startswith("q_")]
            samplepoints=binomial_space.groupby(binomial_space.index).size().max()
            chances=np.linspace(binomial_space.labelling_rate.min(),binomial_space.labelling_rate.max(),samplepoints) 
            
            a=[]
            for n in counts:
                for c in chances:
                    a.append([n,c]+[binom.pmf(r, n, c) for r in occurrences]) 
                    
            a=np.vstack(a)
            a[:,2:]=a[:,2:]/np.sum(a[:,2:],axis=1).reshape(-1,1)
            a_df=pd.DataFrame(a,columns=["element_count","labelling_rate"]+["q_"+str(i) for i in range(a.shape[1]-2)]).set_index("element_count")
            binomial_space=pd.concat([binomial_space,a_df[["q_"+str(i) for i in range(a.shape[1]-2)]]],axis=0)
        
        
        chance_space=binomial_space.loc[to_fit[element].astype(int).drop_duplicates().tolist(),qcols+["labelling_rate"]]
        np.save(precompute_binom_path,binomial_space.reset_index().values) #then update old precomputed binomial file 
        return chance_space
    
    
    import pandas as pd
    import numpy as np
    from scipy.stats import binom
    pd.set_option('mode.chained_assignment', None) #remove annoying warning
    
    
    #def FitBinomial(intensities,binomial_space,isotope_dict):
    
    element,labelled_isotope_natural_abudandance=isotope_dict.get("element"),isotope_dict.get('labelled_isotope_natural_abudandance')
    
    
    ### homogenize columns
    icols=[i for i in intensities.columns if i.startswith("intensity")]
    ricols=[i for i in mdf.columns if i.startswith("theoretical_isotope")]
    theor=mdf[ricols]
    theor.columns=[int(i.replace("theoretical_isotope_","")) for i in theor.columns ]
    measured=intensities.set_index("index")[icols]
    measured.columns=[int(i.replace("intensity_","")) for i in measured.columns ]
    theor[list(set(measured.columns)-set(theor.columns))]=0 #add 0 to theorerical isotopes for each missing column
    theor=theor.loc[measured.index,measured.columns]        #make sure that theor has the same columns as measured
    mv=measured.values
    ztheor=np.where(mv==0,0,theor)                          #
    ztheor=ztheor/ztheor.sum(axis=1).reshape(-1,1)          #renormalize theoretical envelope
    
    
    
 
    
    #subtract theor
    d0=mv-ztheor  #you subtract to find the 
    
    # #What does this remainders line do? is this like a polydiv?
    # remainders=1-(1+(d0.min(axis=1)/ztheor[np.arange(len(d0)),d0.argmin(axis=1)])).reshape(-1,1) #not sure if this needs to be inverted?
    #         #most negative delta   / theoretical value at most negative value
    
    # rdf=pd.DataFrame(mv-ztheor*remainders)
    
    
    #What does this remainders line do? is this like a polydiv?
    remainders=1-(d0.min(axis=1)/ztheor[np.arange(len(d0)),d0.argmin(axis=1)]).reshape(-1,1) #not sure if this needs to be inverted?
            #most negative delta   / theoretical value at most negative value
    
    #26/01/24 Should be an addition not a subtraction?~~~
    #remainders=1+(d0.min(axis=1)/ztheor[np.arange(len(d0)),d0.argmin(axis=1)]).reshape(-1,1) #not sure if this needs to be inverted?
    
    
    rdf=pd.DataFrame(mv*remainders)    
    "i"+1
    #this needs to be normalized nrdf=rdf.divide(rdf.sum(axis=1),axis=0)
    
    # skip binomfit
    
    
    #Simple 1/2 calculation:
    nrdf=rdf.divide(rdf.sum(axis=1),axis=0)
    nrdf[element]=mdf.loc[intensities["index"],element].values
    q=(nrdf.loc[:,[0,1]]>0).all(axis=1)
    d=nrdf[q]
    rs=1/(2*d.C*d.loc[:,0]/d.loc[:,1])*remainders[q].flatten()
    
    rs.to_csv("quick_method_"+Path(IDFile).stem+".tsv",sep="\t")
    
    fig,ax=plt.subplots()
    rs.plot.hist(bins=100,range=[0, 0.1])
    
    from functools import cache
    @cache
    def factorial(n): #do with Sterling approximation

        if n:

            return n * factorial(n-1) 
        
        else:
            return 1
    
    
    
  #%%
    
    
#     #Sterlings approximation
    
#     #%%
    
#     #%%
    
#     qcols=list(rdf.columns)
#     bcols,tcols,thcols,dcols=[[x+str(i) for i in qcols] for x in ["b","t","theoretical_isotope_","fit_deviation_"]]
#     rdf[rdf>1**-5]=0
#     rdf["labelling_coefficient"]=remainders
#     rdf[element]=mdf.loc[intensities["index"],element].values
#     to_fit=rdf.loc[rdf["labelling_coefficient"]<1,qcols+[element]].fillna(0)
#     to_fit[~(to_fit>0)]=0
#     to_fit[qcols]=to_fit[qcols].divide(to_fit[qcols].sum(axis=1),axis=0)

    
 
#     #Precompute binomial space
#     try:
#         chance_space=binomial_space.loc[to_fit[element].astype(int).drop_duplicates().tolist(),qcols+["labelling_rate"]]
#     except:
#         print("precomputed binomial space does not match measured isotopes, they will now be precomputed")
        
#         chance_space=ComputeNewBinomialSpace(binomial_space)
    
    
#     print("fitting binomial")
#     to_fit[bcols]=to_fit[qcols].astype(bool)
#     to_fit=to_fit.sort_values(by=bcols+[element]) #sort on isocolumns
#     pattern_groups=to_fit[bcols].groupby(bcols,sort=False)
#     ul=len(to_fit)
    
#     u,l=0,0
#     p_patterns=[]
#     lrs=[]
#     c=0
#     for pattern,group in pattern_groups:
#         pattern_cols=np.array(qcols)[list(pattern)].tolist()
    
#         #slice 
#         l+=len(group)
#         d=to_fit.iloc[u:l][pattern_cols+[element]]
#         u+=len(group)        
        
#         if u//5000>c:
#             print("fraction completed: "+str(round(u/ul,2)))
#             c+=1
        
#         cs=chance_space.loc[d[element],pattern_cols+["labelling_rate"]]
#         cs[pattern_cols]=cs[pattern_cols].divide(cs[pattern_cols].sum(axis=1),axis=0)
#         cs=cs.dropna()
        
#         cg=cs.groupby(cs.index,sort=False)
#         dg=d.groupby(element,sort=False)[pattern_cols]
    
#         for e in d[element].drop_duplicates():
#             ps=cg.get_group(e)
#             nps=ps[pattern_cols]
            
#             vals=dg.get_group(e)
#             dif=nps.values-vals.values[:,None]
#             ed=abs(dif).sum(axis=2)
#             lr=np.nanargmin(ed,axis=1).tolist()
#             lrs.append(ps["labelling_rate"].iloc[lr])
#             p=nps.iloc[lr,:]
#             p.loc[:,"total_euclidian_distance"]=ed.min()
#             p_patterns.append(p.reset_index().set_index(vals.index))
    
#     lrdf=pd.concat(p_patterns)
#     qcols=[i for i in qcols if i in lrdf.columns]
#     lrdf["mean_euclidian_distance"]=lrdf["total_euclidian_distance"]/lrdf[qcols].notnull().sum(axis=1)
#     lrdf["labelling_rate"]=np.hstack(lrs)
    
    
#     #maybe does need dbscan here?
#     bcols,tcols,thcols,dcols=[[x+str(i) for i in qcols] for x in ["b","t","theoretical_isotope_","fit_deviation_"]]
#     lrdf=lrdf.rename(columns={i:j for i,j in zip(qcols,bcols)})
   

#     rdf=pd.concat([rdf["labelling_coefficient"],lrdf[bcols+["labelling_rate"]]],axis=1).fillna(0)
    
    
#     rdf["isotope_abundance_fraction"]=rdf["labelling_coefficient"]*rdf["labelling_rate"]+(1-rdf["labelling_coefficient"])*labelled_isotope_natural_abudandance

#     #Now it is better!!!
#     #continue later with updating the random forest classifier, but also, bring back dbscan??
#     #rdf.isotope_abundance_fraction.plot.hist()
#     rdf["test"]=rdf["labelling_coefficient"]*rdf["labelling_rate"]
#     rdf.to_csv(Path(IDFile).stem+".tsv",sep="\t")
#     rdf.test.plot.hist(bins=100)

    

# #%%
# res=["/Volumes/One Touch/Sippy/Sippyv3/Run4_Ecoli268_R2a_10_13C_400ng.tsv",
# "/Volumes/One Touch/Sippy/Sippyv3/Run4_Ecoli268_R2a_5_13C_400ng.tsv",
# "/Volumes/One Touch/Sippy/Sippyv3/Run4_Ecoli268_R2a_1_13C_400ng.tsv",
# "/Volumes/One Touch/Sippy/Sippyv3/Run4_Ecoli268_R2a_0-25_13C_400ng.tsv",
# "/Volumes/One Touch/Sippy/Sippyv3/Run4_Ecoli268_R2a_0-1_13C_400ng.tsv",
# "/Volumes/One Touch/Sippy/Sippyv3/Run4_Ecoli268_R2a_0-025_13C_400ng.tsv",
# "/Volumes/One Touch/Sippy/Sippyv3/Run4_Ecoli268_R2a_0-01_13C_400ng.tsv",
# "/Volumes/One Touch/Sippy/Sippyv3/Run4_Ecoli268_R2a_0_13C_400ng.tsv"]
# for i in res:
#     fig,ax=plt.subplots()
#     rdf=pd.read_csv(i,sep="\t")
#     rdf.isotope_abundance_fraction.plot.hist(bins=20)

    #%%
    

    
    # #%%
    
    # #convex combination
    
    # print("convex combination")
    
    # #prep dataframes
    # ix_cols=[int(i.replace("intensity_","")) for i in icols] 
    # bcols,tcols,thcols,dcols=[[x+str(i) for i in ix_cols] for x in ["b","t","theoretical_isotope_","fit_deviation_"]]
    # trc=[i for i in ricols if int(i.replace("theoretical_isotope_","")) in ix_cols]
    # bdf=binomial_space.set_index("labelling_rate",append=True).loc[lrdf[["element_count","labelling_rate"]].apply(tuple).values.tolist(),ix_cols].reset_index().set_index(lrdf.index)#.drop_duplicates()
    # bdf=bdf.rename(columns={i:j for i,j in zip(ix_cols,bcols)})
    # intensities.loc[:,tcols]=mdf.loc[intensities["index"].tolist(),thcols].values
    
    # comb=intensities.merge(bdf,left_on="index",right_index=True,how="left").fillna(0)
    # mv,tv,bv=comb[icols].values,comb[tcols].values,comb[bcols].values #measured, theoretical and binomial values
    # tv,bv=np.where(mv==0,0,tv),np.where(mv==0,0,bv)           
    # tv,bv=tv/tv.sum(axis=1).reshape(-1,1),bv/bv.sum(axis=1).reshape(-1,1) 
    # comb[tcols]=tv
    # comb[bcols]=bv
    
    # coef_step=0.01 
    # coef_range=np.arange(0,1+coef_step,coef_step) 
    # coef_range=np.hstack([np.linspace(0.001,0.009,9),coef_range]) #add more sensitive combinations here? improve sensitivity?
    
    
    # eds=np.vstack([abs(mv-((1-x)*tv+x*bv)).sum(axis=1) for x in coef_range]).T
    # rcoefs=coef_range[eds.argmin(axis=1)]
    # comb["labelled_coefficient"]=rcoefs
    # comb["total_euclidian_distance"]=eds.min(axis=1)
    # comb["mean_euclidian_distance"]=comb["total_euclidian_distance"]/(comb[icols]>0).sum(axis=1).values
    # bcomb=comb.sort_values(by=["index","mean_euclidian_distance"]).groupby(["index"],sort=False).nth(0).fillna(0)
    # bcomb["isotopes_matched"]=(bcomb[icols]>0).sum(axis=1)
    
    # fcomb=bcomb[['labelling_rate',"labelled_coefficient","isotopes_matched","total_euclidian_distance","mean_euclidian_distance"]]
    # fcomb[dcols]=(bcomb[icols]-(bcomb[bcols]*bcomb[["labelled_coefficient"]].values).values-(bcomb[tcols]*(1-bcomb[["labelled_coefficient"]].values)).values).values
    # fcomb["final_labelling_(%)"]=fcomb['labelling_rate'].values*(fcomb["labelled_coefficient"].values)+(1-fcomb["labelled_coefficient"]).values*labelled_isotope_natural_abudandance
    # fcomb["final_labelling_(%)"]=fcomb["final_labelling_(%)"]*100
    
    # #Check if this calculation is correct maybe should add all of background?
    
    # ldf=mdf.merge(fcomb,left_index=True,right_index=True,how="inner")
    # ldf["element_count"]=ldf[element]
    # ldf["total_euclidian_distance"]=abs(ldf[dcols]).sum(axis=1)
    # ldf["mean_euclidian_distance"]=ldf["total_euclidian_distance"]/ldf["isotopes_matched"].values
    
    # [ldf.pop(i) for i in ldf.columns if i.startswith("theoretical_isotope_")]
    
    
    #%% Classification with decoy
    
    # decoy_delimiter="decoy_"
    
    # #decoy clustering
    # ldf["Decoy"]=False
    # ldf.loc[ldf["protein"].astype(str).str.contains(decoy_delimiter),"Decoy"]=True
    
    # score_columns=[i for i in ldf.columns if i.startswith("Score:")]
    
    # #Random forest classification for good and bad isotopic fits. 
    
    
    # #everything in between is classified
    
    
    # #% Random forest
    # from sklearn import svm, datasets
    # from sklearn.model_selection import train_test_split
    # from sklearn.ensemble import RandomForestClassifier
    # from sklearn import metrics
    
    # training_columns=['isotopes_matched','mean_euclidian_distance']+score_columns
    
    
     
    # bad=ldf.loc[ldf["Decoy"],training_columns]
    # bad["training_group"]="bad"
    
    
    # import matplotlib.pyplot as plt
    
    
    
    # good=ldf.loc[(ldf['mean_euclidian_distance']<=ldf['mean_euclidian_distance'].quantile(0.4)) | (ldf["isotopes_matched"]>=ldf["isotopes_matched"].quantile(0.9))
    #               ,training_columns]
    
    
    # good["training_group"]="good"
    
    # training_data=pd.concat([
    #                          bad,  #Fits with high number of aligned peaks and low mean eulcidian distance are seen as good fits
    #                          good #Fits of decoy peptides are taken as "bad" isotopic fits
    #                          ]).dropna()
    # x=training_data[training_columns]
    # y=training_data["training_group"].tolist()
     
    
    # rfX_train, rfX_test, rfY_train, rfY_test = train_test_split(x, y, test_size = 0.25, random_state = 42)
     
    
    # clf=RandomForestClassifier(n_estimators=100)
    # clf.fit(rfX_train,rfY_train)
    # rfY_pred=clf.predict(rfX_test)
    
    # print("Random forest accuracy:",metrics.accuracy_score(rfY_test,rfY_pred))
    # ldf["isotope_fit_classification"]=clf.predict(ldf[training_columns].values)
    # feature_imp = pd.Series(clf.feature_importances_,index=training_columns).sort_values(ascending=False)
    # print(feature_imp)
    
    
    # ldf[dcols]=ldf[dcols].round(3)
    # ldf=ldf.reset_index(drop=True)    
    
    # ldf.to_csv(IDFile.replace(".pepXML","_SIPPY.tsv"),sep="\t")

#%%
#monod fitting, and then wrap it up!
#write to functions