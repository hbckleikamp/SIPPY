# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 14:15:35 2023

@author: ZR48SA
"""


### modules ###
from pathlib import Path
import os
from inspect import getsourcefile
import pandas as pd
import math
### parse isotopic information ###

# change directory to script directory (should work on windows and mac)
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
script_dir=os.getcwd()
print(os.getcwd())
basedir=os.getcwd()




proton_mass    =1.007276466621
neutron_mass   =1.00866491588

na_table_path=    str(Path(Path(basedir).parents[0],"utils","natural_isotope_abundances.tsv"))
isotope_table=pd.read_csv(na_table_path,sep="\t",index_col=[0])
isotope_table["mass_number"]=isotope_table["mass_number"].astype(int)


def IsotopeMassSeries(mdf, #dataframe, should contain the columns scan, charge and mass
                      isotope="C13",
                      isotope_series=10):

    ### parse isotopic data
    element=([i for i in isotope if i.isupper()]+[i for i in isotope if i.islower()])[0]
    atom_number=int("".join([i for i in isotope if i.isdigit()]))
    isotopes=isotope_table[isotope_table["symbol"]==element]
    natural_isotope=isotopes[isotopes["Standard Isotope"]]
    labelled_isotope=isotopes[isotopes["mass_number"]==atom_number]
    isotope_mass=(labelled_isotope['Relative Atomic Mass'].values-natural_isotope['Relative Atomic Mass'].values)[0]
    neutron_count=math.ceil(isotope_mass/neutron_mass)
    #### parse input data ####
    
    if "monoisotopic_mass" in mdf.columns:
        mass_col="monoisotopic_mass"
    else:
        mass_col=[i for i in mdf.columns if "mass" in i][0]
    
    mdf[["charge",mass_col]]=mdf[["charge",mass_col]].astype(float)
    
    edf=pd.DataFrame()
    for i in range(isotope_series):
    
        edf["mass_"+str(i*neutron_count)]=(mdf[mass_col]+proton_mass*mdf["charge"]+isotope_mass*i)/mdf["charge"]
    
    edf["scan"]=mdf["scan"]
    #edf=edf.set_index(mdf.scan)

    return edf
    

