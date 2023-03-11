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
import numpy as np
from collections import Counter

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


def IsotopeMassSeries(file,

                
                      isotope="C13",
                      isotope_series=10):
    

    ### parse isotopic data
    element=([i for i in isotope if i.isupper()]+[i for i in isotope if i.islower()])[0]
    atom_number=int("".join([i for i in isotope if i.isdigit()]))
    isotopes=isotope_table[isotope_table["symbol"]==element]
    natural_isotope=isotopes[isotopes["Standard Isotope"]]
    labelled_isotope=isotopes[isotopes["mass_number"]==atom_number]
    isotope_mass=(labelled_isotope['Relative Atomic Mass'].values-natural_isotope['Relative Atomic Mass'].values)[0]
    #### parse input data ####
    
    if type(file)==str:
        print("Predicting Natural Isotopic Abudance "+str(Path(file).name))
        #dynamic delimiter detection
        if file.endswith('.xlsx') or file.endswith('.xls'): 
            mdf=pd.read_excel(file,engine='openpyxl') 
        else:
            with open(file,"r") as f:
                header=f.readlines()[0]
            delims=[i[0] for i in Counter([i for i in header if not i.isalnum()]).most_common()]
            for delim in delims:
                try:
                    mdf=pd.read_csv(file,sep=delim)
                    if "charge" in mdf.columns:
                        break
                except:
                    pass
    
    elif type(file)==type(pd.DataFrame()):
        mdf=file
        print("Predicting Natural Isotopic Abudance")

    if "monoisotopic_mass" in mdf.columns:
        mass_col="monoisotopic_mass"
    else:
        mass_col=[i for i in mdf.columns if "mass" in i][0]
    
    mdf[["charge",mass_col]]=mdf[["charge",mass_col]].astype(float)
    
    edf=pd.DataFrame()
    for i in range(isotope_series):
        edf["mass_isotope"+str(i)]=(mdf[mass_col]+proton_mass*mdf["charge"]+isotope_mass*i)/mdf["charge"]
    edf=edf.set_index(mdf.scan)


    return edf
    

# icols,ucols,lcols,bcols=[[i.replace("m",x) for i in mcols] for x in ["i","u","l","b"]]
# edf=edf[["scan"]+mcols].drop_duplicates().reset_index()
# edf[lcols]=edf[mcols]*(1-ppm/1000000)
# edf[ucols]=edf[mcols]*(1+ppm/1000000)
# edf[icols]=0 #placeholder