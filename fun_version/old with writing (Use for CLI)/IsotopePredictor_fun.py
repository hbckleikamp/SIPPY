# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 16:39:36 2023

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

na_table_path=    str(Path(Path(basedir).parents[0],"utils","natural_isotope_abundances.tsv"))
isotope_table=pd.read_csv(na_table_path,sep="\t",index_col=[0])
isotope_table["mass_number"]=isotope_table["mass_number"].astype(int)

#linearize isotopic window
isotopic_range=np.arange(isotope_table["delta neutrons"].min(), isotope_table["delta neutrons"].max()+1)
df=pd.DataFrame(np.vstack([np.array(list(zip([i]*len(isotopic_range),isotopic_range))) for i in isotope_table["symbol"].drop_duplicates().values]),columns=["symbol","delta neutrons"])
df["delta neutrons"]=df["delta neutrons"].astype(float).astype(int)
l_iso=isotope_table.merge(df,on=["symbol","delta neutrons"],how="outer").fillna(0)
l_iso=l_iso[['symbol','Isotopic  Composition','delta neutrons']].pivot(index="symbol",columns="delta neutrons",values="Isotopic  Composition")
l_iso[list(range(int(max(l_iso.columns)+1),int(max(l_iso.columns)+1+256-len(l_iso.columns))))]=0 #right-pad with to reach 256 sampling points (needs to be power of 2, but more padding reduces instability) 
l_iso=l_iso[l_iso.columns[l_iso.columns>=0].tolist()+l_iso.columns[l_iso.columns<0].tolist()]    #fftshift
l_iso.columns=l_iso.columns.astype(int)







def IsotopePredictor(file,                     # datafrane or file, requires element columns, any text-like or excel format is accepted
                     write_file   =True,       # should output be written to a file?
                     output_folder=None,       # which folder it should go to
                     filename=None,
                     minimum_value=10**-5 #minimum relative abundance of isotopic peak to be reported (otherwise output woudl have 256 columns). 
                     ):
    
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
                    if mdf.columns.isin(l_iso.index).sum():
                        break
                except:
                    pass
    
    elif type(file)==type(pd.DataFrame()):
        mdf=file
        print("Predicting Natural Isotopic Abudance")
    
    

    #compute fft baseline
    elements=mdf.columns[mdf.columns.isin(l_iso.index)].tolist()
    one=np.ones([len(mdf),len(l_iso.columns)])*complex(1, 0)
    for e in elements:
        one*=np.fft.fft(l_iso.loc[e,:])**mdf[e].values.reshape(-1,1)
    baseline=pd.DataFrame(np.fft.ifft(one).real,columns=[i for i in l_iso.columns])
    baseline=baseline[((baseline>minimum_value).any()[(baseline>minimum_value).any()]).index]
    baseline=baseline[baseline.columns.sort_values()]
    baseline.columns=["isotope_"+str(i) for i in baseline.columns]
    mdf=pd.concat([mdf,baseline],axis=1)


    return mdf