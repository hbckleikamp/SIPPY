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

# change directory to script directory (should work on windows and mac)
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
sippy_dir=os.getcwd()

#Import submodules
os.chdir("fun_version") #ToDo, rename to functions, and put everything together inside a simple script instead of a folder
from IDFileParser_fun import *
from EleCounter_fun import *
from IsotopePredictor_fun import *
from IsotopeMassSeries_fun import *
from MassGrabber_fun import *
try:
    from IPython import get_ipython
    get_ipython().magic('clear')
except:
    pass
os.chdir(sippy_dir)






#mzMLFiles
#IDFiles


### Parameters ###
isotope="C13"
IDFile="C:/Sippy/Datasets/PXD024174 (Ecoli Spiking)/MSFragger/simple_mods_hit1/Run1_MockU2_EcoliR3_5_2000ng.pepXML"
mzMLFile="C:/Sippy/Datasets/PXD024174 (Ecoli Spiking)/mzML/Run1_MockU2_EcoliR3_5_2000ng.mzML"






mdf=IDFileParser(IDFile)
mdf=EleCounter(mdf)
mdf=IsotopePredictor(mdf)
masses=IsotopeMassSeries(mdf,isotope=isotope)
intensities=MassGrabber(masses,mzMLFile)

"i"+1
#%% 

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
isotope_mass=(labelled_isotope['Relative Atomic Mass'].values-natural_isotope['Relative Atomic Mass'].values)[0]
neutron_count=math.ceil(isotope_mass/neutron_mass)


#add unique index
elements=[i for i in mdf.columns if i in isotope_table["symbol"]]
mdf=mdf.merge(mdf[["scan"]+elements].drop_duplicates().reset_index(names="u_ix"))#create unique id for each psm and elemental composition



#homogenize columns

icols=[i for i in intensities.columns if i.startswith("intensity")]
ricols=[i for i in mdf.columns if i.startswith("theoretical_isotope")]
theor=mdf.set_index("u_ix")[ricols]
theor.columns=[int(i.replace("theoretical_isotope_","")) for i in theor.columns ]
measured=intensities[icols]
measured.columns=[int(i.replace("intensity_","")) for i in measured.columns ]





#¤¤ combine based on minimum and maximum of isotope window

#%%


mv=non_unique[icols].values
residuals=np.vstack([np.polydiv(mv[ix],b)[1] for ix,b in enumerate(baseline.values)])
residuals=residuals/residuals.sum(axis=1).reshape(-1,1)
bv=bv[:,:residuals.shape[1]]
mv=mv[:,:residuals.shape[1]]

#linear combination of residuals and baseline
coef_step=0.01
coef_range=np.arange(0,1+coef_step,coef_step)
eds=np.vstack([abs(mv-((1-x)*bv+x*residuals)).sum(axis=1) for x in coef_range])
residual_coefs=coef_range[eds.argmin(axis=0)]

#weighted sum of correlated residuals  
fw=residuals*(non_unique['pattern_total_intensity']/non_unique.groupby("u_ix")['pattern_total_intensity'].transform('sum')).values.reshape(-1,1)
fw=pd.DataFrame(fw).set_index(non_unique["u_ix"])[residual_coefs>0].groupby("u_ix").sum() #retain only correlated residuals
fw=fw.divide(fw.sum(axis=1),axis=0)

#%% Multinomial fitting
fitted_isotopes=[]
binomial_step=0.01
binomial_range=np.arange(0,1+binomial_step,binomial_step)
groups=fw.merge(non_unique[["u_ix",element]],left_index=True,right_on="u_ix")[["u_ix","C"]].set_index("C").drop_duplicates().groupby(element)
for n,g in groups:
    print(n)
    vals=fw.loc[g.u_ix,:].values #binomial fitting 
    ps=pd.DataFrame([[binom.pmf(r, n, c) for r in range(len(fw.columns))]  for c in binomial_range]) #should it be actual fititng?
    nps=ps.divide(ps.sum(axis=1),axis=0).replace(np.nan,np.inf)
    

    #find closest pattern
    dif=nps.values-vals[:,None]
    ed=abs(dif).sum(axis=2)
    cs=np.nanargmin(ed,axis=1)*binomial_step
    fitted_isotopes.append(np.array(list(zip(g.u_ix.tolist(),cs)))) #add euclidian distance

#do it coarse and schwift
ja=pd.DataFrame(np.vstack(fitted_isotopes),columns=["u_ix","isotope(%)"])



r=non_unique.merge(ja,on="u_ix",how="left")
non_unique["coefs"]=residual_coefs
r["total_isotope"]=isotope_natural_abundance*(1-residual_coefs)+r["isotope(%)"]*residual_coefs


#%%

