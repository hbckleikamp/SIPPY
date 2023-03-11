# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 17:11:56 2023

@author: ZR48SA
"""

import pyopenms
import numpy as np
import pandas as pd






def MassGrabber(masses,        # dataframe of masses to search with scan as index
                mzMLFile,      # mzML file to search
                MSLevel="MS1", #MS1,MS2 or both
                neighbouring_scans=20,
                ppm=10):

#%%

    #parse input data
    mcols=[i for i in masses.columns if i.startswith("mass_")]
    if not len(mcols):
        masses.columns=["mass_"+str(i) for i in list(range(len(masses.columns)))]
    icols,ucols,lcols,bcols=[[i.replace("mass_",x) for i in mcols] for x in ["intensity_","u","l","b"]]
    if "scan" in masses.columns:
        id_scans=np.array(masses['scan']).astype(int)
    else:
        id_scans=np.array(masses.index).astype(int)
    
    
    masses[lcols]=masses[mcols]*(1-ppm/1000000)
    masses[ucols]=masses[mcols]*(1+ppm/1000000)
    masses[icols]=0 #placeholder
    

    print("reading mzML")
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(mzMLFile, exp)
    
    if MSLevel=="MS1":
        exp.setSpectra([s for s in exp.getSpectra() if s.getMSLevel()==1]) #load ms1 only
    if MSLevel=="MS2":
        exp.setSpectra([s for s in exp.getSpectra() if s.getMSLevel()==2]) #load ms2 only
    
    
    #1. identify ms1_scan names
    search_scans=np.array([spec.getNativeID().split("scan=")[1] for spec in exp]).astype(int)
    
    #2. find nearest neighbors
    dif=np.argsort(abs(search_scans-id_scans[:,None]),axis=1)
    nearest=dif[:,:neighbouring_scans].astype(int) #pick closest scans
    neighbours=search_scans[nearest.flatten()].reshape(-1,neighbouring_scans)
    uscans=np.unique(neighbours)
    
    #3. set ms1 spectra to only contain neighbors
    exp=[spec for ix,spec in enumerate(exp) if search_scans[ix] in uscans]
    peaks=[np.vstack(spec.get_peaks()).T for spec in exp]
    RTs=[spec.getRT() for spec in exp] 
    #add scan_level
    

    # group ms2 scans by neighbour
    masses["neighbour_scan"]=list(neighbours)
    exdf=masses.copy()
    exdf=exdf.explode("neighbour_scan").groupby("neighbour_scan")
 
    print("extracting isotopes")
    #loop version
    counter=0
    res=[]
    for n,g in exdf:
    
    
    
        p=peaks[counter]
        Mass=p[:,0]
        Intensity=p[:,1]
        counter+=1
        
        if counter%1000==0:
            
            print("fraction completed: "+str(round(counter/len(uscans),3)))
        
    
        l=np.sort(g[lcols].values.flatten())
        f=g[ucols].values.flatten()
        u=np.sort(f)
        ua=np.argsort(f)
        
        inc=0
        ims=[]
        for im,m in enumerate(Mass):
            if m>u[inc]:
                inc+=1
                if inc==len(u):
                    break
                elif (m>l[inc]) & (m<u[inc]):
                    ims.append([inc,im])
    
            
            elif (m>l[inc]) & (m<u[inc]):
                ims.append([inc,im])
    
        ims=np.array(ims)
    
        if len(ims): #now you have to map it back to a 2d array
            s=g[icols].shape
            z0=np.zeros(s)
            incs=ims[:,0]
            ints=Intensity[ims[:,1]]
            z0[np.unravel_index(ua[incs],s)]=ints
            g[icols]=z0
            res.append(g.drop(mcols+lcols+ucols,axis=1))
            
            #res.append(g[icols+["neighbour_scan"]])
    
#%%
    print("done")
    intensities=pd.concat(res).reset_index()
    return intensities.merge(pd.DataFrame(list(zip(search_scans,RTs)),columns=["neighbour_scan","neighbour_scan_RT"]),on="neighbour_scan") #add RTS
