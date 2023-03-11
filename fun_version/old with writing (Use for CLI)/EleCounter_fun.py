# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 16:07:13 2023

@author: ZR48SA
"""

### Modules ###
from pathlib import Path
import os
from inspect import getsourcefile
import pandas as pd
import numpy as np
from collections import Counter

import config

### Load metadata files ###

# change directory to script directory (should work on windows and mac)
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
script_dir=os.getcwd()
print(os.getcwd())
basedir=os.getcwd()

#path construction
aa_formula_path  =str(Path(Path(basedir).parents[0],"utils","AA_formulas.txt"))
unimod_table_path=str(Path(Path(basedir).parents[0],"utils","unimod_parsed.txt"))
na_table_path=    str(Path(Path(basedir).parents[0],"utils","natural_isotope_abundances.tsv"))

#read composition metadata
aa_comp=pd.read_csv(aa_formula_path,sep="\t").set_index("AA") #use dataframe to load in file with masses and compositions
aa_elements=aa_comp.columns[1:].tolist()
for i in ["","N-term","C-term"]:
    aa_comp.loc[i,:]=[0]*len(aa_comp.columns) #add zero 
    
unimod_df=pd.read_csv(unimod_table_path,sep="\t").set_index("title")
unimod_elements=unimod_df.columns[3:].tolist()
unimod_df.loc["",:]=["",0,""]+[0]*len(unimod_elements) #add zero 

element_mass=pd.read_csv(na_table_path,sep="\t")
element_mass=element_mass[element_mass["Standard Isotope"]][["symbol","Relative Atomic Mass"]].set_index("symbol")



def EleCounter(file,                     # datafrane or file,  requires a column titled "peptide", any text-like or excel format is accepted
               write_file   =True,       # should output be written to a file?
               output_folder=None,       # which folder it should go to
               filename=None,            # what should the base name be of the filename
               
               #optional arguments when using only a single peptide column as input
               description="", #how data is looked up in unimod. options: total: total modification mass, delta: delta modification mass, and name: Unimod name
               delimiters=""  #delimiters used to cleave peptide sequences, if empty it will look for delimiters itself (example: ["[","]"])
                ):

    #### parse input data ####
    if type(file)==str:
        print("Calculating Elemental Composition"+str(Path(file).name))
    
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
                    if "peptide" in mdf.columns:
                        break
                except:
                    pass
            
    elif type(file)==type(pd.DataFrame()):
        mdf=file
        print("Calculating Elemental Composition")
        

    
    #parse modifications
    mod_comp=[]
    #from total mass
    if "modification_mass" in mdf.columns: 

        mods=mdf[["peptide","modified_aas","modification_mass"]].set_index("peptide")
        mods=mods.astype(str).applymap(eval).apply(pd.Series.explode).dropna().drop_duplicates()
        mods["modification_mass"]=mods["modification_mass"].astype(float)
        mods["delta"]=mods["modification_mass"]-aa_comp.loc[mods["modified_aas"].values,"Mass"].values
        for n,row in mods.iterrows():
            mdf["modification_mass"]=mdf["modification_mass"].astype(str).str.replace(str(row["modification_mass"]),str(row["delta"]),regex=False)
        mdf["modification_delta_mass"]=mdf["modification_mass"]
    
    #parse from peptide 
    elif 'modification_delta_mass' not in mdf.columns:

        mdf["peptide"]=mdf["Peptide Sequence"]
        mod_peps=mdf.loc[~mdf["peptide"].str.isalnum(),"peptide"]
        if len(mod_peps):
            
            #try to identify delimiters and type of modification info inside peptide
            if not len(description):
                if not mod_peps.str.contains(".",regex=False).sum()==len(mod_peps): #if number of mods does not equal decimals, then it is a name
                    description="name"
                    
            if not len(delimiters):
                delimiters=[k for k,v in Counter(mdf["peptide"].sum()).items() if not k.isalnum() and k!="." and k!=","]
            for d in delimiters: mdf["peptide"]=mdf["peptide"].str.replace(d,"_")   #homogenize delims
            
            mod_peps=mdf["peptide"].str.rsplit("_",expand=True).fillna("")
            
            if not len(description):
                if (mod_peps.apply(pd.to_numeric,errors="coerce")>100).sum().sum():
                    description="total"
                else:
                    description="delta"
                    
            rs=[]
            for n,r in mod_peps.iterrows(): 
                speptide=[c.strip() for c in r if c.isupper()]
                
                if description=="name":
                    desc=[c.strip().strip("'") for c in r if not c.isupper() and len(c)>3]
                else:
                    desc=[c.strip() for c in r if "." if c and len(c)>1]
                aas=[speptide[ix][-1] for ix in range(len(desc))]
                    
                
                peptide="".join(speptide).strip()
                rs.append([n,peptide,aas,desc])
            mdf=pd.DataFrame(rs,columns=["ix","peptide","modified_aas","descriptor"],index=mod_peps.index).set_index("ix").apply(pd.Series.explode).fillna("")   
            
            if description=="name":
                mdf=mdf.rename(columns={"modified_aas":"site","descriptor":"title"})
                
                um=mdf.merge(unimod_df.reset_index(),on=["site","title"]).set_index(mdf.index) 
                mod_comp=um.groupby(um.index)[unimod_elements].sum()
                mod_names=um.groupby(um.index)["title"].apply(list)
                
            elif description=="total":
                mdf["descriptor"]=mdf["descriptor"]-aa_comp.loc[mdf["modified_aas"].tolist(),"Mass"].values
            mdf=mdf.rename(columns={"modified_aas":"site","descriptor":"modification_delta_mass"})




    #### add compositions ####

    #peptide compositions
    aas=mdf["peptide"].apply(list).astype(str).str.strip("[']").str.rsplit("', '",expand=True).fillna("")
    pep_comp=pd.DataFrame(np.sum([aa_comp.loc[aas.iloc[:,i],:][aa_elements].values for i in range(len(aas.columns))],axis=0),columns=aa_elements)  
    
    #modification compositions
    if not len(mod_comp):

        mods=mdf[['modified_aas','modification_delta_mass']]
        mods=mods.astype(str).applymap(eval).apply(pd.Series.explode).fillna("")   
        mods['modification_delta_mass']=mods['modification_delta_mass'].replace("","0").astype(float)
        
        umods=mods[['modified_aas','modification_delta_mass']].drop_duplicates()
        rs=[]
        for n,row in umods.iterrows():

            m=unimod_df[unimod_df["site"]==row["modified_aas"]]
            rs.append(m.iloc[[np.argmin(abs(m["delta_mass"]-row["modification_delta_mass"]))]][unimod_elements].reset_index())
        rs=pd.concat([umods,pd.concat(rs).set_index(umods.index)],axis=1)
        mods=mods.merge(rs,on=["modified_aas","modification_delta_mass"],how="left").set_index(mods.index)
        
        mod_comp=mods.groupby(mods.index)[unimod_elements].sum()
        mod_names=mods.groupby(mods.index)["title"].apply(list)
    
    #### finalize result ####

    if len(mod_comp):
        comp=pep_comp.add(mod_comp,fill_value=0)
    else:
        comp=pep_comp.add(mod_comp,fill_value=0)
    
    comp=comp[comp.sum()[comp.sum()>0].index]
    comp["monoisotopic_mass"]=(comp*element_mass.loc[comp.columns.tolist()].astype(float).values.T).sum(axis=1)
    if len(mod_comp): comp["unimod_name"]=mod_names
    mdf=pd.concat([mdf,comp],axis=1)

    if write_file: 
        if output_folder==None: #write in same folder as input file
            output_file=str(Path(Path(file).parents[0],Path(filename).stem+"_elcounts.tsv"))
        else:
            if not Path(output_folder).is_absolute(): output_folder=str(Path(sippy_dir,output_folder))
            if not os.exists(output_folder): os.mkdir(output_folder,parents=True)
            output_file=str(Path(output_folder,Path(filename).stem+"_elcounts.tsv"))
            
        
            
        mdf.to_csv(output_file,sep="\t")
        return output_file

    return mdf