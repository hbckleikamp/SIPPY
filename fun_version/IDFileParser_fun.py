# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 16:07:13 2023

@author: ZR48SA
"""

#%% set base path

#modules
from pathlib import Path
import os
import pandas as pd
import numpy as np




def parse_fun(iterable, targets):
    s=[]
    for p in targets:
        try:
            s.append(pd.Series([i.split(p+'=')[1].split('" ')[0] for i in iterable[1:]],name=p).str.lstrip('"').str.split('"').apply(lambda x: x[0]))
        except:
            pass
    return pd.concat(s,axis=1)


def IDFileParser(file,                     # should be single file (pepXML, idXML or mzID)
                 decoy_delimiter="decoy_", # case insensitive 
                  ):


    #### assumptions ####
    
    scores=['XCorr','hyperscore','expect','expect_score','nextscore','Posterior Probability_score','Percolator PEP'] #the parsing script will attempt to parse the following scores
    #retention time is in seconds, if not then RT could be missed when parsing from pepXML
    
    #####################
    
    #Disclaimer:
    #this script needs to be tested with nterm and cterm modifications.
    #also has only been tested with Dynamic modifications in Sequest HT
    
    #####################
    #rows that will be aimed to be retrieved:
    #scans,charge,retention time,peptide,modifications,modified_aas, proteins, target_decoy, scores

    
    print("parsing IDfile "+str(Path(file).name))
    
    
    FileType=Path(file).suffix
    if file.endswith("pep.xml"): FileType=".pepXML"
    
    once=True
    with open (file,"r") as f:
        xml_string=f.read()
    
    
    #initial parse
    if FileType==".pepXML":
        ps=[i.split("</spectrum_query>")[0] for i in xml_string.split("<spectrum_query ")]
        parse_targets=["scan","charge","retention_time_sec"]
    if FileType==".idXML":
        ps=[i.split("</PeptideIdentification>")[0] for i in xml_string.split("<PeptideIdentification")]
        parse_targets=["scan","charge","RT","sequence","protein_refs",'name="target_decoy" value']+[i+'" value' for i in scores] 
    if FileType==".mzid":
        ps=[i.split("</SpectrumIdentificationResult")[0] for i in xml_string.split("<SpectrumIdentificationResult")]
        parse_targets=["scan","chargeState",'scan start time" value',"peptideEvidence_ref","peptide_ref"]+[i+'" value' for i in scores]


    
    specdf=parse_fun(ps,parse_targets)
    
    #renaming columns
    rname={"chargeState":"charge",'scan start time" value':"RT","retention_time_sec":"RT",'name="target_decoy" value':"target_decoy"}
    [rname.update({s+'" value':"Score: "+s}) for s in scores]
    specdf=specdf.rename(columns=rname)
     


    #add peptides and proteins 

    #idXML
    if FileType==".idXML":
        
        #add proteins
        prot_df=parse_fun([i.split("</ProteinHit>")[0] for i in xml_string.split("<ProteinHit")],["id","accession"])
        prot_df.columns=["protein_refs","protein"]
        specdf=specdf.merge(prot_df,on="protein_refs")
    
        #add modifications
        smods=specdf["sequence"].str.replace("[","]").str.rsplit("]",expand=True)
        specdf["peptide"]=smods.iloc[:,0::2].fillna("").apply(lambda x: "".join(x),axis=1)
        specdf["modification_delta_mass"]=smods.iloc[:,1::2].fillna(0).apply(lambda x: [float(i) for i in x if i!=0],axis=1)
        specdf["modified_aas"]=smods.iloc[:,0::2].fillna(0).apply(lambda x: [i[-1] for i in x if i!=0][:-1],axis=1)

    #mzid
    if FileType==".mzid":
        
        #add proteins
        specdf["protein"]=specdf["peptideEvidence_ref"].apply(lambda x: x[:len(decoy_delimiter)]).replace(decoy_delimiter,"").replace(decoy_delimiter.upper(),"")+specdf["peptideEvidence_ref"].apply(lambda x: x[len(decoy_delimiter):])
        specdf["protein"]=specdf["protein"].str.split("_").apply(lambda x: "_".join(x[:-2]))
        specdf["target_decoy"]="target"
        specdf.loc[(specdf["peptideEvidence_ref"].str.startswith(decoy_delimiter)) | (specdf["peptideEvidence_ref"].str.startswith(decoy_delimiter.upper())),"target_decoy"]="decoy" 
    
        #add modifications
        ps=[i.split("</Peptide>")[0] for i in xml_string.split("<Peptide ")][1:]
        rs=[]
        for p in ps:
            id=p.split('id="')[1].split('"')[0]
            seq=p.split('<PeptideSequence>')[1].split('<')[0]
            mod_res=[i.split('"')[0] for i in p.split('residues="')[1:]]
            mod_mass=[i.split('"')[0] for i in p.split('monoisotopicMassDelta="')[1:]]
            rs.append([id,seq,mod_mass,mod_res])
        rs=pd.DataFrame(rs,columns=["peptide_ref","peptide","modification_delta_mass","modified_aas"])
        specdf=specdf.merge(rs,on="peptide_ref")
        
    #pepXML
    if FileType==".pepXML":

        r=[]
        for hits in ps[1:]:
            s=hits.split('start_scan="')[1].split('"')[0]
            splithits=hits.split("<search_hit")[1:]
            for hit in splithits:
                
                peptide=hit.split('peptide="')[1].split('"')[0]
                proteins=[i.split('"')[0] for i in hit.split('protein="')[1:]]
                
                #add modification info
                mass=[i.split('"')[0] for i in hit.split('mass="')[1:]]+["0"]
                mod_aas=[]
                if "mod_nterm_mass" in hit: mod_aas+=["N-term"] 
                if 'position="'     in hit: mod_aas+=[peptide[int(i.split('"')[0])-1] for i in hit.split('position="')[1:]]
                if "mod_cterm_mass" in hit: mod_aas+=["C-term"] 
    
                #add score info
                if once:
                    score_names=[i.split('"')[0] for i in hit.split('score name="')[1:]] #determine score columns (do once)
                    once=False
                sc=[hit.split(i+'" value="')[1].split('"')[0] for i in score_names]
                
                r.append([s,peptide,proteins,mass[1:-1],mod_aas]+sc)     
        rs=pd.DataFrame(r,columns=["scan","peptide","protein","modification_mass","modified_aas"]+["Score: "+i for i in score_names])
        specdf=specdf.merge(rs,on="scan")

    score_cols=[c for c in specdf.columns if c.startswith("Score:")]
    if "modification_delta_mass" in specdf.columns:
        specdf=specdf[["scan","charge","RT","peptide","protein","modification_delta_mass","modified_aas"]+score_cols]
    else:
        specdf=specdf[["scan","charge","RT","peptide","protein","modification_mass","modified_aas"]+score_cols]
    
    return specdf
