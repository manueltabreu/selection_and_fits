import argparse
parser = argparse.ArgumentParser(description="")
parser.add_argument("year"      , help = "choose among:2016,2017,2018", default = '2018')
args = parser.parse_args()
year = args.year

import os, sys
import ROOT, math
import numpy as np
import pandas, root_numpy
from ROOT import TLorentzVector
from copy import deepcopy as dc
pandas.options.mode.chained_assignment = None  # default='warn'

muonmass_ = 0.1056583745
kaonmass_ = 0.493677
pionmass_ = 0.139570
JPsiMass_ = 3.096916

from PhysicsTools.HeppyCore.utils.deltar import deltaR

samples = [
           'data',
           'MC_LMNR', 
           'MC_JPSI', 
           'MC_PSI', 
#            'MC_BS', 
#            'MC_BSJPSIPHI', 
#            'MC_BSJPSIKST', 
#            'MC_HBJPSIX',    # 2018

         ]

tkp_lv = TLorentzVector()
tkm_lv = TLorentzVector()
mum_lv = TLorentzVector()
mup_lv = TLorentzVector()

pion_lv = TLorentzVector()
kaon_lv = TLorentzVector()


@np.vectorize
def addDR(
            mumEta,  mumPhi,  
            tkpEta,  tkpPhi
          ):
        
    if mumEta == -99:
        return -99    

    return deltaR(mumEta, mumPhi, tkpEta, tkpPhi )

@np.vectorize
def addPsi2sMass(
            mumPt,  mumEta,  mumPhi,  
            mupPt,  mupEta,  mupPhi,  
            tkmPt,  tkmEta,  tkmPhi,
            tkpPt,  tkpEta,  tkpPhi
          ):
        
    if mumPt == -99:
        return -99    
    
    mum_lv.SetPtEtaPhiM(mumPt, mumEta, mumPhi, muonmass_)
    mup_lv.SetPtEtaPhiM(mupPt, mupEta, mupPhi, muonmass_)
    tkp_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, pionmass_)
    tkm_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, pionmass_)

    opt1 = (mum_lv + mup_lv + tkp_lv + tkm_lv).M()
    return opt1

@np.vectorize
def addpipiMass(
            tkmPt,  tkmEta,  tkmPhi,  
            tkpPt,  tkpEta,  tkpPhi
          ):
        
    if tkmPt == -99:
        return -99    
    
    mum_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, pionmass_)
    tkp_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, pionmass_)

    opt1 = (mum_lv + tkp_lv ).M()
    return opt1


@np.vectorize
def addpiKMass(
            mumPt,  mumEta,  mumPhi,  
            tkpPt,  tkpEta,  tkpPhi
          ):
        
    if mumPt == -99:
        return -99    
    
    mum_lv.SetPtEtaPhiM(mumPt, mumEta, mumPhi, pionmass_)
    tkp_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, kaonmass_)

    opt1 = (mum_lv + tkp_lv ).M()
    return opt1


@np.vectorize
def addmmpi2(
            mumPt,  mumEta,  mumPhi,  
            mupPt,  mupEta,  mupPhi,  
            tkmPt,  tkmEta,  tkmPhi,
            tkpPt,  tkpEta,  tkpPhi,
            tagB0
          ):
        
    if mumPt == -99:
        return -99    
    
    mum_lv.SetPtEtaPhiM(mumPt, mumEta, mumPhi, muonmass_)
    mup_lv.SetPtEtaPhiM(mupPt, mupEta, mupPhi, muonmass_)
    if tagB0 == 1:
      pion_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, pionmass_)
      kaon_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, kaonmass_)
    else:
      pion_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, pionmass_)
      kaon_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, kaonmass_)

    opt1 = (mum_lv + mup_lv + pion_lv ).M()
    opt2 = (mum_lv + mup_lv + kaon_lv ).M()
    return opt1, opt2


@np.vectorize
def addmmpi2Paolo(
            mumuPt,  mumuEta,  mumuPhi, mumuMass, 
            tkmPt,  tkmEta,  tkmPhi,
            tkpPt,  tkpEta,  tkpPhi,
            tagB0
          ):
        
    mum_lv.SetPtEtaPhiM(mumuPt, mumuEta, mumuPhi, mumuMass)
    if tagB0 == 1:
      pion_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, pionmass_)
      kaon_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, kaonmass_)
    else:
      pion_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, pionmass_)
      kaon_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, kaonmass_)

    opt1 = (mum_lv  + pion_lv ).M()
    opt2 = (mum_lv  + kaon_lv ).M()
    return opt1, opt2


@np.vectorize
def addkstarmass(
            tkmPt, tkmEta, tkmPhi,  
            tkpPt, tkpEta, tkpPhi,
            tagB0
          ):
        
    
    if tagB0 == 1:
      kaon_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, kaonmass_)
      pion_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, pionmass_)
    else:
      kaon_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, kaonmass_)
      pion_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, pionmass_)

    opt1 = (kaon_lv + pion_lv ).M()
    return opt1

@np.vectorize
def addbwtmass(
            mumuPt,  mumuEta,  mumuPhi, mumuMass, 
            tkmPt, tkmEta, tkmPhi,  
            tkpPt, tkpEta, tkpPhi,
            tagB0
          ):
        
    
    mum_lv.SetPtEtaPhiM(mumuPt, mumuEta, mumuPhi, mumuMass)
    if tagB0 == 1:
      kaon_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, kaonmass_)
      pion_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, pionmass_)
    else:
      kaon_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, kaonmass_)
      pion_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, pionmass_)

    opt1 = (mum_lv + kaon_lv + pion_lv ).M()
    return opt1




def addmmpiKaon (row):
   if row['tagB0'] == 1 :
      return row['bBarMass']
   else:      
      return row['bMass']
def addkst2 (row):
   if row['tagB0'] == 1 :
      return row['kstBarMass']
   else:      
      return row['kstMass']

def addpi1Pt (row):
   if row['tagB0'] == 1 :
      return row['kstTrkpPt']
   else:      
      return row['kstTrkmPt']
def addpi2Pt (row):
   if row['tagB0'] == 1 :
      return row['kstTrkmPt']
   else:      
      return row['kstTrkpPt']



@np.vectorize
def addmmkkmass(
            mumPt,  mumEta,  mumPhi,  
            mupPt,  mupEta,  mupPhi,  
            tkmPt,  tkmEta,  tkmPhi,
            tkpPt,  tkpEta,  tkpPhi          ):
        
    if mumPt == -99:
        return -99    
    
    mum_lv.SetPtEtaPhiM(mumPt, mumEta, mumPhi, muonmass_)
    mup_lv.SetPtEtaPhiM(mupPt, mupEta, mupPhi, muonmass_)
    pion_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, kaonmass_)
    kaon_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, kaonmass_)

    opt1 = (mum_lv + mup_lv + pion_lv + kaon_lv).M()
    return opt1


@np.vectorize
def addmmpipimass(
            mumPt,  mumEta,  mumPhi,  
            mupPt,  mupEta,  mupPhi,  
            tkmPt,  tkmEta,  tkmPhi,
            tkpPt,  tkpEta,  tkpPhi          ):
        
    if mumPt == -99:
        return -99    
    
    mum_lv.SetPtEtaPhiM(mumPt, mumEta, mumPhi, muonmass_)
    mup_lv.SetPtEtaPhiM(mupPt, mupEta, mupPhi, muonmass_)
    pion_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, pionmass_)
    kaon_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, pionmass_)

    opt1 = (mum_lv + mup_lv + pion_lv + kaon_lv).M()
    return opt1


@np.vectorize
def addxcut(wt_mass,  wt_kstarmass,  kaonPt, pionPt, mmpiMass, mmkMass):
        
    bool1 =  ( (5.27958 - wt_mass) - 0.3 ) / (-0.1-0.3)< (((wt_kstarmass-0.896)--0.4) / (0.6--0.4))
    bool2 = kaonPt > pionPt
    bool3 = (wt_kstarmass-0.896)>0
    bool4 = (mmpiMass > 3.2) & (mmpiMass < 3.6)
    bool5 = (mmkMass >  4.7) & (mmkMass  < 4.9)
    bool6 = ((mmkMass - 3.8) / (4.8 - 3.8)) > ((mmpiMass-3)/(3.6-3))
    
    xcut = bool1 & bool2 & bool3 & bool4 & bool5 & bool6
    return xcut


for str_file in samples:

    input_files = []
    print (str_file)

    for i in range(1):
        input_files = []
        input_files.append('../final_ntuples/%s%s_newphi_punzi_removeTkMu_fixBkg_B0Psicut_fixPres.root'%(args.year, str_file ))        
        ofile = '../final_ntuples/%s%s_newphi_punzi_removeTkMu_fixBkg_B0Psicut_fixPres_addxcutvariable.root'%(year, str_file)
#         input_files.append('../final_ntuples/%s%s_newphi_punzi_noTkMu_B0PsiFlag_part%s.root'%(args.year, str_file, i ))        
#         ofile = '../final_ntuples/%s%s_newphi_punzi_noTkMu_B0PsiFlag_addVars_part%s.root'%(year, str_file,i)

        print ('loading ds...')
        ds = pandas.DataFrame(
                root_numpy.root2array(
                    input_files,
                    'ntuple',
#                     stop = 10000
            )
        )
        
        ds['wt_mass']= ds.apply (lambda row: addmmpiKaon(row), axis=1) 

        ds['wt_kstarmass']= ds.apply (lambda row: addkst2(row), axis=1) 
        
        ds['kaonPt']   = ds.apply (lambda row: addpi1Pt(row), axis=1)
        ds['pionPt']   = ds.apply (lambda row: addpi2Pt(row), axis=1)
        

        ds['mmpiMass'], ds['mmkMass'] = addmmpi2(
                                        ds.mumPt,  ds.mumEta,  ds.mumPhi,  
                                        ds.mupPt,  ds.mupEta,  ds.mupPhi,  
                                        ds.kstTrkmPt,  ds.kstTrkmEta,  ds.kstTrkmPhi,
                                        ds.kstTrkpPt,  ds.kstTrkpEta,  ds.kstTrkpPhi,
                                        ds.tagB0
                                      )
                                      
        ds['mmkkMass']= addmmkkmass(
                                        ds.mumPt,  ds.mumEta,  ds.mumPhi,  
                                        ds.mupPt,  ds.mupEta,  ds.mupPhi,  
                                        ds.kstTrkmPt,  ds.kstTrkmEta,  ds.kstTrkmPhi,
                                        ds.kstTrkpPt,  ds.kstTrkpEta,  ds.kstTrkpPhi,
                                      )
        ds['mmpipiMass']= addmmpipimass(
                                        ds.mumPt,  ds.mumEta,  ds.mumPhi,  
                                        ds.mupPt,  ds.mupEta,  ds.mupPhi,  
                                        ds.kstTrkmPt,  ds.kstTrkmEta,  ds.kstTrkmPhi,
                                        ds.kstTrkpPt,  ds.kstTrkpEta,  ds.kstTrkpPhi,
                                      )
        ds['pipiMass']= addpipiMass(
                                        ds.kstTrkmPt,  ds.kstTrkmEta,  ds.kstTrkmPhi,
                                        ds.kstTrkpPt,  ds.kstTrkpEta,  ds.kstTrkpPhi,
                                      )
        ds['xcut'] = addxcut(ds.wt_mass, ds.wt_kstarmass, ds.kaonPt, ds.pionPt, ds.mmpiMass, ds.mmkMass )
                                      
        import root_pandas
        ds.to_root(ofile, key='ntuple')#, store_index=False)
        
