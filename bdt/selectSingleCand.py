import argparse
parser = argparse.ArgumentParser(description="")
parser.add_argument("year"      , help = "choose among:2016,2017,2018", default = '2018')
parser.add_argument("--tag"      , help = "", default = 'punzi_removeTkMu_fixBkg')
args = parser.parse_args()
year = args.year

import os, sys
import ROOT, math
import numpy as np
import pandas, root_numpy
from ROOT import TLorentzVector
from copy import deepcopy as dc

muonmass_ = 0.1056583745
kaonmass_ = 0.493677
pionmass_ = 0.139570

from PhysicsTools.HeppyCore.utils.deltar import deltaR

tag = args.tag + '_' + year
varName = 'bdt_prob'

print 'BDT tag: ', tag

if year == '2016':
#     BDTCUT = 0.965 
    BDTCUT = 0.99 ## 0.984 previous bdt  
    if 'punzi_removeTkMu_fixBkg' in tag:
        BDTCUT = 0.992 ## was 0.970  
elif year == '2017':
    if tag == 'punzi_noTkMu':
        BDTCUT = 0.97 ## wrong 0.94  
    elif 'punzi_removeTkMu_fixBkg' in tag:
        BDTCUT = 0.994 ## was 0.970  
#     BDTCUT = 0.965  ## as in AN v2
elif year == '2018':
    if tag == 'sign_yesTkMu':
        BDTCUT = 0.955  
    elif tag == 'sign_noTkMu':
        BDTCUT = 0.960  ## update
    elif tag == 'punzi_yesTkMu':
        BDTCUT = 0.955   ## update
    elif tag == 'punzi_noTkMu':
        BDTCUT = 0.975 ## was 0.970  
    elif tag == 'useMinMaxIP':
        BDTCUT = 0.980 ## was 0.970  
    elif 'punzi_removeTkMu_fixBkg' in tag:
        BDTCUT = 0.990 ## was 0.970  
#         BDTCUT = 0.992 ## was 0.970  
    varName = 'bdt_prob__' + tag.replace('no', 'remove')
        
# if year == '2016':
#     BDTCUT = 0.988  
# if year == '2017':
#     BDTCUT = 0.98  
# if year == '2018':
#     BDTCUT = 0.965  
#     BDTCUT = 0.95 ## including higher mumuMass in the optimisation procedure  



samples = [
           'data',
           'MC_JPSI', 
           'MC_LMNR', 
           'MC_PSI', 
#            'MC_BS', 
#            'MC_BSJPSIPHI', 
#            'MC_BSJPSIKST', 
#            'MC_BJPSIK',
# #            'MC_BuJpsiK', 
# #            'MC_LambdaB', 
#            'data_sameSign', 
#            'MC_HBJPSIX'
          ]

tkp_lv = TLorentzVector()
mum_lv = TLorentzVector()
@np.vectorize
def addMuTkMass(
            mumPt,  mumEta,  mumPhi,  
            tkpPt,  tkpEta,  tkpPhi
          ):
        
    if mumPt == -99:
        return -99    
    
    mum_lv.SetPtEtaPhiM(mumPt, mumEta, mumPhi, muonmass_)
    tkp_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, muonmass_)

    opt1 = (mum_lv + tkp_lv ).M()
    return opt1

import root_pandas

for str_file in samples:

    input_files = []
    print (str_file)
    for i in range(11):
        input_files = []
        print 'sample: ' , str_file, '  ', i 
        if 'data' not in str_file:
            input_files.append('sub_samples/sample_%s_%s_%s_newphi_addBDT_punzi_removeTkMu_fixBkg_%s.root'%(args.year, str_file, str(i), args.year))        
        else:
            input_files.append('sub_samples/sample_%s_%s_LMNR_%s_newphi_addBDT_punzi_removeTkMu_fixBkg_%s.root'%(args.year, str_file, str(i), args.year))  
            input_files.append('sub_samples/sample_%s_%s_Charmonium_%s_newphi_addBDT_punzi_removeTkMu_fixBkg_%s.root'%(args.year, str_file, str(i), args.year))  
        print input_files

        ofile  = '../final_ntuples/%s%s_newphi_%s_fixPres_part%s.root'%(year, str_file, tag, i)
        print ofile
        isMC = False
        if 'MC' in str_file:
           isMC = True
       
       
        print ('loading dataset...')
        dataset_all = pandas.DataFrame(
           root_numpy.root2array(
               input_files,
               'ntuple',
           )
        )
        print ('\t...done. n events: ', len(dataset_all))
            
        dataset = dataset_all[ (dataset_all.pass_preselection == 1 ) & ( dataset_all.bdt_prob > BDTCUT) ]
        #     if 'yesTkMu' in tag:
        #         dataset = dataset_all[ (dataset_all.pass_preselection == 1 ) & ( getattr(dataset_all, varName) > BDTCUT) ]
        #     elif 'noTkMu' in tag:
        #         dataset = dataset_all[ (dataset_all.pass_preselectionTkMu == 1 ) & ( getattr(dataset_all, varName) > BDTCUT) ]
        
        print ('\t...passingBDT. n events: ', len(dataset))
        
        dataset['tagged_mass']  = dataset.tagB0*dataset.bMass   + (1- dataset.tagB0)*dataset.bBarMass
        print ('added reco variables...')
            
        
        dataset['mumTrkp'] = addMuTkMass(
                                     dataset.mumPt,      dataset.mumEta,      dataset.mumPhi,  
                                     dataset.kstTrkpPt,  dataset.kstTrkpEta,  dataset.kstTrkpPhi,
                                     )
        dataset['mupTrkm'] = addMuTkMass(
                                     dataset.mupPt,      dataset.mupEta,      dataset.mupPhi,  
                                     dataset.kstTrkmPt,  dataset.kstTrkmEta,  dataset.kstTrkmPhi,
                                     )
    
        clean_trk = dataset[ ( (dataset.kstTrkpTrackerMuon  == 0) & (dataset.kstTrkmTrackerMuon == 0)) ]
        clean_dr  = clean_trk [ ((clean_trk.dR_mum_trkm > 1.E-4) & (clean_trk.dR_mup_trkp > 1.E-4)) ]
        clean_mmk = clean_dr[ (( (clean_dr.mmk1 < 5.158) | (clean_dr.mmk1 > 5.398)) & ( (clean_dr.mmk2 < 5.158) | (clean_dr.mmk2 > 5.398)))]
    
        if isMC and year == '2016':
            thedupl = clean_mmk[clean_mmk.duplicated( ['eventN', 'lumi', 'runN' ],keep=False)]
            thewin  = thedupl.loc[thedupl.groupby(['eventN', 'lumi', 'runN'])['bdt_prob'].idxmax()]
        
        else:
            thedupl = clean_mmk[clean_mmk.duplicated(['eventN', 'runN' ],keep=False)]
            thewin  = thedupl.loc[thedupl.groupby(['eventN', 'runN'])['bdt_prob'].idxmax()]
        
        listd = list(thedupl.index)
        listw = list(thewin.index)
        todel = list(set(listd)-set(listw))
        clean = clean_mmk.drop([x for x in todel])
        clean.to_root(ofile, key='ntuple')#, store_index=False)
        
        
        print ('--- Job Report: Select Best Candidate ---')
        print ('n events analysed:'             , len(dataset_all))
        print ('n events passing BDT and cuts:' , len(dataset))
        print ('n events with duplicates:'      , len(todel))
        print ('final n of events:'             , len(clean))
        print ('\n')
        
        del dataset
        del clean_trk
        del clean_dr
        del clean_mmk
        del clean
