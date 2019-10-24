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

muonmass_ = 0.1056583745
kaonmass_ = 0.493677
pionmass_ = 0.139570

from PhysicsTools.HeppyCore.utils.deltar import deltaR

if year == '2016':
    BDTCUT = 0.988  
if year == '2017':
    BDTCUT = 0.98  
if year == '2018':
    BDTCUT = 0.965  


samples = [
#            'MC_JPSI', 
#            'data',
           'MC_LMNR', 
#            'MC_PSIPRIME', 
#            'MC_JPSIX', 
# #            'MC_BuJpsiK', 
# #            'MC_LambdaB', 
          ]



tk1_lv = TLorentzVector()
tk2_lv = TLorentzVector()
mup_lv = TLorentzVector()
mum_lv = TLorentzVector()


# @np.vectorize
# def addVars(
#             mumPt,  mumEta,  mumPhi,  
#             mupPt,  mupEta,  mupPhi,  
#             tkpPt,  tkpEta,  tkpPhi,
#             tkmPt,  tkmEta,  tkmPhi
#           ):
#           
#           
#     if mumPt == -99:
#         return -99, -99, -99     
#     
#     mmk1 = -99.
#     mmk2 = -99.
#         
#     mum_lv.SetPtEtaPhiM(mumPt, mumEta, mumPhi, muonmass_)
#     mup_lv.SetPtEtaPhiM(mupPt, mupEta, mupPhi, muonmass_)
#             
#     if tkpPt >= tkmPt:                
#       tk1_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, kaonmass_)
#       tk2_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, kaonmass_)
#     else:                
#       tk1_lv.SetPtEtaPhiM(tkmPt, tkmEta, tkmPhi, kaonmass_)
#       tk2_lv.SetPtEtaPhiM(tkpPt, tkpEta, tkpPhi, kaonmass_)
#   
#     mmk1 = (mum_lv+mup_lv+tk1_lv).M()
#     mmk2 = (mum_lv+mup_lv+tk2_lv).M()
#     
#     return mmk1, mmk2
# 
# 
# @np.vectorize
# def addDR(
#            mumEta,      mumPhi,     kstTrkmEta,  kstTrkmPhi
#           ):
#   return deltaR(mumEta, mumPhi, kstTrkmEta, kstTrkmPhi)        





for str_file in samples:

    input_files = []
    for i in range(11):
        if 'data' not in str_file:
            input_files.append('sub_samples/sample_%s_%s_%s_addBDT.root'%(args.year, str_file, str(i)))  
        else:
            input_files.append('sub_samples/sample_%s_%s_LMNR_%s_addBDT.root'%(args.year, str_file, str(i)))  
            input_files.append('sub_samples/sample_%s_%s_Charmonium_%s_addBDT.root'%(args.year, str_file, str(i)))  

    ofile = '../final_ntuples/%s%s.root'%(year, str_file)
    isMC = False
    if 'MC' in str_file:
        isMC = True

    
    print 'loading dataset...'
    dataset_all = pandas.DataFrame(
        root_numpy.root2array(
            input_files,
            'ntuple',
        )
    )
    print '\t...done. n events: ', len(dataset_all)
    
    dataset = dataset_all[ (dataset_all.pass_preselection == 1 ) & (dataset_all.bdt_prob > BDTCUT) ]
    print '\t...passingBDT. n events: ', len(dataset)

    dataset['tagged_mass']  = dataset.tagB0*dataset.bMass   + (1- dataset.tagB0)*dataset.bBarMass
    print 'added reco variables...'
    
    
    list_to_drop = [
                           'bVtxX'               , 'bVtxY'        , 'bVtxZ',
                           'kstVtxX'             , 'kstVtxY'      , 'kstVtxZ',
                           'mumuVtxX'            , 'mumuVtxY'     , 'mumuVtxZ',
                           'mumDCAVtx'           , 'mumDCAVtxE'   , 'mumDCABS', 'mumDCABSE',
                           'mupDCAVtx'           , 'mupDCAVtxE'   , 'mupDCABS', 'mupDCABSE',
                           'mumHighPurity'       , 'mupHighPurity',
                           'mumTMOneStationLoose', 'mupTMOneStationLoose',
                           'mumCL'               , 'mupCL' ,
                           'kstTrkmNormChi2'     , 'kstTrkpNormChi2',
                           'kstTrkmFracHits'             ,'kstTrkpFracHits',
                           'kstTrkmNPixLayers'           ,'kstTrkpNPixLayers',
                           'kstTrkmNTrkLayers'           ,'kstTrkpNTrkLayers',
                           'kstTrkmMuMatch'              ,'kstTrkpMuMatch',
                           'kstTrkmHighPurity'           ,'kstTrkpHighPurity',
                           'kstTrkpGlobalMuonPromptTight','kstTrkmGlobalMuonPromptTight',
                           ]
    if isMC:     
        extension = [      
                           'genB0Mass',
                           'genKstMass',
                           'genKstPx', 'genKstPy', 'genKstPz',
                           'genKstVtxY',
                           'genPsiMass',
                           'genPsiVtxX', 'genPsiVtxY', 'genPsiVtxZ'
                           ]
        list_to_drop = list_to_drop + extension
    
#     dataset = dataset.drop(list_to_drop, axis=1)
#     print '\t...removed unnecessary branches \n'
    
    
#     dataset['mmk1'], dataset['mmk2'] = addVars(
#                                           dataset.mumPt,      dataset.mumEta,      dataset.mumPhi,  
#                                           dataset.mupPt,      dataset.mupEta,      dataset.mupPhi,  
#                                           dataset.kstTrkpPt,  dataset.kstTrkpEta,  dataset.kstTrkpPhi,
#                                           dataset.kstTrkmPt,  dataset.kstTrkmEta,  dataset.kstTrkmPhi
#                                           )
#     dataset['dR_mum_trkm'] = addDR(dataset.mumEta, dataset.mumPhi, dataset.kstTrkmEta, dataset.kstTrkmPhi)
#     dataset['dR_mup_trkp'] = addDR(dataset.mupEta, dataset.mupPhi, dataset.kstTrkpEta, dataset.kstTrkpPhi)

    clean_dr  = dataset [ ((dataset.dR_mum_trkm > 1.E-4) & (dataset.dR_mup_trkp > 1.E-4)) ]
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
    
    
    import root_pandas
    clean.to_root(ofile, key='ntuple')#, store_index=False)
    
    
    print '--- Job Report: Select Best Candidate ---'
    print 'n events analysed:'             , len(dataset)
    print 'n events passing BDT and cuts:' , len(clean_mmk)
    print 'n events with duplicates:'      , len(todel)
    print 'final n of events:'             , len(clean)
    print '\n'
    