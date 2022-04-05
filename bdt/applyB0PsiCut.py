import argparse
parser = argparse.ArgumentParser(description="")
parser.add_argument("year"      , help = "choose among:2016,2017,2018", default = '2018')
# parser.add_argument("tag"      , help = "", default = 'punzi_noTkMu')
args = parser.parse_args()
year = args.year

import os, sys
import ROOT, math
import numpy as np
import pandas, root_numpy
from ROOT import TLorentzVector
from copy import deepcopy as dc

from PhysicsTools.HeppyCore.utils.deltar import deltaR
from os import path
sys.path.append( path.dirname( path.dirname( path.abspath('../utils/utils.py') ) ) )
from utils.utils import *

# tag = args.tag

samples = [
           'data',
           'MC_LMNR', 
           'MC_JPSI', 
           'MC_PSI', 
#            'MC_BS', 
#            'MC_BSJPSIPHI', 
#            'MC_BSJPSIKST', 
#            'MC_BJPSIK',
#            'MC_HBJPSIX'
          ]

@np.vectorize
def addRejectPsi(
            mumuMass, mumuMassE, deltaB0M, deltaJpsiM, deltaPsiPM
          ):

    passSel = ( (abs(mumuMass - JPsiMass_) > 3*mumuMassE) & \
                (abs(mumuMass - PsiPMass_) > 3*mumuMassE) &  \
           (( (mumuMass < JPsiMass_) & ~( (abs(deltaB0M - deltaJpsiM) < 0.19)) ) | \
            ( (mumuMass > PsiPMass_) & ~( (abs(deltaB0M - deltaPsiPM) < 0.08)) ) | \
            ( (mumuMass > JPsiMass_) & (mumuMass < PsiPMass_) & \
              ~( (abs(deltaB0M - deltaJpsiM) < 0.09) | (abs(deltaB0M - deltaPsiPM) < 0.07)) ))) 

    return passSel


def addRejectPsi2016(
            mumuMass, mumuMassE, deltaB0M, deltaJpsiM, deltaPsiPM
          ):

    passSel = ( (abs(mumuMass - JPsiMass_) > 3*mumuMassE) & \
                (abs(mumuMass - PsiPMass_) > 3*mumuMassE) &  \
           (( (mumuMass < JPsiMass_) & ~( (abs(deltaB0M - deltaJpsiM) < 0.20) | (abs(deltaB0M - deltaPsiPM) < 0.0 )) ) | \
            ( (mumuMass > PsiPMass_) & ~( (abs(deltaB0M - deltaJpsiM) < 0.0)  | (abs(deltaB0M - deltaPsiPM) < 0.11)) ) | \
            ( (mumuMass > JPsiMass_) & (mumuMass < PsiPMass_) & ~( (abs(deltaB0M - deltaJpsiM) < 0.10) | (abs(deltaB0M - deltaPsiPM) < 0.08)) ))) 

    return passSel

# was 190
@np.vectorize
def addKeepJpsi(
            mumuMass, mumuMassE
          ):

    passSel = (abs(mumuMass - JPsiMass_) < 3*mumuMassE)
    return passSel

@np.vectorize
def addKeepPsip(
            mumuMass, mumuMassE
          ):

    passSel = (abs(mumuMass - PsiPMass_) < 3*mumuMassE)
    return passSel

nSigma_psiRej = 3.

for str_file in samples:

    input_files = []
    print (str_file)

    for i in range(11):
        ofile = '../final_ntuples/%s%s_newphi_punzi_removeTkMu_fixBkg_B0Psicut_fixPres_part%s.root'%(year, str_file,i)

        input_files = []
        input_files.append('../final_ntuples/%s%s_newphi_punzi_removeTkMu_fixBkg_%s_fixPres_part%s.root'%(args.year, str_file, args.year, i ))        

        print ('loading dataset...')
        dataset_all = pandas.DataFrame(
            root_numpy.root2array(
                input_files,
                'ntuple',
            )
        )
        print ('\t...done. n events: ', len(dataset_all))
    
    
        dataset_all['deltaB0M']   = dataset_all.tagged_mass - B0Mass_  
        dataset_all['deltaJpsiM'] = dataset_all.mumuMass - JPsiMass_   
        dataset_all['deltaPsiPM'] = dataset_all.mumuMass - PsiPMass_   
    
        dataset_all['deltaM_jpsi'] = dataset_all.deltaB0M - dataset_all.deltaJpsiM
        dataset_all['deltaM_psi']  = dataset_all.deltaB0M - dataset_all.deltaPsiPM
    
        dataset_all['passB0Psi_jpsi'] = addKeepJpsi (dataset_all.mumuMass, dataset_all.mumuMassE)
        dataset_all['passB0Psi_psip'] = addKeepPsip (dataset_all.mumuMass, dataset_all.mumuMassE)


    
        if year == '2016':
            dataset_all['passB0Psi_lmnr'] = addRejectPsi2016(dataset_all.mumuMass, dataset_all.mumuMassE, dataset_all.deltaB0M, dataset_all.deltaJpsiM, dataset_all.deltaPsiPM )
        else:
            dataset_all['passB0Psi_lmnr'] = addRejectPsi(dataset_all.mumuMass, dataset_all.mumuMassE, dataset_all.deltaB0M, dataset_all.deltaJpsiM, dataset_all.deltaPsiPM )

        dataset_all['passB0Psi_jpsi'] = dataset_all['passB0Psi_jpsi'].astype(np.int32)
        dataset_all['passB0Psi_psip'] = dataset_all['passB0Psi_psip'].astype(np.int32)
        dataset_all['passB0Psi_lmnr'] = dataset_all['passB0Psi_lmnr'].astype(np.int32)
    
        import root_pandas
        dataset_all.to_root(ofile, key='ntuple')#, store_index=False)
    