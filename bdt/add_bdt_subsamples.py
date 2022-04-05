import argparse
parser = argparse.ArgumentParser(description="")
parser.add_argument("year"      , help = "choose among:2016,2017,2018", default = '2018')
args = parser.parse_args()
year = args.year

import os, sys
import numpy as np

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt 

# sys.path.insert(0, os.environ['HOME'] + '/.local/lib/python2.7/site-packages')
import pandas, root_numpy

import sklearn
from sklearn.externals import joblib

import ROOT
import root_pandas

samples = [
           'data_LMNR',
           'data_Charmonium',
           'MC_LMNR', 
           'MC_JPSI', 
           'MC_PSI', 
#            'MC_BS', 
#            'MC_BSJPSIPHI', 
#            'MC_BSJPSIKST', 
#            'MC_BJPSIK',
#            'MC_HBJPSIX',    # 2018
#            'data_sameSign', 
          ]


tags = [
            '_punzi_removeTkMu_fixBkg_%s'%year,
       ] 

tag = tags[0] ### remember to update pass_preselection

for str_file in samples:
    for i in range(11):  
      
        ifile = 'sub_samples/sample_%s_%s_%s_newphi.root'%(args.year, str_file, str(i))  

        print 'adding bdt score from %s classifier.pkl'%(i)
        classifier = joblib.load('results/classifier_%s_%s.pkl' %(tag,i))
        
        feat_names = [
            'bCosAlphaBS',
            'bLBS/bLBSE',
            'kstTrkmDCABS/kstTrkmDCABSE',
            'kstTrkpDCABS/kstTrkpDCABSE',
            'bVtxCL',
            'bDCABS/bDCABSE',
            'kstTrkmMinIP2D', 
            'kstTrkpMinIP2D', 
        ]
    
        additional = [
#             'kstTrkmMinIP2D', 
#             'kstTrkpMinIP2D', 
            'bMass',
            'bBarMass',
            'bEta',
            'tagB0',
            'mumNTrkLayers',
            'mupNTrkLayers',
            'mumNPixLayers',
            'mupNPixLayers',
            'mupHighPurity',
            'mumHighPurity',
            'mumTMOneStationTight',
            'mupTMOneStationTight',
            'mupPt',
            'mumPt',
            'bPt',
            'kstPt',
            'mumuPt',
            'kstTrkpTrackerMuonArbitrated',
            'kstTrkmTrackerMuonArbitrated',
            'kstTrkpHighPurity',
            'kstTrkmHighPurity',
            'kstTrkpPt',
            'kstTrkmPt',
            'kstTrkmDCABSE',
            'kstTrkmDCABS',
            'kstTrkpDCABS',
            'kstTrkpDCABSE',
            'kstBarMass',
            'kstMass',
            'kkMass',
            'bLBS',
            'bLBSE',    
            'kstTrkmGlobalMuon',
            'kstTrkmNTrkLayers',
            'kstTrkmNPixHits',
            'kstTrkpGlobalMuon',
            'kstTrkpNTrkLayers',
            'kstTrkpNPixHits',
            'eventN',
            'mumuMass',
            'mumIsoPt_dr04',
            'mupIsoPt_dr04',
            'kstTrkmIsoPt_dr04',
            'kstTrkpIsoPt_dr04',
        ]

        if args.year != '2016':
            additional.append('charge_trig_matched')
            
        print 'loading support dataset...'
        dataset_support = pandas.DataFrame(
            root_numpy.root2array(
                ifile,
                'ntuple',
                branches=feat_names + additional,
            )
        )
        print '\t...done'
        
        print 'loading dataset...'
        dataset = pandas.DataFrame(
            root_numpy.root2array(
                ifile,
                'ntuple',
            )
        )
        print '\t...done'
        

        dataset_support['min_trkMinIP2D'] = dataset_support[['kstTrkmMinIP2D', 'kstTrkpMinIP2D']].min(axis=1) 
        dataset_support['max_trkMinIP2D'] = dataset_support[['kstTrkmMinIP2D', 'kstTrkpMinIP2D']].max(axis=1) 
        dataset['min_trkMinIP2D'] = dataset[['kstTrkmMinIP2D', 'kstTrkpMinIP2D']].min(axis=1) 
        dataset['max_trkMinIP2D'] = dataset[['kstTrkmMinIP2D', 'kstTrkpMinIP2D']].max(axis=1) 
#         feat_names.append('min_trkMinIP2D')
#         feat_names.append('max_trkMinIP2D')

        ## define isolation: # tracks with pt in a cone
        dataset_support['isopt_mum_04']  = dataset_support.mumIsoPt_dr04    /dataset_support.mumPt
        dataset_support['isopt_mup_04']  = dataset_support.mupIsoPt_dr04    /dataset_support.mupPt
        dataset_support['isopt_trkm_04'] = dataset_support.kstTrkmIsoPt_dr04/dataset_support.kstTrkmPt
        dataset_support['isopt_trkp_04'] = dataset_support.kstTrkpIsoPt_dr04/dataset_support.kstTrkpPt
        dataset_support['sum_isopt_04']  = dataset_support.isopt_mum_04 + dataset_support.isopt_mup_04 + dataset_support.isopt_trkm_04 + dataset_support.isopt_trkp_04
        
        dataset['isopt_mum_04']  = dataset.mumIsoPt_dr04    /dataset.mumPt
        dataset['isopt_mup_04']  = dataset.mupIsoPt_dr04    /dataset.mupPt
        dataset['isopt_trkm_04'] = dataset.kstTrkmIsoPt_dr04/dataset.kstTrkmPt
        dataset['isopt_trkp_04'] = dataset.kstTrkpIsoPt_dr04/dataset.kstTrkpPt
        dataset['sum_isopt_04']  = dataset.isopt_mum_04 + dataset.isopt_mup_04 + dataset.isopt_trkm_04 + dataset.isopt_trkp_04
        feat_names.append('sum_isopt_04')

        ## define tagged kstar mass
        dataset_support['kstarmass']  = dataset_support.tagB0*dataset_support.kstMass +(1- dataset_support.tagB0)*dataset_support.kstBarMass
        dataset['kstarmass']          = dataset.tagB0*dataset.kstMass +(1- dataset.tagB0)*dataset.kstBarMass
        feat_names.append('kstarmass')


        ## for 2017 and 2018, adding significance of the track wrt BS (required at trigger level)
        if args.year != '2016':
            dataset['trkpDCASign']  = abs(dataset.kstTrkpDCABS/dataset.kstTrkpDCABSE)
            dataset['trkmDCASign']  = abs(dataset.kstTrkmDCABS/dataset.kstTrkmDCABSE)

            dataset_support['trkpDCASign']  = abs(dataset_support.kstTrkpDCABS/dataset_support.kstTrkpDCABSE)
            dataset_support['trkmDCASign']  = abs(dataset_support.kstTrkmDCABS/dataset_support.kstTrkmDCABSE)
      
        
        if args.year == '2016':
            dataset['pass_preselection'] =  ( dataset.mumNTrkLayers >= 6)  & ( dataset.mupNTrkLayers >= 6 ) & \
                                            ( dataset.mumNPixLayers >= 1)  & ( dataset.mupNPixLayers >= 1 ) & \
                                            ( dataset.mumHighPurity == 1 ) & ( dataset.mupHighPurity == 1 ) & \
                                            ( dataset.mumTMOneStationTight == 1 ) & ( dataset.mupTMOneStationTight == 1 ) & \
                                            ( dataset.kstTrkmHighPurity == 1 )    & ( dataset.kstTrkpHighPurity == 1    ) & \
                                            ( dataset.kkMass > 1.035 ) & \
                                            ( dataset.kstTrkmTrackerMuon == 0 ) & ( dataset.kstTrkpTrackerMuon == 0 ) & \
                                            (~((dataset.kstTrkmGlobalMuon == 1) & ( dataset.kstTrkmNTrkLayers > 5 ) & ( dataset.kstTrkmNPixHits > 0))) & \
                                            (~((dataset.kstTrkpGlobalMuon == 1) & ( dataset.kstTrkpNTrkLayers > 5 ) & ( dataset.kstTrkpNPixHits > 0))) 
        
        else:
            dataset['pass_preselection'] = ( dataset.mumTMOneStationTight == 1 ) & ( dataset.mupTMOneStationTight == 1 ) & \
                                            ( dataset.kkMass > 1.035 ) & \
                                            ( dataset.kstTrkmTrackerMuon == 0 ) & ( dataset.kstTrkpTrackerMuon == 0 ) & \
                                            (~((dataset.kstTrkmGlobalMuon == 1) & ( dataset.kstTrkmNTrkLayers > 5 ) & ( dataset.kstTrkmNPixHits > 0))) & \
                                            (~((dataset.kstTrkpGlobalMuon == 1) & ( dataset.kstTrkpNTrkLayers > 5 ) & ( dataset.kstTrkpNPixHits > 0))) & \
                                            (( (dataset.charge_trig_matched ==  1) & (dataset.kstTrkpPt > 1.2) & (dataset.trkpDCASign > 2) ) | \
                                             ( (dataset.charge_trig_matched == -1) & (dataset.kstTrkmPt > 1.2) & (dataset.trkmDCASign > 2) ) )
                
        
        dataset['pass_preselection'] = dataset['pass_preselection'].astype(np.int32)
        ## add branch for BDT 
        print 'computing probabilities...'
        bdt_prob_array = classifier.predict_proba(dataset_support[feat_names])[:,1]
        print '\t...done'
        
        print 'adding new column to the dataset...'
        dataset['bdt_prob'] = bdt_prob_array
        print '\t...done'	
        
        # https://github.com/scikit-hep/root_pandas
        dataset.to_root( ifile.replace('.root', '_addBDT%s.root'%tag), key='ntuple', store_index=False)
